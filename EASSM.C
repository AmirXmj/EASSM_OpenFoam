/*---------------------------------------------------------------------------*\
EASSM - Implementation of the dynamic Smagorinsky
		     SGS model.
    
Copyright Information
    Copyright (C) 1991-2009 OpenCFD Ltd.
    Copyright (C) 2010-2014 Alberto Passalacqua 
    
License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "EASSM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(EASSM, 0);
addToRunTimeSelectionTable(LESModel, EASSM, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void EASSM::updateSubGridScaleFields
(
    const volSymmTensorField& D
)
{
    nuSgs_ =0*cD(D)*sqr(delta())*sqrt(magSqr(D));
//    nuSgs_.correctBoundaryConditions();
}


volScalarField EASSM::cD(const volSymmTensorField& D) const
{
    volSymmTensorField LL = dev(filter_(sqr(U())) - (sqr(filter_(U()))));

    volSymmTensorField MM =
        sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D));

    // Locally averaging MMMM on cell faces
    volScalarField MMMM = fvc::average(magSqr(MM));

    MMMM.max(VSMALL);

    volScalarField LLMM = LL && MM;

    // Performing local average on cell faces on return
    return fvc::average(LLMM)/MMMM;
}


volScalarField EASSM::cI(const volSymmTensorField& D) const
{
    volScalarField KK = 0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));

    volScalarField mm =
        sqr(delta())*(4*sqr(mag(filter_(D))) - filter_(sqr(mag(D))));

    // Locally averaging mmmm on cell faces
    volScalarField mmmm = fvc::average(magSqr(mm));

    mmmm.max(VSMALL);

    volScalarField KKmm = KK*mm;

    // Performing local average on cell faces on return
    return fvc::average(KKmm)/mmmm;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

EASSM::EASSM
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(typeName, U, phi, transport),
    GenSGSStress(U, phi, transport),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    B_
    (
        IOobject
        (
            "B",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    c1_
     (
     IOobject
     (
     "c1",
     runTime_.timeName(),
     mesh_,
     IOobject::NO_READ,
     IOobject::AUTO_WRITE
     ),
     mesh_,
     dimensionedScalar("c1", dimless, scalar(1.0))
     ),
//    BBR_
//     (
//     IOobject
//     (
//     "BBR",
//     runTime_.timeName(),
//     mesh_,
//     IOobject::MUST_READ,
//     IOobject::AUTO_WRITE
//     ),
//     mesh_
////     symmetricTensor("BB", dimless, scalar(1.0))
//     ),
//    BB_
//     (
//     IOobject
//     (
//     "BB",
//     runTime_.timeName(),
//     mesh_,
//     IOobject::MUST_READ,
//     IOobject::AUTO_WRITE
//     ),
//     mesh_
//     symmetricTensor("BBR", dimless, scalar(1.0))
//     ),
//    BBN_
//     (
//     IOobject
//     (
//     "BBN",
//     runTime_.timeName(),
//     mesh_,
//     IOobject::MUST_READ,
//     IOobject::AUTO_WRITE
//     ),
//     mesh_
////     symmetricTensor("BBN", dimless, scalar(1.0))
//     ),
    cD_
     (
     IOobject
     (
     "cD",
     runTime_.timeName(),
     mesh_,
     IOobject::NO_READ,
     IOobject::AUTO_WRITE
     ),
     mesh_,
     dimensionedScalar("cD", dimless, scalar(1.0))
     ),


    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    updateSubGridScaleFields(dev(symm(fvc::grad(U))));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void EASSM::correct
(
    const tmp<volTensorField>& gradU
)
{
    LESModel::correct(gradU);

    const volSymmTensorField D(dev(symm(gradU)));

    k_ = cI(D)*sqr(delta())*magSqr(D);
    bound(k_,  kMin_);

// to define dynamic Smagorinsky coef
    volSymmTensorField LL = dev(filter_(sqr(U())) - (sqr(filter_(U()))));
    volSymmTensorField MM =
        sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D));

    // Locally averaging MMMM on cell faces
    volScalarField MMMM = fvc::average(magSqr(MM));

    MMMM.max(VSMALL);

    volScalarField LLMM = LL && MM;

    // Performing local average on cell faces on return
    cD_= fvc::average(LLMM)/MMMM;







//Define identity Tensor
    tensor Ei(1,0,0,0,1,0,0,0,1);
//    cD_=0.5;
//beta1_=beta_1(c1(cD(D)),beta_4(c1(cD(D))));
//    tmp<volTensorField> Bp
   c1_ =166.372847337*pow(mag(cD_),0.25)+1;
//   Info<< "B_ = " << c1_ <<  endl;

  volScalarField beta4_
          = 1.65*(-1)/(pow(2.25*c1_,2)+magSqr(skew(fvc::grad(U_))*
//                                             //to normalize vorticity tensor
                                             (8.358883747/mag(symm(fvc::grad(U_))))));
 volScalarField beta1_
             =2.25*c1_*beta4_;
//   Info<< "B_ = " << beta4_ <<  endl;
 volTensorField a_ = symm(fvc::grad(U_)) & skew(fvc::grad(U_));
 volTensorField b_=   skew(fvc::grad(U_)) & symm(fvc::grad(U_));

    B_
//            =- nuSgs_*twoSymm(fvc::grad(U()))



            =  k_*(
            (2.0/3.0)
            *symm(Ei)

       //   (2*symm(fvc::grad(U_) & inv(fvc::grad(U_)))/3

              +beta1_
              *symm(fvc::grad(U_))
//                    to normalize
                *(8.358883747/mag(symm(fvc::grad(U_))))

              +beta4_
//                       to normalize
              *(69.870937496/magSqr(symm(fvc::grad(U_)))*

              symm(a_ - b_
//                    symm(fvc::grad(U_)) ^ skew(fvc::grad(U_))
//                    - skew(fvc::grad(U_)) ^ symm(fvc::grad(U_))
               )

)
)
            ;

//    BB_=-k_*(
//                (2.0/3.0)
//                *symm(Ei)

//           //   (2*symm(fvc::grad(U_) & inv(fvc::grad(U_)))/3

//                  +beta1_
//                  *symm(fvc::grad(U_))
//    //                    to normalize
//                    *(8.358883747/mag(symm(fvc::grad(U_))))

//                  +beta4_
//    //                       to normalize
//                  *(69.870937496/magSqr(symm(fvc::grad(U_)))*

//                  symm(a_ - b_
//    //                    symm(fvc::grad(U_)) ^ skew(fvc::grad(U_))
//    //                    - skew(fvc::grad(U_)) ^ symm(fvc::grad(U_))
//                   )
//))
//              ;
//    BBN_ =- nuSgs_*twoSymm(fvc::grad(U()))
            ;
//    Info<< "B_ = " << B() <<  endl;
//    BBR_=B_;

//B_=symm(Bp);
    // Update SGS Stress Tensor  at the wall
    B_.boundaryField().updateCoeffs();
//   B().boundaryManipulate(B_.boundaryField());
    B_.correctBoundaryConditions();
//B().boundaryManipulate(B_.boundaryField());

//    label wallPatchID = mesh_.boundaryMesh().findPatchID("wall");


    updateSubGridScaleFields(D);
}

bool EASSM::read()
{
    if (GenSGSStress::read())
    {
        filter_.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


tmp<volSymmTensorField> EASSM::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
                  B_
                    - nu()*dev(twoSymm(fvc::grad(U())))
        )
    );
}

tmp<fvVectorMatrix> EASSM::divDevReff
(
    volVectorField& U
) const
{
    return
    (
        fvc::div(B_)
      - fvm::laplacian(nu(), U)


//      + fvc::laplacian(nuSgs_, U, "laplacian(nuEff,U)")
//      - fvm::laplacian(nuEff(), U)
    );
}

tmp<fvVectorMatrix> EASSM::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
//    volScalarField muEff("muEff", rho*nuEff());

//    if (couplingFactor_.value() > 0.0)
//    {
//        return
//        (
//            fvc::div(rho*B_ + couplingFactor_*rho*nuSgs_*fvc::grad(U))
//          + fvc::laplacian
//            (
//                (1.0 - couplingFactor_)*rho*nuSgs_, U, "laplacian(muEff,U)"
//            )
//          - fvm::laplacian(muEff, U)
//        );
//    }
//    else
//    {
        return
        (
            fvc::div(rho*B_)
          + fvm::laplacian(rho*nu(), U)

//          + fvc::laplacian(rho*nuSgs_, U, "laplacian(muEff,U)")
//          - fvm::laplacian(muEff, U)
        );
//    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
