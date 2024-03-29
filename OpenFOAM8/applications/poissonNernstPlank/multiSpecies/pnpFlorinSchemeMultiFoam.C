/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 \\    /   O peration     | Website:  https://openfoam.org
  \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
   \\/     M anipulation  |
-------------------------------------------------------------------------------
License
This file is part of OpenFOAM.

OpenFOAM is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
pnpFlorinSchemeFoam

Description
Transient solver for a single species transport with Nernst-Planck forcing 
and Poisson potential (electrostatic approximation).
Coupling is implemented using th L scheme (from Florin Radu notes).

Authors:
Federico Municchi, Matteo Icardi, Roberto Nuca,  Nottingham (2021)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "constrainFluxes.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define specie from dictionary
typedef dictionary specie;

int main(int argc, char *argv[])
{
  #include "setRootCaseLists.H"
  #include "createTime.H"
  #include "createMesh.H"

  // Add pimple coupling controls
  pimpleControl pimple(mesh);

  #include "createFields.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nStarting time loop\n" << endl;

  while (runTime.loop())
  {
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // --- Coupling loop
    while (pimple.loop())
    {
      //-Update charge density
      {
        rho *= scalar(0); //- Simply reset to zero

        forAll(species,sp)
        {
          dimensionedScalar Z(species[sp].lookup("Z"));
          rho += cPtrL[sp]*Z;
        }
      }

      //- Solve Poisson
      while (pimple.correctNonOrthogonal())
      {

        // V.storePrevIter();
        fvScalarMatrix VEqn
        (
          - fvm::laplacian(epsilon_,V)
          ==
          - rho
        );

        // VEqn.setReference(VRefCell, VRefValue);
        VEqn.relax();
        VEqn.solve();
        V.relax();
        V.correctBoundaryConditions();
      }

      //- Solve species
      forAll(species,sp)
      {

        //- Collect pointers and coefficients
        volScalarField&   C(cPtrL[sp]);
        const dimensionedScalar D(species[sp].lookup("D"));
        const dimensionedScalar Dphi(species[sp].lookup("Dphi"));
        const dimensionedScalar Z(species[sp].lookup("Z"));

        //- Evaluate specific Nerst-Plank flux
        const surfaceScalarField phiNP("phiNP",-fvc::flux(fvc::grad(V))*(Dphi)*Z);
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;

        {
            scalarField sumPhi
            (
                fvc::surfaceSum(mag(phi+phiNP))().primitiveField()
            );

            CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

            meanCoNum =
                0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
        }

        Info<< "Courant Number mean: " << meanCoNum
            << " max: " << CoNum << endl;

        //- Non-orthogonal correction loop
        while (pimple.correctNonOrthogonal())
        {

          C.storePrevIter();
          fvScalarMatrix CEqn
          (
              fvm::ddt(C)
            + fvm::div(phi,C,"div(phi,C)")
            + fvm::div(phiNP,C,"div(phiNP,C)")
            - fvm::laplacian(D,C,"laplacian(D,C)")
          );

          //constrainFluxes(CEqn);

          CEqn.relax();
          CEqn.solve();

          C.relax();
          C.correctBoundaryConditions();

          // //- Regularise solution
          // {
          //     //- Force solution lower bound
          //     C = (mag(C) + C)/scalar(2.0);
          //     C.correctBoundaryConditions();
          //
          //     //- Relax
          //     C.relax();
          // }
        }

      }
    }

    forAll(species,sp)
    {
      Info << "Average concentration specie " << sp << " = "
           << Foam::gSum(cPtrL[sp]().field()*mesh.V())/Foam::gSum(mesh.V())
     << endl;
    }

    runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;
  }

  Info<< "End\n" << endl;

  return 0;
}


// ************************************************************************* //
