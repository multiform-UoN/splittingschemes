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

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Add pimple coupling controls
    pimpleControl pimple(mesh);

    #include "createFields.H"

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Coupling loop
        while ( pimple.loop() )
        {
            //-Update charge density
            rho *= scalar(0); //- Simply reset to zero
            rho += Z*C;

            //- Solve Poisson
            while ( pimple.correctNonOrthogonal() )
            {
                // V.storePrevIter();
                fvScalarMatrix VEqn( -fvm::laplacian(epsilon_, V) == -rho );

                // VEqn.setReference(VRefCell, VRefValue);
                // VEqn.relax();
                VEqn.solve();
                // V.relax();
                V.correctBoundaryConditions();
            }

            //- Evaluate specific Nerst-Plank flux
            phiNP = -fvc::flux( Dphi * Z * fvc::grad(V) );
            // const surfaceScalarField phiNP("phiNP", -fvc::flux(fvc::grad(V))*(Dphi)*Z);

            // scalar CoNum = 0.0;

            // scalar meanCoNum = 0.0;

            // scalarField sumPhi( fvc::surfaceSum( mag( phi + phiNP ) )().primitiveField() );

            // CoNum = 0.5 * gMax( sumPhi / mesh.V().field() ) * runTime.deltaTValue();

            // meanCoNum = 0.5 * ( gSum( sumPhi ) / gSum( mesh.V().field() ) ) * runTime.deltaTValue();

            // Info<< "Courant Number mean: " << meanCoNum << " max: " << CoNum << endl;

            //- Non-orthogonal correction loop
            while ( pimple.correctNonOrthogonal() )
            {
                C.storePrevIter();

                fvScalarMatrix CEqn(
                    fvm::ddt(C)
                  + fvm::div(phi, C, "div(phi,C)")
                  + fvm::div(phiNP, C, "div(phiNP,C)")
                  - fvm::laplacian(D, C, "laplacian(D,C)")
                  + fvm::Sp(L,C)
                  == 
                  L*C
                );

                // CEqn.relax();
                CEqn.solve();
                // C.relax();
                C.correctBoundaryConditions();
            }

            Info << "\n" << endl;
        }

        Info << "Concentration  = " << Foam::gSum(C().field()*mesh.V())/Foam::gSum(mesh.V()) << endl;
        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
