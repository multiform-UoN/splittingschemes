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

Application: chFoam

Description: Cahn-Hilliard solver

Authors: Roberto Nuca, Nottingham (2021)
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

    pimpleControl pimple(mesh);

    #include "createFields.H"

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        // dimensionedScalar L("rho", dimensionSet(0,0,1,0,0,0,0), scalar(0));

        // fvScalarMatrix Aeqn( fvm::laplacian(mA, u) );
        // fvScalarMatrix Beqn( fvm::laplacian(mB, v) );
        // fvScalarMatrix Ceqn( fvm::laplacian(mC, u) );
        // fvScalarMatrix Deqn( fvm::laplacian(mD, v) );

        // volScalarField CinvA(Ceqn.A()*(scalar(1)/Aeqn.A()));


        while ( pimple.loop() )
        {
            // solve(fvm::laplacian(mB, v) == - fvc::laplacian(mA, u));
            // solve(fvm::laplacian(mC, u) == - fvc::laplacian(mD, u));


            solve(fvm::laplacian(mA, u) == - fvc::laplacian(mB, v));
            solve(fvm::laplacian(mD, v) == - fvc::laplacian(mC, u));            

            // V.correctBoundaryConditions();
            // V.storePrevIter();


            // Info << "Concentration  = " << Foam::gSum(U().field()*mesh.V())/Foam::gSum(mesh.V()) << endl;
            // Info << "\n" << endl;
        }

        // Info << "Concentration  = " << Foam::gSum(C().field()*mesh.V())/Foam::gSum(mesh.V()) << endl;
        runTime.write();
        // Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
