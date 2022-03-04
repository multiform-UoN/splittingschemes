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

    // fvScalarMatrix Beqn( fvm::laplacian(mB, v) );
    // fvScalarMatrix Deqn( fvm::laplacian(mD, v) );

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while ( pimple.loop() )
        {
            fvScalarMatrix Aeqn( fvm::laplacian(mA, u) );
            
            fvScalarMatrix Ceqn( fvm::laplacian(mC, u) );

            volScalarField CinvA(Ceqn.A()*(scalar(1.0)/Aeqn.A()));

            solve( fvm::laplacian(mD, v) - fvm::laplacian(mB*CinvA, v) == Ceqn.H() - CinvA*Aeqn.H() );

            solve( fvm::laplacian(mA, u) == - fvc::laplacian(mB*v) );

            Info << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
