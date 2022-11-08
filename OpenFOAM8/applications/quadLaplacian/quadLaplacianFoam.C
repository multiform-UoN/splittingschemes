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

        if (method=="Block-Jacobi")
        {
            volScalarField uOld(u);
            volScalarField vOld(v);
            while ( pimple.loop() )
            {
                // Au = Bv
                // Dv = Cu
                solve(-fvm::laplacian(muu, u) ==  fvc::laplacian(muv, vOld));
                solve(-fvm::laplacian(mvv, v) ==  fvc::laplacian(mvu, uOld));
                uOld = u;
                vOld = v;
                Info << endl;
            }
        }
        else if (method=="Block-Gauss-Seidel")
        {
            while ( pimple.loop() )
            {
                solve(-fvm::laplacian(muu, u) ==  fvc::laplacian(muv, v));
                solve(-fvm::laplacian(mvv, v) ==  fvc::laplacian(mvu, u));
                Info << endl;
            }
        }
        else if (method=="S2PJ-alternate")
        {
            while ( pimple.loop() )
            {
                fvScalarMatrix Aeqn(-fvm::laplacian(muu, u) );
                fvScalarMatrix Ceqn( fvm::laplacian(mvu, u) );
                volScalarField CinvA(Ceqn.A()*(scalar(1.0)/Aeqn.A()));
                solve(-fvm::laplacian(mvv, v) - CinvA*fvm::laplacian(muv, v) == - Ceqn.H() + CinvA*Aeqn.H() );
                fvScalarMatrix Beqn( fvm::laplacian(muv, v) );
                fvScalarMatrix Deqn(-fvm::laplacian(mvv, v) );
                volScalarField BinvD(Beqn.A()*(scalar(1.0)/Deqn.A()));
                solve(-fvm::laplacian(muu, u) - BinvD*fvm::laplacian(mvu, u) == -Beqn.H() + BinvD*Deqn.H() );
                Info << endl;
            }
        }
        else if (method=="S2PJ-a")
        {
            while ( pimple.loop() )
            {
                // Au = Bv -> AAu - AH = Bv -> u = (AH + Bv)/AA
                // Dv = Cu -> Dv - CAu = -CH -> Dv - (CA/AA)*Bv = -CH + (CA/AA)*AH
                fvScalarMatrix Aeqn(-fvm::laplacian(muu, u) );
                fvScalarMatrix Ceqn( fvm::laplacian(mvu, u) );
                volScalarField CinvA(Ceqn.A()*(scalar(1.0)/Aeqn.A()));
                solve(-fvm::laplacian(mvv, v) - CinvA*fvm::laplacian(muv, v) == - Ceqn.H() + CinvA*Aeqn.H() );
                solve(-fvm::laplacian(muu, u) == fvc::laplacian(muv*v) );
                Info << endl;
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
