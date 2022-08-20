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

Description: Dual Porosity solver

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

    Info << "\nStarting time loop\n" << endl;

    Info << "METHOD: " << method << "\n" << endl;


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (method=="Block-Jacobi")
        {
            volScalarField uOld(u);
            volScalarField vOld(v);
            while ( pimple.loop() )
            {
                uOld = u;
                vOld = v;
                solve( fvm::Sp(beta, u) - fvm::laplacian(mu, u) == beta*vOld );
                solve( fvm::Sp(beta, v) - fvm::laplacian(mv, v) == beta*uOld );
                Info << endl;
            }
        }
        else if (method=="Block-Gauss-Seidel")
        {
            while ( pimple.loop() )
            {
                // Au = beta*v --> Au = beta*(-Cu/D.A + D.H/D.A) = 
                solve( fvm::Sp(beta, u) - fvm::laplacian(mu, u) == beta*v ); // pEqn == tau*p_fr
                // Dv = -Cu  --> D.Av - D.H = -Cu -->  v = -Cu/D.A + D.H/D.A
                solve( fvm::Sp(beta, v) - fvm::laplacian(mv, v) == beta*u ); // p_frEqn == tau*p
                Info << endl;
            }
        }
        else if (method=="S3PJ-alternate")
        {
            while ( pimple.loop() )
            {
                {
                    #include "ABCDEqn.H"
                    solve(
                        Aeqn - fvm::Sp(Beqn.A()*Ceqn.A()/Deqn.A(), u)
                        ==
                        Beqn.H() - Beqn.A()*(Deqn.H()+Ceqn.H())/Deqn.A()
                        );
                }
                {
                    #include "ABCDEqn.H"
                    solve(
                        Deqn - fvm::Sp(Ceqn.A()*Beqn.A()/Aeqn.A(), v)
                        ==
                        Ceqn.H() - Ceqn.A()*(Aeqn.H()+Beqn.H())/Aeqn.A()
                        );
                }
                Info << endl;
            }
        }
        // TO BE DONE, AS IT'S WRITTEN NOW IT'S WRONG BUT PERHAPS WE CAN USE HERE H1() TO FIX IT
        else if (method=="S2PJ-alternate")
        {
            while ( pimple.loop() )
            {
                // #include "ABCDEqn.H"
                // solve( Aeqn - Beqn*Ceqn.A()/Deqn.A() == Beqn.H() - Beqn.A()*Deqn.H()/Deqn.A() );
                // #include "ABCDEqn.H"
                // solve( Deqn - Ceqn*Beqn.A()/Aeqn.A() == Ceqn.H() - Ceqn.A()*Aeqn.H()/Aeqn.A() );
                // Info << endl;
            }
        }
        else if (method=="S2PJ-a")
        {
            while ( pimple.loop() )
            {
                fvScalarMatrix Aeqn( fvm::Sp( beta, u) - fvm::laplacian(mu, u) );
                fvScalarMatrix Ceqn( fvm::Sp(-beta, u) );
                volScalarField CinvA(Ceqn.A()*(scalar(1.0)/Aeqn.A()));
                solve(
                    fvm::Sp(beta, v) - fvm::laplacian(mv, v) - fvm::Sp(-beta*CinvA, v)
                    ==
                    Ceqn.H() - CinvA*Aeqn.H()
                    );
                solve(
                    fvm::Sp(beta, u) - fvm::laplacian(mu, u)
                    ==
                    beta*v
                    );
                Info << endl;
            }
        }
        else if (method=="S2PJ-b")
        {
            while ( pimple.loop() )
            {
                fvScalarMatrix Aeqn( fvm::Sp(-beta, u) );
                fvScalarMatrix Ceqn( fvm::Sp( beta, u) - fvm::laplacian(mu, u) );
                volScalarField CinvA(Ceqn.A()*(scalar(1.0)/Aeqn.A()));
                solve(
                    fvm::Sp(-beta, v) - fvm::Sp(mv*CinvA, v) + fvm::laplacian(mv*CinvA, v) 
                    == 
                    Ceqn.H() - CinvA*Aeqn.H()
                    );
                solve( fvm::Sp(-beta, u) == - (beta*v - fvc::laplacian(mv*v)) );
                Info << endl;
            }
        }

        runTime.write();

        Info 
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" 
        << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
