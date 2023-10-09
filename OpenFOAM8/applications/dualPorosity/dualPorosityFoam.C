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
        Info << "Time = " << runTime.timeName() << nl << endl;

        if (method == "Block-Jacobi")
        {
            volScalarField uOld(u);
            volScalarField vOld(v);
            while (pimple.loop())
            {
                // solve(fvm::Sp(beta, u) - fvm::laplacian(mu, u) == beta * vOld);
                // solve(fvm::Sp(beta, v) - fvm::laplacian(mv, v) == beta * uOld);
                {
                    #include "ABCDEqn.H"
                    solve(Deqn + Ceqn.A()*(uOld) - Ceqn.H());
                    solve(Aeqn + Beqn.A()*(vOld) - Beqn.H());
                }
                uOld = u;
                vOld = v;
                Info << endl;
            }
        }
        else if (method == "Block-Gauss-Seidel")
        {
            while (pimple.loop())
            {
                // // Au = beta*v --> Au = beta*(-Cu/D.A + D.H/D.A) =
                // solve(fvm::Sp(beta, u) - fvm::laplacian(mu, u) == beta * v); // pEqn == tau*p_fr
                // // Dv = -Cu  --> D.Av - D.H = -Cu -->  v = -Cu/D.A + D.H/D.A
                // solve(fvm::Sp(beta, v) - fvm::laplacian(mv, v) == beta * u); // p_frEqn == tau*p
                {
                    #include "ABCDEqn.H"
                    solve(Aeqn + Beqn.A()*(v) - Beqn.H());
                }
                {
                    #include "ABCDEqn.H"
                    solve(Deqn + Ceqn.A()*(u) - Ceqn.H());
                }
                Info << endl;
            }
        }
        else if (method == "S3PJ-alternate")
        {
            while (pimple.loop())
            {
                {
                    #include "ABCDEqn.H"
                    solve(Aeqn - fvm::Sp(Beqn.A() * Ceqn.A() / Deqn.A(), u) == Beqn.H() - Beqn.A() * (Deqn.H() + Ceqn.H()) / Deqn.A());
                }
                {
                    #include "ABCDEqn.H"
                    solve(Deqn - fvm::Sp(Ceqn.A() * Beqn.A() / Aeqn.A(), v) == Ceqn.H() - Ceqn.A() * (Aeqn.H() + Beqn.H()) / Aeqn.A());
                }
                Info << endl;
            }
        }
        else if (method == "S2PJ-alternate")
        {
            while (pimple.loop())
            {
                {
                    #include "ABCDEqn.H"
                    // solve(Aeqn - (Beqn.A() / Deqn.A()) * Ceqn == Beqn.H() - Beqn.A() * Deqn.H() / Deqn.A()); // u equation
                    fvScalarMatrix Leqn
                    (
                        -(Beqn.A() / Deqn.A()) * Ceqn
                    );
                    Info << fvc::domainIntegrate(Foam::magSqr(Aeqn.A()*u-Aeqn.H())) << endl;
                    Info << fvc::domainIntegrate(Foam::magSqr(Leqn.A()*u-Leqn.H())) << endl;
                    Info << fvc::domainIntegrate(Foam::magSqr(beta*v)) << endl;
                    solve(Aeqn + Beqn.A()*(v) - Beqn.H() + Leqn == Leqn.A()*u - Leqn.H()); // u equation
                }
                {
                    #include "ABCDEqn.H"
                    // solve(Deqn - (Ceqn.A() / Aeqn.A()) * Beqn == Ceqn.H() - Ceqn.A() * Aeqn.H() / Aeqn.A()); // v equation
                    fvScalarMatrix Leqn
                    (
                        -(Ceqn.A() / Aeqn.A()) * Beqn
                    );
                    Info << fvc::domainIntegrate(Foam::magSqr(Deqn.A()*v-Deqn.H())) << endl;
                    Info << fvc::domainIntegrate(Foam::magSqr(Leqn.A()*v-Leqn.H())) << endl;
                    Info << fvc::domainIntegrate(Foam::magSqr(beta*u)) << endl;
                    solve(Deqn + Ceqn.A()*(u) - Ceqn.H() + Leqn == Leqn.A()*v - Leqn.H()); // v equation
                    }
                Info << endl;
            }
        }
        else if (method == "S2PJ-v")
        {
            while (pimple.loop())
            {
                // NB: x/A.A() = x/AA = (AA^-1)x
                // Au + Bv = 0 -> AAu + Bv = AH -> u = AH/AA - Bv/AA
                // Cu + Dv = 0 -> CCu + Dv = CH 
                //             -> CC*AH/AA -CC*Bv/AA + Dv = CH
                //             -> Dv - (CC/AA*B)v = CH - CC/AA*AH
                fvScalarMatrix Aeqn(fvm::Sp( beta, u) - fvm::laplacian(mu, u));
                fvScalarMatrix Ceqn(fvm::Sp(-beta, u));
                volScalarField CCinvAA(Ceqn.A()/Aeqn.A());
                solve(
                    fvm::Sp(beta, v) - fvm::laplacian(mv, v) // Dv
                    - CCinvAA*fvm::Sp(-beta, v) // equivalent to - fvm::Sp(beta*beta/Aeqn.A(),v)
                    ==
                    Ceqn.H() - CCinvAA*Aeqn.H() // equivalent to beta*Aeqn.H()/Aeqn.A()
                );
                solve(fvm::Sp(beta, u) - fvm::laplacian(mu, u) == beta*v);
                Info << endl;
            }
        }
        else if (method == "S2PJ-u")
        {
            while (pimple.loop())
            {
                fvScalarMatrix Beqn( fvm::Sp(-beta, v) );
                fvScalarMatrix Deqn( fvm::Sp( beta, v) - fvm::laplacian(mv, v) );
                volScalarField BBinvDD(Beqn.A()/Deqn.A()); // BB*DD^-1 = BB/DD
                solve(fvm::Sp( beta, u) - fvm::laplacian(mu, u) - BBinvDD*fvm::Sp(-beta, u) == Beqn.H() - BBinvDD*Deqn.H());  // compute u
                solve(fvm::Sp( beta, v) - fvm::laplacian(mv, v) == beta*u); // compute v
                Info << endl;
                Info << fvc::domainIntegrate(Foam::magSqr(beta*u - fvc::laplacian(mu, u))) << endl;
                Info << fvc::domainIntegrate(Foam::magSqr(BBinvDD*beta*u)) << endl;
                Info << fvc::domainIntegrate(Foam::magSqr(beta*v)) << endl;
            }
        }


        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    Info << "End\n"
         << endl;
    return 0;
}

// ************************************************************************* //
