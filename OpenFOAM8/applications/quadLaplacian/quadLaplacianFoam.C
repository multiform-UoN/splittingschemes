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

    Info << "\nStarting time loop\n"
         << endl;

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        if (method == "Block-Jacobi")
        {
            volScalarField uOld(u);
            volScalarField vOld(v);
            while (pimple.loop())
            {
                // Au = Bv
                // Dv = Cu
                solve(fvm::laplacian(-muu, u) == fvc::laplacian(muv, vOld));
                solve(fvm::laplacian(-mvv, v) == fvc::laplacian(mvu, uOld));
                uOld = u;
                vOld = v;
                Info << endl;
            }
        }
        else if (method == "Block-Gauss-Seidel")
        {
            while (pimple.loop())
            {
                solve(fvm::laplacian(-mvv, v) == fvc::laplacian(mvu, u));
                solve(fvm::laplacian(-muu, u) == fvc::laplacian(muv, v));
                Info << endl;
            }
        }
        else if (method == "S2PJ-alternate")
        {
            while (pimple.loop())
            {
                {
                    #include "ABCDEqn.H"
                    solve(Deqn - (Ceqn.A() / Aeqn.A()) * Beqn == Ceqn.H() - Ceqn.A() * Aeqn.H() / Aeqn.A());
                }
                {
                    #include "ABCDEqn.H"
                    solve(Aeqn - (Beqn.A() / Deqn.A()) * Ceqn == Beqn.H() - Beqn.A() * Deqn.H() / Deqn.A());
                }
                Info << endl;
            }
        }
        else if (method == "S2PJ-alternateExpl")
        {
            while (pimple.loop())
            {
                fvScalarMatrix Aeqn(fvm::laplacian(-muu, u));
                fvScalarMatrix Ceqn(fvm::laplacian(-mvu, u));
                volScalarField CCinvAA(Ceqn.A()/Aeqn.A());
                solve(fvm::laplacian(-mvv, v) - CCinvAA * fvm::laplacian(-muv, v) == Ceqn.H() - CCinvAA * Aeqn.H());
                fvScalarMatrix Beqn(fvm::laplacian(-muv, v));
                fvScalarMatrix Deqn(fvm::laplacian(-mvv, v));
                volScalarField BBinvDD(Beqn.A()/Deqn.A());
                solve(fvm::laplacian(-muu, u) - BBinvDD * fvm::laplacian(-mvu, u) == Beqn.H() - BBinvDD * Deqn.H());
                Info << endl;
            }
        }
        else if (method == "S2PJ-v")
        {
            while (pimple.loop())
            {
                // Au = Bv -> AAu - AH = Bv -> u = (AH + Bv)/AA
                // Dv = Cu -> Dv - CAu = -CH -> Dv - (CA/AA)*Bv = -CH + (CA/AA)*AH
                fvScalarMatrix Aeqn(fvm::laplacian(-muu, u));
                fvScalarMatrix Ceqn(fvm::laplacian(-mvu, u));
                volScalarField CCinvAA(Ceqn.A()/Aeqn.A());
                solve(fvm::laplacian(-mvv, v) - CCinvAA*fvm::laplacian(-muv, v) == Ceqn.H() - CCinvAA * Aeqn.H());
                solve(fvm::laplacian(-muu, u) == fvc::laplacian(muv, v));
                Info << endl;
            }
        }
        else if (method == "S2PJ-u")
        {
            while (pimple.loop())
            {
                fvScalarMatrix Beqn(fvm::laplacian(-muv, v));
                fvScalarMatrix Deqn(fvm::laplacian(-mvv, v));
                volScalarField BBinvDD(Beqn.A()/Deqn.A());
                solve(fvm::laplacian(-muu, u) - BBinvDD*fvm::laplacian(-muv, u) == Beqn.H() - BBinvDD * Deqn.H());
                solve(fvm::laplacian(-mvv, v) == fvc::laplacian(mvu, u));
                Info << endl;
            }
        }
        else if (method == "S3PJ-alternate") // CONTROLLARE SE È TUTTO OK
        {
            while (pimple.loop())
            {
                {
                    #include "ABCDEqn.H"
                    solve(Aeqn - fvm::Sp(Beqn.A()*Ceqn.A() / Deqn.A(), u) == Beqn.H() - Beqn.A()*(Deqn.H() + Ceqn.H()) / Deqn.A());
                }
                {
                    #include "ABCDEqn.H"
                    solve(Deqn - fvm::Sp(Ceqn.A()*Beqn.A() / Aeqn.A(), v) == Ceqn.H() - Ceqn.A()*(Aeqn.H() + Beqn.H()) / Aeqn.A());
                }
                Info << endl;
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
