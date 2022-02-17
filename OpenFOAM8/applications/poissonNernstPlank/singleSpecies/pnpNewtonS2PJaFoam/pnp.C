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

        // dimensionedScalar aaa("aaa",dimensionSet(0,2,0,0,0,0,0), 1.0);

        while ( pimple.loop() ) // coupling loop
        {
            // rho *= scalar(0); //- Simply reset to zero
            // rho += z1*v;

            fvScalarMatrix Aeqn( fvm::laplacian(m, u) );
            fvScalarMatrix Ceqn( -fvm::laplacian(z2*v, u) );

            volScalarField CinvA(Ceqn.A()*(scalar(1)/Aeqn.A()));

            solve( fvm::ddt(v) - fvm::div(fvc::flux(z2*fvc::grad(u)), v) - fvm::laplacian(a, v) + fvm::Sp(CinvA*z1, v) == Ceqn.H() - CinvA*Aeqn.H() - fvc::div(fvc::flux(z2*v*fvc::grad(u))));

            solve( fvm::laplacian(m, u) == -z1*v );

            Info << "Concentration  = " << Foam::gSum(v().field()*mesh.V())/Foam::gSum(mesh.V()) << endl;
            Info << "\n" << endl;
        }

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
