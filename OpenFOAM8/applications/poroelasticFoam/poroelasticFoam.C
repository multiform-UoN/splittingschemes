/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    poroelasticFoam

Description
    Solves the equations of poroelasticity (coupled Darcy and elasticity).
    A total lagrangian formulation is employed for the stress.

Contributors
    Federico Municchi, Nottingham (2020)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "reverseLinear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"

    pimpleControl pimple(mesh);

    #include "readPoroelasticProperties.H"
    #include "createFields.H"
    
    //#include "readStressControls.H"
    const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");
    int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

    scalar convergenceTolerance(readScalar(stressControl.lookup("tolerance")));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while(runTime.run())
    {

        runTime++;
        Info<< "Time: " << runTime.value() << nl << endl;

        while(pimple.loop())
        {
            while (pimple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::ddt(rhof*stor,p) -  fvm::laplacian(Mf,p)
                    ==
                    - fvc::div(phiG)
                    - fvc::ddt(alpha,J)
                    + dims*fvOptions(p)
                );

                fvOptions.constrain(pEqn); //////////////////////////////////////////////////////
                pEqn.relax();
                pEqn.solve();
                fvOptions.correct(p); //////////////////////////////////////////////////////

                if (pimple.finalNonOrthogonalIter())
                {
                    U = fvc::reconstruct(pEqn.flux() + phiG);
                }
            }

            while(pimple.correct())
            {
                int iCorr = 0;
                scalar initialResidual = 0;

                do
                {
                    //- Create stress flux
                    surfaceVectorField phiStress
                    (
                        fvc::interpolate
                        (
                            (mu*T(gradD) + lambda*tr(gradD)*I),
                            "interpolate(stress)"
                        )&mesh.Sf()
                    );

                    while(pimple.correctNonOrthogonal())
                    {
                        fvVectorMatrix DEqn
                        (
                        fvm::d2dt2(rhos,D)
                        - fvm::laplacian(mu,D)
                        ==
                            fvc::div(phiStress)
                        + rhos*g
                        - (alpha*fvc::grad(p))
                        );

                        DEqn.relax();
                        initialResidual = DEqn.solve().max().initialResidual();
                    }

                    gradD = fvc::grad(D);
                    J = scalar(1.0) + fvc::div(D);
                } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);
            }

            Sigma = mu*dev(twoSymm(gradD + T(gradD))) + lambda*tr(gradD)*I;
        }


        runTime.write();

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
