//- V is the electric potential
Info<< "Reading field V\n" << endl;
volScalarField V
(
    IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// label VRefCell = 0;
// scalar VRefValue = 0.0;
// setRefCell(
//     V,
//     mesh.solutionDict().subDict("couplingControls"),
//     VRefCell,
//     VRefValue
// );



//- U is the velocity field
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



//#include "createPhi.H" 
Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);



Info<< "Reading/calculating Nernst-Planck face flux field phiNP\n" << endl;
surfaceScalarField phiNP
(
    IOobject
    (
        "phiNP",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    phi*scalar(0)
);



Info<< "Reading electrokineticProperties\n" << endl;
IOdictionary electrokineticProperties
(
    IOobject
    (
        "electrokineticProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);



//- Read fluid relative dielectric
const dimensionedScalar epsilonR_
(
    electrokineticProperties.lookup("epsilonR")
);



//- Calculate effective dielectric
const dimensionedScalar epsilon_(epsilonR_);



//- Read chemical species
PtrList<dictionary> species;

PtrList<entry> specEntries
(
    electrokineticProperties.lookup("species")
);

species.setSize(specEntries.size());

//- Create concentration fields
PtrList<volScalarField> cPtrL(species.size());



Info<< "Reading ions concentrations\n" << endl;
forAll(species,sp)
{
    species.set(sp, new dictionary(specEntries[sp].dict()));
    const dimensionedScalar D(species[sp].lookup("D"));
    const dimensionedScalar Dphi(species[sp].lookup("Dphi"));
    const dimensionedScalar Z(species[sp].lookup("Z"));

    const surfaceScalarField phiNP("phiNP",-fvc::flux(fvc::grad(V))*(Dphi)*Z);

    cPtrL.set(
        sp,
        new volScalarField
        (
            IOobject
            (
                specEntries[sp].keyword(),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )

    );

}



//- Create charge density field
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(
        "rho",
        dimensionSet(0,0,-1,0,0,0,0),
        scalar(0)
    )
);
