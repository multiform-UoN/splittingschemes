Info<< "Reading field V\n" << endl; // electric potential
volScalarField u
(
    IOobject
    (
        "u",
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

Info<< "Reading field v\n" << endl; // concentration
volScalarField v
(
    IOobject
    (
        "v",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

Info<< "Reading/calculating face flux field phi\n" << endl; // #include "createPhi.H" 
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

/*----------------------------------------------------------------------------------------*/

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


const dimensionedScalar LC(electrokineticProperties.lookup("L"));

//- Read the Cprop dictionary
const dictionary Cprop(electrokineticProperties.subDict("Cprop"));
const dimensionedScalar m(Cprop.lookup("m"));
const dimensionedScalar a(Cprop.lookup("a"));
const dimensionedScalar z1(Cprop.lookup("z1"));
const dimensionedScalar z2(Cprop.lookup("z2"));


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
        dimensionSet(0,-2,0,0,0,0,0),
        scalar(0)
    )
);
