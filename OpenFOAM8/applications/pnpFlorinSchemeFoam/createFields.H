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

Info<< "Reading field C\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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


const dimensionedScalar L(electrokineticProperties.lookup("L"));

//- Read fluid relative dielectric
const dimensionedScalar epsilonR_(electrokineticProperties.lookup("epsilonR"));

//- Calculate effective dielectric
const dimensionedScalar epsilon_(epsilonR_);

//- Read the Cprop dictionary
const dictionary Cprop(electrokineticProperties.subDict("Cprop"));

const dimensionedScalar D(Cprop.lookup("D"));
const dimensionedScalar Dphi(Cprop.lookup("Dphi"));
const dimensionedScalar Z(Cprop.lookup("Z"));


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
