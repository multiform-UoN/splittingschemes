Info<< "Reading field V\n" << endl; // electric potential
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

volScalarField U
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

volScalarField dW
(
    IOobject
    (
        "dW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(
        "dW",
        dimensionSet(0,0,0,0,0,0,0),
        scalar(0)
    )
);

volScalarField ddW
(
    IOobject
    (
        "ddW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(
        "ddW",
        dimensionSet(0,0,0,0,0,0,0),
        scalar(0)
    )
);

// label VRefCell = 0;
// scalar VRefValue = 0.0;
// setRefCell(
//     V,
//     mesh.solutionDict().subDict("couplingControls"),
//     VRefCell,
//     VRefValue
// );
/*----------------------------------------------------------------------------------------*/

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar m(transportProperties.lookup("m"));
const dimensionedScalar epsilonsquare(transportProperties.lookup("epsilonsquare"));
