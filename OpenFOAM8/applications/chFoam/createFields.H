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
