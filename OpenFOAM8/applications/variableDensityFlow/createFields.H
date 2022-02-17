Info<< "Reading field u\n" << endl; // electric potential
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
const dimensionedScalar a(transportProperties.lookup("a"));

dimensionedVector g( "g", dimensionSet(0,1,-2,0,0,0,0), vector(0.0, -9.81, 0.0));
