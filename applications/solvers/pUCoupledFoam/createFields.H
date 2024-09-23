PtrList<volScalarField> scalarVars(1);
PtrList<volVectorField> vectorVars(1);


Info<< "Reading field p\n" << endl;
scalarVars.set
(
    0,
    new volScalarField
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

volScalarField& p = scalarVars[0];

Info<< "Reading field U\n" << endl;
vectorVars.set
(
    0,
    new volVectorField
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
    )
);

volVectorField& U = vectorVars[0];

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);


#include "createFvOptions.H"
