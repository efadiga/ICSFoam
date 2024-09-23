    // Momentum predictor

  //  MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
    //    fvm::ddt(U)
       fvm::div(phi, U)
      + turbulence->divDevReff(U)
	  - fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    volScalarField rAU(1.0/UEqn.A());
    tmp<volScalarField> rAtU(rAU);
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    adjustPhi(phiHbyA, U, p);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU());


    UEqn.relax();


//    fvOptions.constrain(UEqn);

//    if (simple.momentumPredictor())
//    {
//        solve(UEqn == -fvc::grad(p));
//
//        fvOptions.correct(U);
//    }
