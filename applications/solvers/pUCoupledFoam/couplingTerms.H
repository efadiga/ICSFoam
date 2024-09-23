    // Calculate grad p coupling matrix. Needs to be here if one uses
    // gradient schemes with limiters.  VV, 9/June/2014
//    BlockLduSystem<vector, vector> pInU(fvm::grad(p));

    blockFvMatrix<vector, vector> dUByp(fvm::grad(p));
    // Calculate div U coupling.  Could be calculated only once since
    // it is only geometry dependent.  VV, 9/June/2014
    blockFvMatrix<scalar, vector> dpByU(fvm::div(U));
//    BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));

    // Last argument in insertBlockCoupling says if the column direction
    // should be incremented. This is needed for arbitrary positioning
    // of U and p in the system. This could be better. VV, 30/April/2014

    // Insert momentum equation
    eqSystem.dSByS(0,0).insertEquation(pEqn);
    eqSystem.dSByV(0,0) = dpByU;
    eqSystem.dVByS(0,0) = dUByp;
    eqSystem.dVByV(0,0).insertEquation(UEqn);

