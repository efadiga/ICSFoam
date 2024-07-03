/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
    Copyright (C) 2022-2023 Cineca
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

\*---------------------------------------------------------------------------*/

#include "AMGXSolver.H"

#include "direction.H"
#include "AmgXLinearSolverContext.H"
#include "linearSolverContextTable.H"
#include "coupledCsrMatrix.H"

#include "globalIndex.H"

// #include <iostream>
// #include <fstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(AMGXSolver, 0);

	coupledMatrix::solver::adddictionaryConstructorToTable<AMGXSolver>
    	addAMGSolverDictionaryConstructorToTable_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AMGXSolver::AMGXSolver
(
    const dictionary& dict,
    const coupledMatrix& matrix
)
:
	coupledMatrix::solver(typeName, dict, matrix)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- construct the matrix for AmgX
/*void Foam::AmgXSolver::buildAndApplyMatrixPermutation
(
    deviceCsrMatrix* csrMatrix,
    label& nRowsGlobal
) const
{
    const lduInterfacePtrsList interfaces(this->matrix_.mesh().openFoamMesh().interfaces());

    // Local degrees-of-freedom i.e. number of local rows
    const label nLocalRows = this->matrix_.mesh().nCells();
    label nRowsLocal = nLocalRows;
    nRowsGlobal = returnReduce(nRowsLocal, sumOp<label>());

    // Number of internal faces (connectivity)
    const label nIntFaces = this->matrix_.mesh().nInternalFaces();

    const globalIndex globalNumbering(nLocalRows);

    const label diagIndexGlobal = globalNumbering.toGlobal(0);
    const label firstLowerInd = this->matrix_.mesh().openFoamMesh().owner()[0]; //lower.cdata()[0];
    label lowOffGlobal = globalNumbering.toGlobal(firstLowerInd) - firstLowerInd;
    const label firstUpperInd = this->matrix_.mesh().openFoamMesh().neighbour()[0]; // upper.cdata()[0];
    label uppOffGlobal = globalNumbering.toGlobal(firstUpperInd) - firstUpperInd;

    labelList globalCells
    (
        identity
        (
            globalNumbering.localSize(),
            globalNumbering.localStart()
        )
    );

    // Connections to neighbouring processors
    const label nReq = Pstream::nRequests();

    label nProcValues = 0;

    // Initialise transfer of global cells
    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi))
        {
            nProcValues += interfaces[patchi].faceCells().size(); //lduAddr.patchAddr(patchi).size();

            interfaces[patchi].initInternalFieldTransfer
            (
                Pstream::commsTypes::nonBlocking,
                globalCells
            );
        }
    }

    if (Pstream::parRun())
    {
        Pstream::waitRequests(nReq);
    }

    deviceField<label> procRows(nProcValues, 0);
    deviceField<label> procCols(nProcValues, 0);
    deviceField<scalar> procVals(nProcValues, Foam::Zero);
    nProcValues = 0;

    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi))
        {
            // Processor-local values
            const label len = interfaces[patchi].faceCells().size();
            const deviceField<label> faceCells(len, interfaces[patchi].faceCells().cdata());
            const deviceField<scalar>& bCoeffs = this->interfaceBouCoeffs_[patchi];

            labelList nbrCells
            (
                interfaces[patchi].internalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    globalCells
                )
            );

            if (faceCells.size() != nbrCells.size())
            {
                FatalErrorInFunction
                    << "Mismatch in interface sizes (AMI?)" << nl
                    << "Have " << faceCells.size() << " != "
                    << nbrCells.size() << nl
                    << exit(FatalError);
            }

            procRows.copy(faceCells, nProcValues, 0);
            procCols.copyIn(nbrCells, len, nProcValues);
            procVals.copy(bCoeffs, nProcValues, 0);

            nProcValues += len;
        }
    }

    procVals.negate();  // Change sign for entire field (see previous note)

    csrMatrix->applyPermutation
    (
        this->matrix_,
        diagIndexGlobal,
        lowOffGlobal,
        uppOffGlobal,
        procRows,
        procCols,
        procVals
    );

    DebugInfo<< "Converted LDU matrix to CSR format" << nl;

}


//- construct the matrix for AmgX
void Foam::AmgXSolver::applyMatrixPermutation
(
    deviceCsrMatrix* csrMatrix,
    label& nRowsGlobal
) const
{
    const UPtrList<const devicelduInterfaceField>& interfaces = this->interfaces_;

    // Local degrees-of-freedom i.e. number of local rows
    const label nLocalRows = this->matrix_.mesh().nCells();
    label nRowsLocal = nLocalRows;
    nRowsGlobal = returnReduce(nRowsLocal, sumOp<label>());

    //- Number of internal faces (connectivity)
    const label nIntFaces = this->matrix_.mesh().nInternalFaces();

    //- Connections to neighbouring processors
    const label nReq = Pstream::nRequests();

    label nProcValues = 0;

    // Initialise transfer of global cells
    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi)) nProcValues += this->interfaceBouCoeffs_[patchi].size();
    }

    if (Pstream::parRun())
    {
        Pstream::waitRequests(nReq);
    }

    deviceField<scalar> procVals(nProcValues, Foam::Zero);
    nProcValues = 0;

    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi))
        {
            //- Processor-local values
            const deviceField<scalar>& bCoeffs = this->interfaceBouCoeffs_[patchi];
            const label len = bCoeffs.size();

            procVals.copy(bCoeffs, nProcValues, 0);

            nProcValues += len;
        }
    }

    procVals.negate();  // Change sign for entire field (see previous note)

    csrMatrix->applyPermutation
    (
        this->matrix_,
        procVals
    );

    DebugInfo<< "Converted LDU matrix values to CSR format" << nl;

}*/

residualsIO AMGXSolver::solve
(
    PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
    const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource
) const
{
    const int nScalar = this->matrix().nScal();
    const int nVector = this->matrix().nVect();
	residualsIO solverPerf(nScalar,nVector);
	NotImplemented;
return solverPerf;}
//    const fvMesh& fvm = dynamicCast<const fvMesh>(this->matrix_.mesh().thisDb());
//
//    label nCells = this->matrix_.mesh().nCells();
//
//    const linearSolverContextTable<AmgXLinearSolverContext<coupledCsrMatrix>>& contexts =
//        linearSolverContextTable<AmgXLinearSolverContext<coupledCsrMatrix>>::New(fvm);
//
//    AmgXLinearSolverContext<coupledCsrMatrix>& ctx = contexts.getContext("flowSolver");
//
//    if (!ctx.loaded())
//    {
//        FatalErrorInFunction
//            << "Could not initialize AMGx" << nl << abort(FatalError);
//    }
//
//    AmgXWrapper& amgx = ctx.amgx_;
//
//    coupledCsrMatrix& Amat = ctx.Amat_;
//
//    label nGlobalCells;
//
//    if(!Pstream::parRun())
//    {
//        nGlobalCells = nCells;
//        Amat.applyPermutation(this->matrix_);
//    }
//    /*else
//    {
//        if(!Amat.hasPermutation()) buildAndApplyMatrixPermutation(&Amat, nGlobalCells);
//        else applyMatrixPermutation(&Amat, nGlobalCells);
//    }*/
//
//    label nnz = Amat.nnz();
//
//    if(!ctx.initialized())
//    {
//        Info<< "Initializing AmgX Linear Solver " << eqName_ << nl;
//
//        amgx.setOperator(nGlobalCells, &Amat);
//
//        // Amat.clearAddressing();
//
//        ctx.initialized() = true;
//    }
//    else
//    {
//        amgx.updateOperator(&Amat);
//    }
//
//    amgx.solve(psi.data(), source.cdata(), &Amat);
//
//    scalarField iNorm(nScalar+3*nVector, 0.0);
//    amgx.getResidual(0, iNorm);
//    solverPerf.initialResidual() = iNorm[0];
//
//    label nIters = 0;
//    amgx.getIters(nIters);
//    solverPerf.nIterations() = nIters;
//
//    scalarField fNorm(nScalar+3*nVector, 0.0);
//    amgx.getResidual(nIters, fNorm);
//    solverPerf.finalResidual() = fNorm[0];
//
//    return solverPerf;
//
//}

residualsIO AMGXSolver::solveDelta
(
    PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
    const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource
) const
{
    const int nScalar = this->matrix().nScal();
    const int nVector = this->matrix().nVect();

	residualsIO solverPerf(nScalar,nVector);
	NotImplemented;
    return solverPerf;}
	//
//    // Allocate variables to hold solved increment
//    PtrList<volScalarField> dsW(nScalar);
//    PtrList<volVectorField> dvW(nVector);
//
//    const fvMesh& fvm = dynamicCast<const fvMesh>(this->matrix_.mesh().thisDb());
//
//    label nCells = this->matrix_.mesh().nCells();
//
//    const linearSolverContextTable<AmgXLinearSolverContext<coupledCsrMatrix>>& contexts =
//        linearSolverContextTable<AmgXLinearSolverContext<coupledCsrMatrix>>::New(fvm);
//
//    AmgXLinearSolverContext<coupledCsrMatrix>& ctx = contexts.getContext("flowSolver");
//
//    if (!ctx.loaded())
//    {
//        FatalErrorInFunction
//            << "Could not initialize AMGx" << nl << abort(FatalError);
//    }
//
//    AmgXWrapper& amgx = ctx.amgx_;
//
//    coupledCsrMatrix& Amat = ctx.Amat_;
//
//    label nGlobalCells;
//
//    if(!Pstream::parRun())
//    {
//        nGlobalCells = nCells;
//        Amat.applyPermutation(this->matrix_);
//    }
//    /*else
//    {
//        if(!Amat.hasPermutation()) buildAndApplyMatrixPermutation(&Amat, nGlobalCells);
//        else applyMatrixPermutation(&Amat, nGlobalCells);
//    }*/
//
//    label nnz = Amat.nnz();
//
//    if(!ctx.initialized())
//    {
//        Info<< "Initializing AmgX Linear Solver " << eqName_ << nl;
//
//        amgx.setOperator(nGlobalCells, &Amat);
//
//        // Amat.clearAddressing();
//
//        ctx.initialized() = true;
//    }
//    else
//    {
//        amgx.updateOperator(&Amat);
//    }
//
//    amgx.solve(psi.data(), source.cdata(), &Amat);
//
//    scalarField iNorm(nScalar+3*nVector, 0.0);
//    amgx.getResidual(0, iNorm);
//    solverPerf.initialResidual() = iNorm[0];
//
//    label nIters = 0;
//    amgx.getIters(nIters);
//    solverPerf.nIterations() = nIters;
//
//    scalarField fNorm(nScalar+3*nVector, 0.0);
//    amgx.getResidual(nIters, fNorm);
//    solverPerf.finalResidual() = fNorm[0];
//
//    forN(nScalar,i)
//	{
//		sW[i].primitiveFieldRef() += dsW[i].primitiveFieldRef();
//	}
//	forN(nVector,i)
//	{
//		vW[i].primitiveFieldRef() += dvW[i].primitiveFieldRef();
//	}
//
//    return solverPerf;
//
//}


residualsIO AMGXSolver::solveDelta
(
    PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
    const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource,
	PtrList<volScalarField>& dsW, PtrList<volVectorField>& dvW
) const
{
    const int nScalar = this->matrix().nScal();
    const int nVector = this->matrix().nVect();

	residualsIO solverPerf(nScalar,nVector);

    forN(nScalar,i)
    {
    	dsW[i].primitiveFieldRef() = sW[i].primitiveField();
    	dsW[i].boundaryFieldRef() = sW[i].boundaryField();
    }
    forN(nVector,i)
    {
    	dvW[i].primitiveFieldRef() = vW[i].primitiveField();
    	dvW[i].boundaryFieldRef() = vW[i].boundaryField();
    }

    // Allocate temp storage
    PtrList<scalarField> sTmp(nScalar);
    PtrList<vectorField> vTmp(nVector);

    forN(nScalar,j)
    {
        sTmp.set(j, new scalarField(sW[j].size(), 0.0));
    }

    forN(nVector,j)
    {
    	vTmp.set(j, new vectorField(vW[j].size(), Zero));
    }

    scalarList sNormFactor(nScalar);
    scalarList vNormFactor(nVector);

    forN(nScalar,i)
    {
        scalar avg = gAverage(sW[i]);
        dsW[i].primitiveFieldRef() -= avg;
        dsW[i].boundaryFieldRef() -= avg;
    }
    forN(nVector,i)
    {
        vector avg = gAverage(vW[i]);
        dvW[i].primitiveFieldRef()-= avg;
        dvW[i].boundaryFieldRef() -= avg;
    }

    // Matrix multiplication
    this->matrix_.matrixMul(dsW, dvW, sTmp, vTmp);

    forN(nScalar,i)
    {
        dsW[i].primitiveFieldRef() = 0.0;
    }
    forN(nVector,i)
    {
        dvW[i].primitiveFieldRef() = vector::zero;
    }

    forAll(sNormFactor, i)
    {
        sNormFactor[i] = gSum(mag(sTmp[i]) + mag(sSource[i])) + VSMALL;
        solverPerf.sInitRes()[i] = gSumMag(sSource[i]) / sNormFactor[i];
        solverPerf.sFinalRes()[i] = solverPerf.sInitRes()[i];
    }

    // Use vector magnitude for normalisation
    forAll(vNormFactor, i)
    {
        vNormFactor[i] = gSum(mag(vTmp[i]) + mag(vSource[i])) + VSMALL;
        solverPerf.vInitRes()[i] = gSumCmptMag(vSource[i])/vNormFactor[i];
        solverPerf.vFinalRes()[i] = solverPerf.vInitRes()[i];
    }


    /// AMGX SOLUTION START ///

    const fvMesh& fvm = dynamicCast<const fvMesh>(this->matrix_.mesh().thisDb());

    label nCells = this->matrix_.mesh().nCells();

    const linearSolverContextTable<AmgXLinearSolverContext<coupledCsrMatrix>>& contexts =
        linearSolverContextTable<AmgXLinearSolverContext<coupledCsrMatrix>>::New(fvm);

    AmgXLinearSolverContext<coupledCsrMatrix>& ctx = contexts.getContext("flowSolver",typeName);

    if (!ctx.loaded())
    {
        FatalErrorInFunction
            << "Could not initialize AMGx" << nl << abort(FatalError);
    }

    AmgXWrapper& amgx = ctx.amgx_;

    coupledCsrMatrix& Amat = ctx.Amat_;

    label nGlobalCells;

    if(!ctx.initialized())
    {
        if(Pstream::parRun()) amgx.initialiseMatrixComms(&Amat);
        else Amat.setGpuProc(true);
    }
    nGlobalCells = Amat.applyPermutation(this->matrix_);
    //}
    /*else
    {
        if(!Amat.hasPermutation()) buildAndApplyMatrixPermutation(&Amat, nGlobalCells);
        else applyMatrixPermutation(&Amat, nGlobalCells);
    }*/

    if(!ctx.initialized())
    {
        Info<< "Initializing AmgX Linear Solver (Coupled)"  << nl;

        amgx.setOperator(nGlobalCells, &Amat);

        Amat.clearAddressing();

        ctx.initialized() = true;
        amgx.updateOperator(&Amat);
    }
    else
    {
        amgx.updateOperator(&Amat);
    }

    Info << "Allocating source and variables on the device" << endl;
    //Amat.allocSource(this->matrix());
    //Amat.allocVariables(this->matrix());

    Amat.fillSource(this->matrix(),sSource,vSource);
    Amat.fillVariables(this->matrix(),dsW,dvW);

    Info << "Solving " << endl;
    amgx.solve(&Amat);

    //Amat.clearSource();

    Amat.transferVariables(this->matrix(),dsW,dvW);
    //Amat.clearVariables();

    scalarField iNorm(nScalar+3*nVector, 0.0);
    amgx.getResidual(0, iNorm, Amat.nBlocks());
    //solverPerf.initialResidual() = iNorm[0];

    Info << "AMGX initial residual: " << iNorm << endl;

    label nIters = 0;
    amgx.getIters(nIters);
    //solverPerf.nIterations() = nIters;
    Info << "AMGX nIterations: " << nIters << endl;

    scalarField fNorm(nScalar+3*nVector, 0.0);
    amgx.getResidual(nIters, fNorm, Amat.nBlocks());
    //solverPerf.finalResidual() = fNorm[0];

    Info << "AMGX final residual: " << fNorm << endl;

    /// AMGX SOLUTION END ///


    // Re-calculate the residual
	this->matrix_.matrixMul(dsW, dvW, sTmp, vTmp);

    forN(nScalar,j)
    {
        forAll(sTmp[j], iCell)
        {
            sTmp[j][iCell] = sSource[j][iCell] - sTmp[j][iCell];
        }
    }
    forN(nVector,j)
    {
        forAll(vTmp[j], iCell)
        {
            vTmp[j][iCell] = vSource[j][iCell] - vTmp[j][iCell];
        }
    }

    forN(nScalar,i)
    {
    	solverPerf.sFinalRes()[i] = gSumMag(sTmp[i]) / sNormFactor[i];
    }

    forN(nVector,i)
    {
    	vector res = gSumCmptMag(vTmp[i])/vNormFactor[i];

        // Zero the residual in non-solved directions
        vector::labelType validComponents = this->matrix_.mesh().solutionD(); //-1 for empty directions

        forAll(validComponents, cmpt)
        {
            if (validComponents[cmpt] == -1)
            {
                res[cmpt] = 0.0;
            }
        }

        solverPerf.vFinalRes()[i] = res;
    }


    return solverPerf;
}



} // End namespace Foam

// * * * * * * * * * * * * * Explicit instantiations  * * * * * * * * * * * //


// ************************************************************************* //
