/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "coupledCsrMatrix.H"
#include "coupledMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::coupledCsrMatrix::coupledCsrMatrix(word mode)
:
    csrMatrix(mode),
	source_(nullptr),
	variables_(nullptr)
{
    if (mode.starts_with("h"))
    {
        coupledCsrMatExec_ = cpuCoupledCsrMatrixExecutor();
	}
    else if (mode.starts_with("d"))
    {
        coupledCsrMatExec_ = cudaCoupledCsrMatrixExecutor();
	}
    else
    {
        FatalErrorInFunction
            << "'" << mode << "' is not a valid AMGx execution mode"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * *  Public Member Functions * * * * * * * * * * * *  //

// * * * * * * * * * * * * * * * * Operations * * * * * * * * * * * * * * * //

//- Apply permutation to LDU values (no permutation)
void Foam::coupledCsrMatrix::applyPermutation(const Foam::coupledMatrix& matrix)
{
    // Verify that the permutation has already been computed
    if(!ldu2csrPerm_)
    {
        computePermutation(&(matrix.mesh().lduAddr()));
    }


    // 1) reorder dXbydY pointers according to row major order:
    //         - effects on the first scalar
    //         - effects on the second scalar
    //         - effects on the n scalar
    //         - effects on the first vector
    //         - effects on the second vector
    //         - effects on the n vector
    // 2) Move every Ldu array on the device
    // 3) Use single CUDA kernel to access coefficients in the order given by 1)
    //    and reorder using ldu2Csr permutation while filling values array
    // 4) Free Ldu arrays
    // 5) solve using the new values array

    int nScal = matrix.nScal();
    int nVect = matrix.nVect();
    int nCells = matrix.mesh().nCells();
    int nIntFaces = matrix.mesh().nInternalFaces();
    int totNnz = nCells + 2*nIntFaces;

    this->setNblocks(nScal+nVect*3);


    if(!valuesPtr_)
    {
        std::visit([this, totNnz](const auto& exec)
               { this->valuesPtr_ = exec.template allocZero<scalar>(totNnz*nBlocks_*nBlocks_); },
               this->csrMatExec_);
    }

    // Initialize valuesTmp = [(diag), (upper), (lower)] * nBlocks
    // AoS order
    //scalar* valuesTmp = nullptr;
    //    std::visit([this, &valuesTmp, totNnz](const auto& exec)
    //    		{ valuesTmp = exec.template alloc<scalar>(totNnz*nBlocks_*nBlocks_); },
    //           this->csrMatExec_);


    for(int ds=0; ds<nScal; ds++)
    {
    	for(int bys=0; bys<nScal; bys++)
    	{
            if(matrix.dSBySExists(ds,bys))
    		{
                labelList offsets(1);
                offsets[0] = ds*nBlocks_ + bys;
                const label* offPtr = nullptr;
                std::visit([offsets, &offPtr](const auto& exec)
                           { offPtr = exec.template copyFromFoam<label>(1, offsets.cdata()); },
                           csrMatExec_);

                const blockFvMatrix<scalar,scalar>& blockMat = matrix.dSByS(ds,bys);
                const scalar* diagPtr = nullptr;
                const scalar* upperPtr = nullptr;
                const scalar* lowerPtr = nullptr;

                this->copyLDUPtrs(blockMat,nCells,nIntFaces,&diagPtr,&upperPtr,&lowerPtr);

                this->initializeAndApplyValue
                (
                	nBlocks_,
                    nCells,
                    nIntFaces,
	            	1,
	            	offPtr,
                    ldu2csrPerm_,
                    diagPtr,
                    upperPtr,
                    lowerPtr,
                    valuesPtr_
                );
                std::visit([diagPtr, upperPtr, lowerPtr](const auto& exec)
                           {exec.template clear<scalar>(diagPtr);
                            exec.template clear<scalar>(upperPtr);
                            exec.template clear<scalar>(lowerPtr);}, csrMatExec_);
    		}
    	}
    	for(int byv=0; byv<nVect; byv++)
    	{
    		if(matrix.dSByVExists(ds,byv))
    		{
                labelList offsets(3);
                offsets[0] = ds*nBlocks_ + nScal + 3*byv;
                offsets[1] = ds*nBlocks_ + nScal + 3*byv + 1;
                offsets[2] = ds*nBlocks_ + nScal + 3*byv + 2;
                const label* offPtr = nullptr;
                std::visit([offsets, &offPtr](const auto& exec)
                           { offPtr = exec.template copyFromFoam<label>(3, offsets.cdata()); },
                           csrMatExec_);

                const blockFvMatrix<scalar,vector>& blockMat = matrix.dSByV(ds,byv);
                const scalar* diagPtr = nullptr;
                const scalar* upperPtr = nullptr;
                const scalar* lowerPtr = nullptr;

                this->copyLDUPtrs(blockMat,3*nCells,3*nIntFaces,&diagPtr,&upperPtr,&lowerPtr);

                this->initializeAndApplyValue
                (
                	nBlocks_,
                    nCells,
                    nIntFaces,
	            	3,
	            	offPtr,
                    ldu2csrPerm_,
                    diagPtr,
                    upperPtr,
                    lowerPtr,
                    valuesPtr_
                );
                std::visit([diagPtr, upperPtr, lowerPtr](const auto& exec)
                           {exec.template clear<scalar>(diagPtr);
                            exec.template clear<scalar>(upperPtr);
                            exec.template clear<scalar>(lowerPtr);}, csrMatExec_);
    		}
    	}
    }

    for(int dv=0; dv<nVect; dv++)
    {
    	for(int bys=0; bys<nScal; bys++)
    	{
    		if(matrix.dVBySExists(dv,bys))
    		{
                labelList offsets(3);
                offsets[0] = (nScal+3*dv)*nBlocks_ + bys;
                offsets[1] = (nScal+3*dv)*nBlocks_ + bys + nBlocks_;
                offsets[2] = (nScal+3*dv)*nBlocks_ + bys + 2*nBlocks_;
                const label* offPtr = nullptr;
                std::visit([offsets, &offPtr](const auto& exec)
                           { offPtr = exec.template copyFromFoam<label>(3, offsets.cdata()); },
                           csrMatExec_);

                const blockFvMatrix<vector,vector>& blockMat = matrix.dVByS(dv,bys);
                const scalar* diagPtr = nullptr;
                const scalar* upperPtr = nullptr;
                const scalar* lowerPtr = nullptr;

                this->copyLDUPtrs(blockMat,3*nCells,3*nIntFaces,&diagPtr,&upperPtr,&lowerPtr);

                this->initializeAndApplyValue
                (
                	nBlocks_,
                    nCells,
                    nIntFaces,
	            	3,
	            	offPtr,
                    ldu2csrPerm_,
                    diagPtr,
                    upperPtr,
                    lowerPtr,
                    valuesPtr_
                );
                std::visit([diagPtr, upperPtr, lowerPtr](const auto& exec)
                           {exec.template clear<scalar>(diagPtr);
                            exec.template clear<scalar>(upperPtr);
                            exec.template clear<scalar>(lowerPtr);}, csrMatExec_);
    		}
    	}
    	for(int byv=0; byv<nVect; byv++)
    	{
    		if(matrix.dVByVExists(dv,byv))
    		{
                labelList offsets(9);
                offsets[0] = (nScal+3*dv)*nBlocks_ + nScal + 3*byv;
                offsets[1] = (nScal+3*dv)*nBlocks_ + nScal + 3*byv + 1;
                offsets[2] = (nScal+3*dv)*nBlocks_ + nScal + 3*byv + 2;
                offsets[3] = (nScal+3*dv)*nBlocks_ + nScal + 3*byv + nBlocks_;
                offsets[4] = (nScal+3*dv)*nBlocks_ + nScal + 3*byv + nBlocks_ + 1;
                offsets[5] = (nScal+3*dv)*nBlocks_ + nScal + 3*byv + nBlocks_ + 2;
                offsets[6] = (nScal+3*dv)*nBlocks_ + nScal + 3*byv + 2*nBlocks_;
                offsets[7] = (nScal+3*dv)*nBlocks_ + nScal + 3*byv + 2*nBlocks_ + 1;
                offsets[8] = (nScal+3*dv)*nBlocks_ + nScal + 3*byv + 2*nBlocks_ + 2;
                const label* offPtr = nullptr;
                std::visit([offsets, &offPtr](const auto& exec)
                           { offPtr = exec.template copyFromFoam<label>(9, offsets.cdata()); },
                           csrMatExec_);

                const blockFvMatrix<vector,tensor>& blockMat = matrix.dVByV(dv,byv);
                const scalar* diagPtr = nullptr;
                const scalar* upperPtr = nullptr;
                const scalar* lowerPtr = nullptr;

                this->copyLDUPtrs(blockMat,9*nCells,9*nIntFaces,&diagPtr,&upperPtr,&lowerPtr);

                this->initializeAndApplyValue
                (
                	nBlocks_,
                    nCells,
                    nIntFaces,
	            	9,
	            	offPtr,
                    ldu2csrPerm_,
                    diagPtr,
                    upperPtr,
                    lowerPtr,
                    valuesPtr_
                );
                std::visit([diagPtr, upperPtr, lowerPtr](const auto& exec)
                           {exec.template clear<scalar>(diagPtr);
                            exec.template clear<scalar>(upperPtr);
                            exec.template clear<scalar>(lowerPtr);}, csrMatExec_);
    		}
    	}
    }
}

template<class sourceType, class blockType>
void Foam::coupledCsrMatrix::copyLDUPtrs
(
    const blockFvMatrix<sourceType,blockType>& matrix,
	const label nCells,
	const label nIntFaces,
    const scalar** diagPtr,
    const scalar** upperPtr,
    const scalar** lowerPtr
)
{

    std::visit([matrix, nCells, nIntFaces, &diagPtr, &upperPtr, &lowerPtr](const auto& exec)
               {
    	           if (matrix.hasDiag())
    	           {
    	               *(diagPtr) = exec.template copyFromFoam<scalar>
                       (
                           nCells,
				    	   reinterpret_cast<const scalar*>(matrix.diag().cdata())
				       );
    	           }

    	           if (matrix.hasUpper())
    	           {
    	               *(upperPtr) = exec.template copyFromFoam<scalar>
                       (
                           nIntFaces,
				    	   reinterpret_cast<const scalar*>(matrix.upper().cdata())
				       );
    	           }
    	           if (matrix.hasLower())
    	           {
    	               *(lowerPtr) = exec.template copyFromFoam<scalar>
                       (
                           nIntFaces,
				    	   reinterpret_cast<const scalar*>(matrix.lower().cdata())
				       );
    	           }
               }, csrMatExec_);
}


void Foam::coupledCsrMatrix::fillSource
(
    const coupledMatrix& matrix,
    const PtrList<scalarField>& sSource,
	const PtrList<vectorField>& vSource
)
{
    const label nCells = matrix.mesh().nCells();
    const label nScalar = matrix.nScal();
    const label nVector = matrix.nVect();

    forN(nScalar,j)
    {
    	fillField<scalar>
    	(
            nCells,
			j,
			sSource[j],
			source_
    	);
    }

    forN(nVector,j)
    {
    	fillField<vector>
    	(
            nCells,
			nScalar + 3*j,
			vSource[j],
			source_
    	);
    }
}

void Foam::coupledCsrMatrix::fillVariables
(
    const coupledMatrix& matrix,
    const PtrList<volScalarField>& sVolField,
	const PtrList<volVectorField>& vVolField
)
{
    const label nCells = matrix.mesh().nCells();
    const label nScalar = matrix.nScal();
    const label nVector = matrix.nVect();

    forN(nScalar,j)
    {
    	fillField<scalar>
    	(
            nCells,
			j,
			sVolField[j].primitiveField(),
			variables_
    	);
    }
    forN(nVector,j)
    {
    	fillField<vector>
    	(
            nCells,
			nScalar + 3*j,
			vVolField[j].primitiveField(),
			variables_
    	);
    }
}

void Foam::coupledCsrMatrix::transferVariables
(
    const coupledMatrix& matrix,
    PtrList<volScalarField>& sVolField,
	PtrList<volVectorField>& vVolField
)
{
    const label nCells = matrix.mesh().nCells();
    const label nScalar = matrix.nScal();
    const label nVector = matrix.nVect();

    scalar* hostPtr = new scalar[nCells*nBlocks_];

    std::visit([this, &hostPtr, nCells](const auto& exec)
               { exec.template copyToFoam<scalar>(nCells*this->nBlocks_, this->variables_, &hostPtr);},
               csrMatExec_);

    forN(nScalar,j)
    {
        for(label i=0; i<nCells;i++)
        {
            sVolField[j].primitiveFieldRef()[i] = hostPtr[i*nBlocks_ + j];
        }
    }
    forN(nVector,j)
    {
        for(label i=0; i<nCells;i++)
        {
            vVolField[j].primitiveFieldRef()[i].x() = hostPtr[i*nBlocks_ + nScalar + 3*j];
            vVolField[j].primitiveFieldRef()[i].y() = hostPtr[i*nBlocks_ + nScalar + 3*j + 1];
            vVolField[j].primitiveFieldRef()[i].z() = hostPtr[i*nBlocks_ + nScalar + 3*j + 2];
        }
    }
}

////- Apply permutation from LDU to CSR considering the interface values
//void Foam::csrMatrix:: applyPermutation
//(
//    const lduMatrix& lduMatrix,
//    const label diagIndexGlobal,
//    const label lowOffGlobal,
//    const label uppOffGlobal,
//    const labelList& extRows,
//    const labelList& extCols,
//    const scalarField& extVals
//)
//{
//    // Verify that the permutation has already been computed
//    if(!ldu2csrPerm_)
//    {
//        computePermutation(
//            &(lduMatrix.lduAddr()),
//            diagIndexGlobal,
//            lowOffGlobal,
//            uppOffGlobal,
//            extRows,
//            extCols
//        );
//    }
//
//    const scalarField& diag = lduMatrix.diag();
//    const scalarField& upper = lduMatrix.upper();
//    const scalarField& lower = lduMatrix.lower();
//
//    label nIntFaces = upper.size();
//    label nCells = diag.size();
//    label nnzExt = extVals.size();
//    label totNnz = nCells + 2*nIntFaces + nnzExt;
//
//    if(!valuesPtr_)
//    {
//        valuesPtr_ = new scalarField(totNnz);
//    }
//
//    // Initialize valuesTmp = [(diag), (upper), (lower), (extValues)]
//    scalarField valuesTmp(totNnz);
//
//    initializeValueExt
//    (
//        nCells,
//        nIntFaces,
//        nnzExt,
//        diag.cdata(),
//        upper.cdata(),
//        lower.cdata(),
//        extVals.cdata(),
//        valuesTmp.data()
//    );
//
//    // Apply permutation
//    applyValuePermutation
//    (
//        totNnz,
//        ldu2csrPerm_->cdata(),
//        valuesTmp.cdata(),
//        valuesPtr_->data()
//    );
//}
//
////- Apply permutation from LDU to CSR considering the interface values
//void Foam::csrMatrix:: applyPermutation
//(
//    const lduMatrix& lduMatrix,
//    const scalarField& extVals
//)
//{
//    const scalarField& diag = lduMatrix.diag();
//    const scalarField& upper = lduMatrix.upper();
//    const scalarField& lower = lduMatrix.lower();
//
//    label nIntFaces = upper.size();
//    label nCells = diag.size();
//    label nnzExt = extVals.size();
//    label totNnz = nCells + 2*nIntFaces + nnzExt;
//
//    if(!valuesPtr_)
//    {
//        valuesPtr_ = new scalarField(totNnz);
//    }
//
//    // Initialize valuesTmp = [(diag), (upper), (lower), (extValues)]
//    scalarField valuesTmp(totNnz);
//
//    initializeValueExt
//    (
//        nCells,
//        nIntFaces,
//        nnzExt,
//        diag.cdata(),
//        upper.cdata(),
//        lower.cdata(),
//        extVals.cdata(),
//        valuesTmp.data()
//    );
//
//    // Apply permutation
//    applyValuePermutation
//    (
//        totNnz,
//        ldu2csrPerm_->cdata(),
//        valuesTmp.cdata(),
//        valuesPtr_->data()
//    );
//}

// * * * * * * * * * * * * * Explicit instantiations  * * * * * * * * * * * //

// ************************************************************************* //
