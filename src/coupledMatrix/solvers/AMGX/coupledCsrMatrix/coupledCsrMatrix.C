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
    csrMatrix(mode)
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
Foam::label Foam::coupledCsrMatrix::applyPermutation(const Foam::coupledMatrix& matrix)
{
    label nScal = matrix.nScal();
    label nVect = matrix.nVect();

    this->setNblocks(nScal+nVect*3);
	label nnzExt = 0;
    const lduInterfacePtrsList& interfaces = matrix.mesh().interfaces();
    // Verify that the permutation has already been computed
    if(!hasPermutation())
    {
    	const lduAddressing& lduAddr = matrix.mesh().lduAddr();
        computePermutation(lduAddr, interfaces, nnzExt);
    }
    else
    {
        forAll(interfaces, patchi)
        {
            if (interfaces.set(patchi))
            {
            	// Processor-local values
                const labelUList& faceCells = matrix.mesh().lduAddr().patchAddr(patchi);
                const label len = faceCells.size();
                nnzExt += len;
            }
        }
    }


    //nnzExt *= nBlocks_;
    label nCells = matrix.mesh().nCells();
    label nIntFaces = matrix.mesh().nInternalFaces();
    label totNnz;

    //- Compute global number of equations
    label nGlobalCells = returnReduce(nCells, sumOp<label>());

    if(consolidationStatus_ == ConsolidationStatus::initialized)
    {
        totNnz = nConsRows_ + 2*nConsIntFaces_ + nConsExtNz_;
    }
    else
    {
        totNnz = nCells + 2*nIntFaces + nnzExt;
    }

    if(gpuProc_)
    {
        if(!valuesPtr_)
        {
            std::visit([this, totNnz](const auto& exec)
                   { this->valuesPtr_ = exec.template allocZero<scalar>(totNnz*nBlocks_*nBlocks_); },
                   this->csrMatExec_);
        }
    }
// Applying permutation block by block

    for(int ds=0; ds<nScal; ds++)
    {
    	for(int bys=0; bys<nScal; bys++)
    	{
            if(matrix.dSBySExists(ds,bys))
    		{
                const blockFvMatrix<scalar,scalar>& blockMat = matrix.dSByS(ds,bys);
                labelList offsets(1);
                offsets[0] = ds*nBlocks_ + bys;
                blockApplyPermutation(matrix,blockMat,offsets,nCells,nIntFaces,nnzExt);
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
                const blockFvMatrix<scalar,vector>& blockMat = matrix.dSByV(ds,byv);
                blockApplyPermutation(matrix,blockMat,offsets,nCells,nIntFaces,nnzExt);
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
                const blockFvMatrix<vector,vector>& blockMat = matrix.dVByS(dv,bys);
                blockApplyPermutation(matrix,blockMat,offsets,nCells,nIntFaces,nnzExt);
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
                const blockFvMatrix<vector,tensor>& blockMat = matrix.dVByV(dv,byv);
                blockApplyPermutation(matrix,blockMat,offsets,nCells,nIntFaces,nnzExt);
    		}
    	}
    }

    return nGlobalCells;
}

template<class sourceType, class blockType>
void Foam::coupledCsrMatrix::blockApplyPermutation
(
	const coupledMatrix& cpldMatrix,
    const blockFvMatrix<sourceType,blockType>& matrix,
	const labelList offsets,
	const label cells,
	const label intFaces,
	const label nnzExt
)
{
    const label* offPtr = nullptr;
    const label offSize = offsets.size();
    if (gpuProc_)
    {
        std::visit([offsets, &offPtr](const auto& exec)
                   { offPtr = exec.template copyFromFoam<label>(offsets.size(), offsets.cdata()); },
                   csrMatExec_);
    }
    scalar* diagPtr = nullptr;
    scalar* upperPtr = nullptr;
    scalar* lowerPtr = nullptr;
    scalar* extPtr = nullptr;

    if(consolidationStatus_ == ConsolidationStatus::initialized)
    {
        blkConsInit
 	    (
 	       cpldMatrix,
           matrix,
 	 	   cells,
 	 	   intFaces,
 	 	   nnzExt,
 	 	   offSize,
 	 	   &diagPtr,
 	 	   &upperPtr,
 	 	   &lowerPtr,
 	 	   &extPtr
 	    );
    }
    else
    {
        copyLDUPtrs
 	    (
 	       cpldMatrix,
           matrix,
 	 	   cells,
 	 	   intFaces,
 	 	   nnzExt,
 	 	   offSize,
 	 	   &diagPtr,
 	 	   &upperPtr,
 	 	   &lowerPtr,
 	 	   &extPtr
 	    );
    }
    if(gpuProc_)
    {
        this->initializeAndApplyValue
        (
            nBlocks_,
            nConsRows_,
            nConsIntFaces_,
			nConsExtNz_,
        	offSize,
        	offPtr,
            ldu2csrPerm_,
            diagPtr,
            upperPtr,
            lowerPtr,
			extPtr,
            valuesPtr_
        );
    }
    std::visit([diagPtr, upperPtr, lowerPtr, offPtr, extPtr](const auto& exec)
               {exec.template clear<scalar>(diagPtr);
                exec.template clear<scalar>(upperPtr);
                exec.template clear<scalar>(lowerPtr);
                exec.template clear<label>(offPtr);
                exec.template clear<scalar>(extPtr);}, csrMatExec_);

}

template<class sourceType, class blockType>
void Foam::coupledCsrMatrix::copyLDUPtrs
(
	const coupledMatrix& cpldMatrix,
    const blockFvMatrix<sourceType,blockType>& matrix,
	const label nCells,
	const label nIntFaces,
	const label nnzExt,
	const label offSize,
    scalar** diagPtr,
    scalar** upperPtr,
    scalar** lowerPtr,
    scalar** extPtr
)
{
    Field<blockType> foamExtVals(nnzExt,Foam::Zero);
    label localNnz = 0;
    forAll(cpldMatrix.mesh().boundary(), patchi)
    {
		if (cpldMatrix.mesh().boundary()[patchi].coupled() && matrix.interfacesUpper().set(patchi))
		{
            //- Processor-local values
            const Field<blockType>& bCoeffs = matrix.interfacesUpper()[patchi];
            const label len = bCoeffs.size();

            SubList<blockType>(foamExtVals, len, localNnz) = bCoeffs;
            localNnz += len;
        }
    }
    if (localNnz != nnzExt && !matrix.diagonal())
    {
    	Info << "Warning: dimension mismatch in interface consolidation" << endl;
    }

    std::visit([matrix, foamExtVals, offSize, nCells, nIntFaces, nnzExt, &diagPtr, &upperPtr, &lowerPtr, &extPtr]
			   (const auto& exec)
               {
    	           if (matrix.hasDiag())
    	           {
    	               *(diagPtr) = const_cast<scalar*>(exec.template copyFromFoam<scalar>
                       (
                           offSize*nCells,
				    	   reinterpret_cast<const scalar*>(matrix.diag().cdata())
				       ));
    	           }
    	           if (matrix.hasUpper())
    	           {
    	               *(upperPtr) = const_cast<scalar*>(exec.template copyFromFoam<scalar>
                       (
                           offSize*nIntFaces,
				    	   reinterpret_cast<const scalar*>(matrix.upper().cdata())
				       ));
    	           }
    	           if (matrix.hasLower())
    	           {
    	               *(lowerPtr) = const_cast<scalar*>(exec.template copyFromFoam<scalar>
                       (
                           offSize*nIntFaces,
				    	   reinterpret_cast<const scalar*>(matrix.lower().cdata())
				       ));
    	           }
    	           if (foamExtVals.size() > 0)
    	           {
    	               *(extPtr) = const_cast<scalar*>(exec.template copyFromFoam<scalar>
                       (
                           offSize*nnzExt,
				    	   reinterpret_cast<const scalar*>(foamExtVals.cdata())
				       ));
    	           }
               }, csrMatExec_);
}

template<class sourceType, class blockType>
void Foam::coupledCsrMatrix::blkConsInit
(
	const coupledMatrix& cpldMatrix,
    const blockFvMatrix<sourceType,blockType>& matrix,
	const label nCells,
	const label nIntFaces,
	const label nnzExt,
	const label offSize,
    scalar** diagPtr,
    scalar** upperPtr,
    scalar** lowerPtr,
    scalar** extPtr
)
{
    const label nC = pTraits<blockType>::nComponents;
    List<List<blockType>> diagLst(gpuWorldSize_);
    if (matrix.hasDiag())
    {
        diagLst[myGpuWorldRank_] = matrix.diag();
    }
    Pstream::gatherList(diagLst, UPstream::msgType(), gpuWorld_);

    List<List<blockType>> upperLst(gpuWorldSize_);
    if (matrix.hasUpper())
    {
        upperLst[myGpuWorldRank_] = matrix.upper();
    }
    Pstream::gatherList(upperLst, UPstream::msgType(), gpuWorld_);

    List<List<blockType>> lowerLst(gpuWorldSize_);
    if (matrix.hasLower())
    {
        lowerLst[myGpuWorldRank_] = matrix.lower();
    }
    Pstream::gatherList(lowerLst, UPstream::msgType(), gpuWorld_);

    Field<blockType> foamExtVals(nnzExt,Foam::Zero);
    List<List<blockType>> extValLst(gpuWorldSize_);
    label localNnz = 0;
    forAll(cpldMatrix.mesh().boundary(), patchi)
    {
	   if (cpldMatrix.mesh().boundary()[patchi].coupled() && matrix.interfacesUpper().set(patchi))
	   {
           //- Processor-local values
           const Field<blockType>& bCoeffs = matrix.interfacesUpper()[patchi];
           const label len = bCoeffs.size();

           SubList<blockType>(foamExtVals, len, localNnz) = bCoeffs;
           localNnz += len;
        }
    }
    if (localNnz != nnzExt && !matrix.diagonal())
    {
    	Info << "Warning: dimension mismatch in interface consolidation" << endl;
    }
    extValLst[myGpuWorldRank_] = foamExtVals;
    Pstream::gatherList(extValLst, UPstream::msgType(), gpuWorld_);

    if(gpuProc_)
    {
        std::visit([this, &diagPtr](const auto& exec)
                    { *(diagPtr) = exec.template alloc<scalar>(nC*this->nConsRows_); },
                    csrMatExec_);
        std::visit([this, &diagLst, &diagPtr](const auto& exec)
                    { exec.template concatenate<blockType>(this->nConsRows_, diagLst, *(diagPtr)); },
                    coupledCsrMatExec_);

        std::visit([this, &upperPtr](const auto& exec)
                    { *(upperPtr) = exec.template alloc<scalar>(nC*this->nConsIntFaces_); },
                    csrMatExec_);
        std::visit([this, &upperLst, &upperPtr](const auto& exec)
                    { exec.template concatenate<blockType>(this->nConsIntFaces_, upperLst, *(upperPtr)); },
                    coupledCsrMatExec_);

        std::visit([this, &lowerPtr](const auto& exec)
                    { *(lowerPtr) = exec.template alloc<scalar>(nC*this->nConsIntFaces_); },
                    csrMatExec_);
        std::visit([this, &lowerLst, &lowerPtr](const auto& exec)
                    { exec.template concatenate<blockType>(this->nConsIntFaces_, lowerLst, *(lowerPtr)); },
                    coupledCsrMatExec_);

        std::visit([this, &extPtr](const auto& exec)
                    { *(extPtr) = exec.template alloc<scalar>(nC*this->nConsExtNz_); },
                    csrMatExec_);
        std::visit([this, &extValLst, &extPtr](const auto& exec)
                    { exec.template concatenate<blockType>(this->nConsExtNz_, extValLst, *(extPtr)); },
                    coupledCsrMatExec_);
    }


    Pstream::barrier(gpuWorld_);
}

void Foam::coupledCsrMatrix::fillSource
(
    const coupledMatrix& matrix,
    const PtrList<scalarField>& sSource,
	const PtrList<vectorField>& vSource
)
{
    const label nScalar = matrix.nScal();
    const label nVector = matrix.nVect();
    const label nCells = matrix.mesh().nCells();
    //scalar* sourcePtr = nullptr;
    //label consDispl = 0;
    if(!isConsolidated())
    {
    //}
    //else
    //{
        std::visit([this,nCells](const auto& exec)
               { this->rhsCons_ = exec.template allocZero<scalar>(nBlocks_*nCells); },
               this->csrMatExec_);
    }

    forN(nScalar,j)
    {
        if(isConsolidated())
        {
            List<scalarField> sourceLst(gpuWorldSize_);
            sourceLst[myGpuWorldRank_] = sSource[j];
            Pstream::gatherList(sourceLst, UPstream::msgType(), gpuWorld_);

            if(gpuProc_)
            {
            	for (label nCons = 0; nCons < gpuWorldSize_; nCons++)
            	{
                     label consDispl = rowsConsDispPtr_->cdata()[nCons];
    	            fillField<scalar>
    	            (
                        nCells,
		            	j,
		            	sourceLst[nCons],
		            	&rhsCons_[consDispl*nBlocks_]
    	            );
            	}
            }
        }
        else
        {
    	    fillField<scalar>
    	    (
                nCells,
		    	j,
		    	sSource[j],
		    	rhsCons_
    	    );
        }
    }

    forN(nVector,j)
    {
        if(isConsolidated())
        {
            List<vectorField> sourceLst(gpuWorldSize_);
            sourceLst[myGpuWorldRank_] = vSource[j];
            Pstream::gatherList(sourceLst, UPstream::msgType(), gpuWorld_);

            if(gpuProc_)
            {
            	for (label nCons = 0; nCons < gpuWorldSize_; nCons++)
            	{
                    label consDispl = rowsConsDispPtr_->cdata()[nCons];
    	            fillField<vector>
    	            (
                        nCells,
		            	nScalar + 3*j,
		            	sourceLst[nCons],
		            	&rhsCons_[consDispl*nBlocks_]
    	            );
            	}
            }
        }
        else
        {
    	       fillField<vector>
    	       (
                nCells,
		       	nScalar + 3*j,
		       	vSource[j],
		       	rhsCons_
    	       );
        }
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

    if(!isConsolidated())
    {
    //    consDispl = rowsConsDispPtr_->cdata()[myGpuWorldRank_];
    //}
    //else
    //{
        std::visit([this,nCells](const auto& exec)
               { this->psiCons_ = exec.template allocZero<scalar>(nBlocks_*nCells); },
               this->csrMatExec_);
    }

    forN(nScalar,j)
    {
        if(isConsolidated())
        {
            List<scalarField> psiLst(gpuWorldSize_);
            psiLst[myGpuWorldRank_] = sVolField[j].primitiveField();
            Pstream::gatherList(psiLst, UPstream::msgType(), gpuWorld_);

            if(gpuProc_)
            {
            	for (label nCons = 0; nCons < gpuWorldSize_; nCons++)
            	{
                    label consDispl = rowsConsDispPtr_->cdata()[nCons];
    	            fillField<scalar>
    	            (
                        nCells,
		            	j,
		            	psiLst[nCons],
		            	&psiCons_[consDispl*nBlocks_]
    	            );
            	}
            }
        }
        else
        {
    	    fillField<scalar>
    	    (
                nCells,
		    	j,
		    	sVolField[j].primitiveField(),
		    	psiCons_
    	    );
        }
    }

    forN(nVector,j)
    {
        if(isConsolidated())
        {
            List<vectorField> psiLst(gpuWorldSize_);
            psiLst[myGpuWorldRank_] = vVolField[j].primitiveField();
            Pstream::gatherList(psiLst, UPstream::msgType(), gpuWorld_);

            if(gpuProc_)
            {
            	for (label nCons = 0; nCons < gpuWorldSize_; nCons++)
            	{
                    label consDispl = rowsConsDispPtr_->cdata()[nCons];
    	            fillField<vector>
    	            (
                        nCells,
		            	nScalar + 3*j,
		            	psiLst[nCons],
		            	&psiCons_[consDispl*nBlocks_]
    	            );
            	}
            }
        }
        else
        {
    	       fillField<vector>
    	       (
                nCells,
		       	nScalar + 3*j,
		       	vVolField[j].primitiveField(),
		       	psiCons_
    	       );
        }
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

    scalarField hostFld(nCells*nBlocks_);
    scalar* hostPtr = hostFld.data();
    scalar* devPtr = nullptr;
    if(isConsolidated())
    {
        scalarField psiLst(nConsRows_*nBlocks_);
        labelList nLocalRows(gpuWorldSize_);
        forAll(nLocalRows, proci)
        {
            nLocalRows[proci] = rowsConsDispPtr_->cdata()[proci+1] - rowsConsDispPtr_->cdata()[proci];
        }

        if(gpuProc_)
        {
            scalar* psiPtr = psiLst.data();
            scalar** psiPPtr = &psiPtr;
            std::visit([this, &psiPPtr](const auto& exec)
                { exec.template copyToFoam<scalar>(this->nConsRows_*this->nBlocks_, this->psiCons_, psiPPtr); },
                csrMatExec_);
        }

        UPstream::scatter(psiLst.cdata(), nLocalRows*nBlocks_, *rowsConsDispPtr_*nBlocks_, hostPtr, nBlocks_*nLocalRows[myGpuWorldRank_], gpuWorld_);
    }
    else
    {
    	devPtr = psiCons_;
        std::visit([this, &hostPtr, &devPtr, nCells](const auto& exec)
                   { exec.template copyToFoam<scalar>(nCells*this->nBlocks_, devPtr, &hostPtr);},
                   csrMatExec_);
    }


    forN(nScalar,j)
    {
        for(label i=0; i<nCells;i++)
        {
            sVolField[j].primitiveFieldRef()[i] = hostFld[i*nBlocks_ + j];
        }
    }
    forN(nVector,j)
    {
        for(label i=0; i<nCells;i++)
        {
            vVolField[j].primitiveFieldRef()[i].x() = hostFld[i*nBlocks_ + nScalar + 3*j];
            vVolField[j].primitiveFieldRef()[i].y() = hostFld[i*nBlocks_ + nScalar + 3*j + 1];
            vVolField[j].primitiveFieldRef()[i].z() = hostFld[i*nBlocks_ + nScalar + 3*j + 2];
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
