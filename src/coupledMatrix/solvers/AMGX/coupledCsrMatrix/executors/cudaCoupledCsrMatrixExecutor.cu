/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

// ************************************************************************* //

#include "cudaCoupledCsrMatrixExecutor.H"

#include "scalar.H"
#include "csrMatrix.H"
#include "global.cuh"
#include <cub/cub.cuh>

// * * * * * * * * * * * * * * * * CUDA Kernels  * * * * * * * * * * * * * * //

__global__
void cudaInitializeValueD
(
	const int   nBlocks,
    const int   nCells,
	const int   nOffsets,
	const int * const offsets,
	const int * const ldu2csr,
    const double * const diag,
          double * values
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < nCells)
    {
    	for(int j=0; j<nOffsets;j++)
    	{
            values[ldu2csr[i]*nBlocks*nBlocks + offsets[j]] = diag[i*nOffsets + j];
    	}
    }
}

__global__
void cudaInitializeValueUL
(
	const int   nBlocks,
    const int   nCells,
    const int   nIntFaces,
	const int   nOffsets,
	const int * const offsets,
	const int * const ldu2csr,
    const double * const upper,
    const double * const lower,
          double * values
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < nIntFaces)
    {
    	for(int j=0; j<nOffsets;j++)
    	{
            values[nBlocks*nBlocks*ldu2csr[nCells + i]+offsets[j]] =
                upper[i*nOffsets + j];
            values[nBlocks*nBlocks*ldu2csr[nCells + nIntFaces + i]+offsets[j]] =
                lower[i*nOffsets + j];
    	}
    }
}

__global__
void cudaInitializeValueExt
(
    const int   nCells,
    const int   nIntFaces,
    const int   nnzExt,
    const double * const extValues,
          double * valuesTmp
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < nnzExt)
    {
        valuesTmp[nCells + 2*nIntFaces + i] = extValues[i];
    }
}

__global__
void cudaApplyValuePermutation 
(
    const int      length,
    const int      blockLen,
    const int    * const permArray,
    const double * const srcArray,
          double * dstArray
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < length)
    {
        dstArray[permArray[i]] = srcArray[i];
    }
} 
//NOTA: this function (when csrAdressing will be joined back to csrMatrix) will 
//      become e template on the array type to be used both for adressing and 
//      values permutaiton

template <int nComps>
__global__
void cudaFillField
(
	const int nCells,
	const int position,
	const int nBlocks,
	const double * const input,
		  double * output
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < nCells)
    {
    	for(int j=0; j<nComps; j++)
    	{
            output[i*nBlocks + position + j] = input[nComps*i + j];
    	}
    }
}

// * * * * * * * * * * * * * *  Wrapper functions * * * * * * * * * * * * * * //

template<class Type>
void Foam::cudaCoupledCsrMatrixExecutor::fillField
(
	const label nCells,
	const label position,
	const label nBlocks,
	const Field<Type>& input,
		  scalar * output
) const
{
    const label nComps = pTraits<Type>::nComponents;

	void* inPtr;
    label err = CHECK_CUDA_ERROR(cudaMalloc((void**)&inPtr, (size_t) nCells*nBlocks*sizeof(scalar)));
    if (err != 0)
    {
        FatalErrorInFunction << "ERROR: cudaMalloc returned " << err << abort(FatalError);
    }

    err = CHECK_CUDA_ERROR(cudaMemcpy(inPtr, reinterpret_cast<const scalar*>(input.cdata()),
    		                             (size_t) nCells*nBlocks*sizeof(scalar), cudaMemcpyHostToDevice));
    if (err != 0)
    {
        FatalErrorInFunction << "ERROR: cudaMemcpy returned " << err << abort(FatalError);
    }

    int numBlocks = (nCells + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK;
    cudaFillField<nComps><<<numBlocks, NUM_THREADS_PER_BLOCK>>>
    (
    	nCells,
		position,
		nBlocks,
		static_cast<const scalar*>(inPtr),
		output
    );

    cudaDeviceSynchronize();
    CHECK_LAST_CUDA_ERROR();
    err = CHECK_CUDA_ERROR(cudaFree(inPtr));
    CHECK_LAST_CUDA_ERROR();
}



void Foam::cudaCoupledCsrMatrixExecutor::initializeAndApplyValue
(
	const label nBlocks,
    const label nCells,
    const label nIntFaces,
    const label nOffsets,
    const label * const offsets,
    const label  * const ldu2csr,
    const scalar * const diag,
    const scalar * const upper,
    const scalar * const lower,
          scalar * values
) const
{
    label numBlocks;
    
    // Initialize valuesTmp = [(diag), (upper), (lower)]
    if (diag)
    {
        numBlocks = (nCells + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK;
        cudaInitializeValueD<<<numBlocks, NUM_THREADS_PER_BLOCK>>>
        (
            nBlocks,
            nCells,
            nOffsets,
	    	offsets,
	    	ldu2csr,
            diag,
            values
        );
    }
    if (upper)
    {
        numBlocks = (nIntFaces + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK;
        cudaInitializeValueUL<<<numBlocks, NUM_THREADS_PER_BLOCK>>>
        (
            nBlocks,
            nCells,
            nIntFaces,
	    	nOffsets,
	    	offsets,
	    	ldu2csr,
            upper,
            lower,
            values
        );
    }
    cudaDeviceSynchronize();

    CHECK_LAST_CUDA_ERROR();
    return;
}


//void Foam::cudaCoupledCsrMatrixExecutor::initializeValueExt
//(
//    const label    nCells,
//    const label    nIntFaces,
//    const label    nnzExt,
//    const scalar * const diag,
//    const scalar * const upper,
//    const scalar * const lower,
//    const scalar * const extValue,
//          scalar * valuesTmp
//) const
//{
//    // Initialize valuesTmp = [(diag), (upper), (lower), (extValues)]
//    initializeValue
//    (
//        nCells,
//        nIntFaces,
//        diag,
//        upper,
//        lower,
//        valuesTmp
//    );
//
//    int numBlocks = (nnzExt + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK;
//    cudaInitializeValueExt<<<numBlocks, NUM_THREADS_PER_BLOCK>>>
//    (
//        nCells,
//        nIntFaces,
//        nnzExt,
//        extValue,
//        valuesTmp
//    );
//
//    cudaDeviceSynchronize();
//
//    CHECK_LAST_CUDA_ERROR();
//    return;
//}
//
//
//void Foam::cudaCoupledCsrMatrixExecutor::applyValuePermutation
//(
//    const label    totNnz,
//    const label  * const ldu2csr,
//    const scalar * const valuesTmp,
//          scalar * values,
//    const label    nBlocks
//) const
//{
//    int blockLen = nBlocks*nBlocks;
//    int numBlocks = (totNnz + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK;
//    cudaApplyValuePermutation<<<numBlocks, NUM_THREADS_PER_BLOCK>>>
//    (
//        totNnz,
//        blockLen,
//        ldu2csr,
//        valuesTmp,
//        values
//    );
//
//    cudaDeviceSynchronize();
//
//    CHECK_LAST_CUDA_ERROR();
//    return;
//}

// * * * * * * * * * * * * * Explicit instantiations  * * * * * * * * * * *  //

#define makecudaCoupledCsrMatrixExecutor(Type)                                 \
    template void Foam::cudaCoupledCsrMatrixExecutor::fillField<Type>          \
    (                                                                         \
	const Foam::label nCells,                                                       \
	const Foam::label position,                                                     \
	const Foam::label nBlocks,                                                      \
	const Foam::Field<Type>& input,                                               \
		  Foam::scalar * output                                                     \
    ) const;

makecudaCoupledCsrMatrixExecutor(Foam::scalar)
makecudaCoupledCsrMatrixExecutor(Foam::vector)
// ************************************************************************* //
