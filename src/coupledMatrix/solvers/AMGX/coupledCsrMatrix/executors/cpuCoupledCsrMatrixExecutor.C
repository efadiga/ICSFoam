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

#include <cmath>
#include <bits/stdc++.h>
#include "cpuCoupledCsrMatrixExecutor.H"
#include "vector.H"
// * * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * //

void Foam::cpuCoupledCsrMatrixExecutor::initializeAndApplyValue
(
    const label nBlocks,
    const label nCells,
    const label nIntFaces,
	const label nnzExt,
    const label nOffsets,
    const label * const offsets,
    const label  * const ldu2csr,
    const scalar * const diag,
    const scalar * const upper,
    const scalar * const lower,
    const scalar * const ext,
          scalar * values
) const
{
    for(label i=0; i<nCells; ++i)
    {
        values[i] = diag[i];
    }

    for(label i=0; i<nIntFaces; ++i)
    {
        values[nCells + i] = upper[i];
        values[nCells + nIntFaces + i] = lower[i];
    }
}

template<class Type>
void Foam::cpuCoupledCsrMatrixExecutor::fillField
(
	const label nCells,
	const label position,
	const label nBlocks,
	const Field<Type>& input,
		  scalar * output
) const
{
    const label nComps = pTraits<Type>::nComponents;

}

template<class Type>
void Foam::cpuCoupledCsrMatrixExecutor::concatenate
(
    label globSize,
    List<List<Type>> lst,
    scalar * ptr
) const
{
	NotImplemented;
    //label newStart = 0;
    //label size;
    //const label nC = pTraits<Type>::nComponents;

    //for(label i=0; i<lst.size(); ++i)
    //{
    //    size = lst[i].size();
    //    label err = CHECK_CUDA_ERROR(
    //                    cudaMemcpy
	//					(
    //                        &ptr[newStart],
	//						reinterpret_cast<const scalar*>(lst[i].cdata()),
	//						(size_t) size*sizeof(scalar)*nC,
	//						cudaMemcpyHostToDevice
	//					)
    //                );
    //    if (err != 0)
    //    {
    //        FatalErrorInFunction << "ERROR: cudaMemcpy returned " << err << abort(FatalError);
    //    }
    //    newStart += size;
    //    if(newStart > globSize)
    //    {
    //        FatalErrorInFunction << "Concatenate size mismatch" << nl;
    //    }
    //}
}

//void Foam::cpuCoupledCsrMatrixExecutor::initializeValueExt
//(
//    const label   nCells,
//    const label   nIntFaces,
//    const label   nnzExt,
//    const scalar * const diag,
//    const scalar * const upper,
//    const scalar * const lower,
//    const scalar * const extValue,
//          scalar * valuesTmp
//) const
//{
//    // Initialize valuesTmp = [(diag), (upper), (lower), (extValues)]
//
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
//    for(label i=0; i<nnzExt; ++i)
//    {
//        valuesTmp[nCells + 2*nIntFaces + i] = extValue[i];
//    }
//}
//
//void Foam::cpuCoupledCsrMatrixExecutor::applyValuePermutation
//(
//    const label    totNnz,
//    const label *  const ldu2csr,
//    const scalar * const valuesTmp,
//          scalar * values,
//    const label    nBlocks
//) const
//{
//    label blockLen = nBlocks * nBlocks;
//
//    for(label i=0; i<totNnz; ++i)
//    {
//        for(label j=0; j<blockLen; ++j)
//        {
//            values[ldu2csr[i]*blockLen + j] = valuesTmp[i*blockLen + j];
//        }
//    }
//}
//

// ************************************************************************* //
#define makecpuCoupledCsrMatrixExecutor(Type)                                 \
    template void Foam::cpuCoupledCsrMatrixExecutor::fillField<Type>          \
    (                                                                         \
	const Foam::label nCells,                                                       \
	const Foam::label position,                                                     \
	const Foam::label nBlocks,                                                      \
	const Foam::Field<Type>& input,                                               \
		  Foam::scalar * output                                                     \
    ) const;

makecpuCoupledCsrMatrixExecutor(Foam::scalar)
makecpuCoupledCsrMatrixExecutor(Foam::vector)
template void Foam::cpuCoupledCsrMatrixExecutor::concatenate<Foam::scalar>
(
    Foam::label globSize,
    Foam::List<Foam::List<Foam::scalar>> lst,
    Foam::scalar * ptr
) const;

template void Foam::cpuCoupledCsrMatrixExecutor::concatenate<Foam::vector>
(
    Foam::label globSize,
    Foam::List<Foam::List<Foam::vector>> lst,
    Foam::scalar * ptr
) const;

template void Foam::cpuCoupledCsrMatrixExecutor::concatenate<Foam::tensor>
(
    Foam::label globSize,
    Foam::List<Foam::List<Foam::tensor>> lst,
    Foam::scalar * ptr
) const;
