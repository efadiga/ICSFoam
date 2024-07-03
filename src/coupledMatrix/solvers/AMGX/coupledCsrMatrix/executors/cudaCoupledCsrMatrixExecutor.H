/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
    Copyright (C) 2024 Cineca
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


Class
    cudaCoupledCsrMatrixExecutor

Description
    cuda visitor of the csrMatrix class

SourceFiles
    cudaCoupledCsrMatrixExecutor.H
    cudaCoupledCsrMatrixExecutor.cu

\*---------------------------------------------------------------------------*/

#ifndef cuda_CoupledCsrMatrix_Executor_H
#define cuda_CoupledCsrMatrix_Executor_H

#include "scalar.H"
#include "Field.H"

namespace Foam
{

// Forward Declarations
class coupledCsrMatrix;

class cudaCoupledCsrMatrixExecutor
{
public:

    cudaCoupledCsrMatrixExecutor() {};

    ~cudaCoupledCsrMatrixExecutor() {};

    void initializeAndApplyValue
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
    ) const;


    template<class Type>
    void fillField
	(
	    const label nCells,
		const label position,
		const label nBlocks,
		const Field<Type>& input,
		      scalar * output
	) const;

    template<class Type>
    void concatenate
    (
        label globSize,
        List<List<Type>> lst,
        scalar * ptr
    ) const;
    //void initializeValueExt
    //(
    //    const label    nCells,
    //    const label    nIntFaces,
    //    const label    nnzExt,
    //    const scalar * const diag,
    //    const scalar * const upper,
    //    const scalar * const lower,
    //    const scalar * const extValue,
    //          scalar * valuesTmp
    //) const;

    //void applyValuePermutation
    //(
    //    const label    totNnz,
    //    const label  * const ldu2csr,
    //    const scalar * const valuesTmp,
    //          scalar * values,
    //    const label    nBlocks = 1
    //) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

