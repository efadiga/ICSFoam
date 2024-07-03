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

#include "coupledCsrMatrix.C"

// * * * * * * * * * * * * * * * * * Kernels * * * * * * * * * * * * * * * * //

void Foam::coupledCsrMatrix::initializeAndApplyValue
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
)
{
    std::visit
	([
	    nBlocks,
        nCells,
        nIntFaces,
		nOffsets,
		nnzExt,
		&offsets,
		&ldu2csr,
        &diag,
        &upper,
        &lower,
		&ext,
        &values
	 ]
	 (const auto& exec){ exec.initializeAndApplyValue
                           (
                        	   nBlocks,
                               nCells,
                               nIntFaces,
							   nnzExt,
		                       nOffsets,
		                       offsets,
							   ldu2csr,
                               diag,
                               upper,
                               lower,
							   ext,
                               values
					       );
                       },
     coupledCsrMatExec_);
}

template<class Type>
void Foam::coupledCsrMatrix::fillField
(
    const label nCells,
    const label position,
    const Field<Type>& input,
    scalar * output
) const
{
    std::visit
	([
	    this,
        nCells,
        position,
		input,
		&output
	 ]
	 (const auto& exec){ exec.template fillField<Type>
                           (
                               nCells,
							   position,
							   this->nBlocks_,
                               input,
							   output
					       );
                       },
     coupledCsrMatExec_);
}


//inline void Foam::csrMatrix::initializeValueExt
//(
//    const label nCells,
//    const label nIntFaces,
//    const label nnzExt,
//    const scalar * const diag,
//    const scalar * const upper,
//    const scalar * const lower,
//    const scalar * const extValue,
//          scalar * valuesTmp
//)
//{
//    std::visit
//	([
//        nCells,
//        nIntFaces,
//        nnzExt,
//        &diag,
//        &upper,
//        &lower,
//        &extValue,
//        &valuesTmp
//	 ]
//	 (const auto& exec){ exec.initializeValueExt
//                           (
//                                nCells,
//                                nIntFaces,
//                                nnzExt,
//                                diag,
//                                upper,
//                                lower,
//                                extValue,
//                                valuesTmp
//					       );
//                       },
//     coupledCsrMatExec_);
//}
//
//
//inline void Foam::csrMatrix::applyValuePermutation
//(
//    const label totNnz,
//    const label  * const ldu2csr,
//    const scalar * const valuesTmp,
//          scalar * values
//)
//{
//    std::visit
//	([
//        totNnz,
//        &ldu2csr,
//        &valuesTmp,
//        &values
//	 ]
//	 (const auto& exec){ exec.applyValuePermutation
//                            (
//                                totNnz,
//                                ldu2csr,
//                                valuesTmp,
//                                values
//					        );
//                       },
//     coupledCsrMatExec_);
//}

// * * * * * * * * * * * *  Public Member Functions * * * * * * * * * * * *  //

// * * * * * * * * * * * * * Explicit instantiations  * * * * * * * * * * *  //

// ************************************************************************* //
