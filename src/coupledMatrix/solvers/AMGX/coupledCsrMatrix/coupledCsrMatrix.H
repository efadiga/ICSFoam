/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
    Copyright (C) 2023      Cineca
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
    Foam::coupledCsrMatrix

Description
     csrMatrix class for coupled problems.

SourceFiles
    coupledCsrMatrix.H
    coupledCsrMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coupledCsrMatrix_H
#define Foam_coupledCsrMatrix_H

#include "Enum.H"
#include "scalar.H"
#include "coupledMatrix.H"
#include "coupledCsrMatrixExecutor.H"
#include "blockFvMatrix.H"
#include "csrMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class coupledCsrMatrix;


/*---------------------------------------------------------------------------*\
                           Class coupledCsrMatrix Declaration
\*---------------------------------------------------------------------------*/

class coupledCsrMatrix
:
    public csrMatrix
{
    // Private Data

        scalar * source_;
        scalar * variables_;
        //- Data executor
        coupledCsrMatrixExecutor coupledCsrMatExec_;

    // Private methods
           
        //void initializeValue
        //(
        //    const int   nCells,
        //    const int   nIntFaces,
        //    const double * const diag,
        //    const double * const upper,
        //    const double * const lower,
        //          double * valuesTmp
        //);

        //void initializeValueExt
        //(
        //    const int      nCells,
        //    const int      nIntFaces,
        //    const int      nnzExt,
        //    const double * const diag,
        //    const double * const upper,
        //    const double * const lower,
        //    const double * const extValue,
        //          double * valuesTmp
        //);

public:

    //- Public Types

    //- Static Data

        // Declare name of the class and its debug switch
        ClassName("coupledCsrMatrix");

    //- Constructors

        //- Construct default.
        //  The coefficients are initialized empty.
        coupledCsrMatrix(word mode);

    //- Destructor
        ~coupledCsrMatrix() = default;


    //- Member Functions

            //- Apply permutation from LDU to CSR
            void applyPermutation(const coupledMatrix& matrix);

            //- Helper function that copies LDU pointers
            template<class sourceType, class blockType>
            void copyLDUPtrs
			(
				const blockFvMatrix<sourceType,blockType>& matrix,
			    const label nCells,
				const label nIntFaces,
                const scalar** diagPtr,
                const scalar** upperPtr,
                const scalar** lowerPtr
			);


            //- Fills source array with re-ordered values (block AoS)
            void fillSource
			(
				const coupledMatrix& matrix,
			    const PtrList<scalarField>& sSource,
				const PtrList<vectorField>& vSource
			);

            //- Fills source array with re-ordered values (block AoS)
            void fillVariables
			(
				const coupledMatrix& matrix,
			    const PtrList<volScalarField>& sVolField,
				const PtrList<volVectorField>& vVolField
			);

            //- Fills source array with re-ordered values (block AoS)
            void transferVariables
			(
				const coupledMatrix& matrix,
			    PtrList<volScalarField>& sVolField,
				PtrList<volVectorField>& vVolField
			);

            void initializeAndApplyValue
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
            );

            template<class Type>
            void fillField
            (
                const label nCells,
                const label position,
                const Field<Type>& input,
                scalar * output
            ) const;
            ////- Apply permutation from LDU to CSR considering the interface values
            //virtual void applyPermutation
            //(
            //    const lduMatrix& lduMatrix,
            //    const label diagIndexGlobal,
            //    const label lowOffGlobal,
            //    const label uppOffGlobal,
            //    const labelList& extRows,
            //    const labelList& extCols,
            //    const scalarField& extVals
            //);

            ////- Apply permutation from LDU to CSR considering the interface values
            //virtual void applyPermutation
            //(
            //    const lduMatrix& lduMatrix,
            //    const scalarField& extVals
            //);

            scalar* source()
            {
            	return source_;
            };

            void allocSource(const coupledMatrix& matrix)
            {
                std::visit([this, matrix](const auto& exec)
                           { this->source_ = exec.template alloc<scalar>(matrix.mesh().nCells()*nBlocks_); },
                           csrMatExec_);
            }
            void clearSource()
            {
                std::visit([this](const auto& exec)
                       {exec.template clear<scalar>(this->source_); }, csrMatExec_);
            }

            scalar* variables()
            {
            	return variables_;
            };

            void allocVariables(const coupledMatrix& matrix)
            {
                std::visit([this, matrix](const auto& exec)
                           { this->variables_ = exec.template alloc<scalar>(matrix.mesh().nCells()*nBlocks_); },
                           csrMatExec_);
            }
            void clearVariables()
            {
                std::visit([this](const auto& exec)
                       {exec.template clear<scalar>(this->variables_); }, csrMatExec_);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
