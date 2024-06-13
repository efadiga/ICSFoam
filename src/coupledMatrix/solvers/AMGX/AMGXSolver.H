/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::AMGXSolver

Description
    Interface to the external package AmgX (https://github.com/NVIDIA/AMGX).

SourceFiles
    AMGXSolver.C

\*---------------------------------------------------------------------------*/

#ifndef AMGXSolver_H
#define AMGXSolver_H

#include "coupledMatrix.H"
#include "csrMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class AMGXSolver Declaration
\*---------------------------------------------------------------------------*/

class AMGXSolver
:
    public coupledMatrix::solver
{
    // Private Member Functions

        //- No copy construct
        AMGXSolver(const AMGXSolver&) = delete;

        //- No copy assignment
        void operator=(const AMGXSolver&) = delete;

        //- Compute matrix addressing and value permutation for AmgX use
        /*void buildAndApplyMatrixPermutation
        (
            deviceCsrMatrix* csrMatrix,
            label& nRowsGlobal
        ) const;*/

        //- Compute matrix value permutation for AmgX use
        /*void applyMatrixPermutation
        (
            deviceCsrMatrix* csrMatrix,
            label& nRowsGlobal
        ) const;*/


public:

    //- Runtime type information
    TypeName("AmgX");


    // Constructors

        //- Construct from matrix components and solver controls
        AMGXSolver
        (
            const dictionary& dict,
            const coupledMatrix& matrix
        );


    //- Destructor
        virtual ~AMGXSolver() = default;

        // Main solve functions:

        //- Returns solver performance
        virtual residualsIO solve
        (
            PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
            const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource
        ) const;

        //- Returns solver performance
        virtual residualsIO solveDelta
        (
            PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
            const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource
        ) const;

        //- Returns solver performance
        virtual residualsIO solveDelta
        (
            PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
            const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource,
			PtrList<volScalarField>& dsW, PtrList<volVectorField>& dvW
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
