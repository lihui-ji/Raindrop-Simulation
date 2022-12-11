/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::solidBodyMotionFunctions::dynamicMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Constant velocity displacement.

SourceFiles
    dynamicMotion.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotion_H
#define dynamicMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class dynamicMotion Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Linear velocity
        vector velocity_;
        //- Whether to move along x axis
        bool ifMoveX_;
        //- Whether to move along y axis
        bool ifMoveY_;
        //- Whether to move along z axis
        bool ifMoveZ_;

    // Private Member Functions

        //- Disallow copy construct
        dynamicMotion(const dynamicMotion&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicMotion&);


public:

    //- Runtime type information
    TypeName("dynamicMotion");


    // Constructors

        //- Construct from components
        dynamicMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new dynamicMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~dynamicMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

	    virtual septernion transformation(vector displacement);

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
