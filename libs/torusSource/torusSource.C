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

\*---------------------------------------------------------------------------*/

#include "torusSource.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(torusSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        torusSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::torusSource::writeProps
(
    const scalar gradP
) const
{
    // Only write on output time
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::torusSource::torusSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(sourceName, modelType, dict, mesh),
    Ubar_(readScalar(coeffs_.lookup("Ubar"))),
    gradP0_(0.0),
    dGradP_(0.0),
    flowDir_(NULL),	//initialized below ...
    relaxation_(coeffs_.lookupOrDefault<scalar>("relaxation", 1.0)),
    rAPtr_(NULL)
{
    coeffs_.lookup("fields") >> fieldNames_;
    outputFlowdir_ = coeffs_.lookupOrDefault<bool>("writeForceDirVecs", true);
    axisOfR_ = coeffs_.lookupOrDefault<vector>("axisOfRotation", vector(0,0,0));
	
    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);

    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timePath()/"uniform"/(name_ + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("gradient") >> gradP0_;
    }

    Info<< "    Initial pressure gradient = " << gradP0_ << nl << endl;
    
    
    flowDir_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "dir",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
			vector::zero
        )
     );
     
    
    computeDirVecs(mesh);
    
    if (outputFlowdir_){
    	flowDir_->write();
    }
   
     
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::torusSource::computeDirVecs(const fvMesh& mesh){

	//cell centers
	const volVectorField& ccs = mesh.C();

	vectorField& dirvecs = flowDir_->primitiveFieldRef(); 

	
	int component1 = 100; int component2 = 100;
	
	if (axisOfR_ == vector(0, 0, 1)){
		component1 = 0; //x
		component2 = 1; //y
	
	}
	
	else if (axisOfR_ == vector(0, 1, 0)){
		component1 = 0; //x
		component2 = 2; //z
	
	}
	
	else if (axisOfR_ == vector(1, 0, 0)){
		component1 = 1; //x
		component2 = 2; //z
	
	}
	
	else {
		Info << endl << endl << endl;
		Info << "Axis of rotation for torusSource has to be one of the cartesian axes" 
		<< endl;
		NotImplemented;
		
		
	
	}
	
	
	forAll(ccs, i)
    {
    	// Make sure that this is what you want to do. 
    	
    	
        vector cc = ccs[i];
        scalar R1 = cc.component(component1);
        scalar R2 = cc.component(component2);  
        
        vector new_vec(0,0,0);
        new_vec.component(component1) = R2;
        new_vec.component(component2) = -R1;
        
        new_vec /=mag(new_vec);
        
        dirvecs[i] = new_vec;
        
        
    }


}



Foam::scalar Foam::fv::torusSource::magUbarAve
(
    const volVectorField& U
) const
{
    scalar magUbarAve = 0.0;
	//vectorField& dirvecs = flowDir_->primitiveFieldRef();
	const vectorField& dirvecs = flowDir_;


    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label celli = cells_[i];
        scalar volCell = cv[celli];
        //magUbarAve += (flowDir_ & U[celli])*volCell;
        magUbarAve += (dirvecs[i] & U[celli])*volCell;
        
    }

    reduce(magUbarAve, sumOp<scalar>());

    magUbarAve /= V_;

    return magUbarAve;
}


void Foam::fv::torusSource::correct(volVectorField& U)
{
    const scalarField& rAU = rAPtr_();

    // Integrate flow variables over cell set
    scalar rAUave = 0.0;
    const scalarField& cv = mesh_.V();
    
    
    
    forAll(cells_, i)
    {
        label celli = cells_[i];
        scalar volCell = cv[celli];
        rAUave += rAU[celli]*volCell;
    }

    // Collect across all processors
    reduce(rAUave, sumOp<scalar>());

    // Volume averages
    rAUave /= V_;

    scalar magUbarAve = this->magUbarAve(U);

    // Calculate the pressure gradient increment needed to adjust the average
    // flow-rate to the desired value
    //dGradP_ = relaxation_*(mag(Ubar_) - magUbarAve)/rAUave;
	dGradP_ = relaxation_*(Ubar_ - magUbarAve)/rAUave;
	
	const vectorField& dirvecs = flowDir_;
    // Apply correction to velocity field
    forAll(cells_, i)
    {
        label celli = cells_[i];
        U[celli] += dirvecs[i]*rAU[celli]*dGradP_;
    }

    scalar gradP = gradP0_ + dGradP_;

    Info<< "Pressure gradient source: uncorrected Ubar = " << magUbarAve
        << ", pressure gradient = " << gradP << endl;

	

    writeProps(gradP);
}


void Foam::fv::torusSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{

	
	
	
    DimensionedField<vector, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", eqn.dimensions()/dimVolume, Zero)
    );

    scalar gradP = gradP0_ + dGradP_;


	const vectorField& dirvecs = flowDir_;
	
	/*forAll(cells_, i){
	
		UIndirectList<vector>(Su, cells_)[i] = gradP * dirvecs[i];
	}*/
	
	UIndirectList<vector>(Su, cells_) = gradP * dirvecs;

    //UIndirectList<vector>(Su, cells_) = flowDir_*gradP;
    

    eqn += Su;
}


void Foam::fv::torusSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    this->addSup(eqn, fieldi);
}


void Foam::fv::torusSource::constrain
(
    fvMatrix<vector>& eqn,
    const label
)
{
    if (rAPtr_.empty())
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + ":rA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                1.0/eqn.A()
            )
        );
    }
    else
    {
        rAPtr_() = 1.0/eqn.A();
    }

    gradP0_ += dGradP_;
    dGradP_ = 0.0;
}


// ************************************************************************* //
