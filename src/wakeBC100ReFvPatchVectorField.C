/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

//#include <torch/torch.h>
//#include <torch/script.h>
#include "wakeBC100ReFvPatchVectorField.H"

#include <iostream>
#include <fstream>
#include <string>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wakeBC100ReFvPatchVectorField::wakeBC100ReFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    // NOTE: call the default constructor to make sure everything gets initialised properly
    fixedValueFvPatchVectorField(p, iF),
    flowSpeed_(0.),
	locationY_(0.),
	locationStreamwise_(0.),
	D_(0.),
	streamwise_("x"),
	spanwise_("y")
	//centrepoint_(vector::zero)
{}

Foam::wakeBC100ReFvPatchVectorField::wakeBC100ReFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    // NOTE: this constructor reads all of the control parameters from the boundary
    // condition definition specified in the time folder U file, imported here
    // as a dictionary reference.
    fixedValueFvPatchVectorField(p, iF),
    flowSpeed_(0.),
	locationY_(0.),
	locationStreamwise_(0.),
	D_(0.),
	streamwise_("x"),
	spanwise_("y")
	//centrepoint_(vector::zero)
{
    // NOTE: calls the = operator to assign the value to the faces held by this BC
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    // NOTE: looks up the necessary paramters
    //approximationType_ = dict.lookupOrDefault<word>("approximationType","exponential");
    dict.lookup("flowSpeed") >> flowSpeed_;
	dict.lookup("locationY") >> locationY_;
	dict.lookup("locationStreamwise") >> locationStreamwise_;
	dict.lookup("D") >> D_;
	streamwise_ = dict.lookupOrDefault<word>("streamwise","x");
	spanwise_ = dict.lookupOrDefault<word>("spanwise","y");

    // NOTE: calls the .updateCoeffs() method to calculate the inlet profile in
    // accordance with the controls which have just been read.
	updateCoeffs();
}

Foam::wakeBC100ReFvPatchVectorField::wakeBC100ReFvPatchVectorField
(
    const wakeBC100ReFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    // NOTE: this constructor, and the two subsequent ones, transfer data to the
    // instance being created from another one.
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    flowSpeed_(ptf.flowSpeed_),
	locationY_(ptf.locationY_),
	locationStreamwise_(ptf.locationStreamwise_),
	D_(ptf.D_),
	streamwise_(ptf.streamwise_),
	spanwise_(ptf.spanwise_)
{}

Foam::wakeBC100ReFvPatchVectorField::wakeBC100ReFvPatchVectorField
(
    const wakeBC100ReFvPatchVectorField& rifvpvf
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    flowSpeed_(rifvpvf.flowSpeed_),
    locationY_(rifvpvf.locationY_),
	locationStreamwise_(rifvpvf.locationStreamwise_),
    D_(rifvpvf.D_),
	streamwise_(rifvpvf.streamwise_),
	spanwise_(rifvpvf.spanwise_)
{}

Foam::wakeBC100ReFvPatchVectorField::wakeBC100ReFvPatchVectorField
(
    const wakeBC100ReFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF),
    flowSpeed_(rifvpvf.flowSpeed_),
    locationY_(rifvpvf.locationY_),
	locationStreamwise_(rifvpvf.locationStreamwise_),
    D_(rifvpvf.D_),
	streamwise_(rifvpvf.streamwise_),
	spanwise_(rifvpvf.spanwise_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wakeBC100ReFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
}


void Foam::wakeBC100ReFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

// NOTE: this is the key method which implements the actual maths for calculating
// the inlet profiles.
void Foam::wakeBC100ReFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// assign inlet velocity normal to the patch
	// by convention, patch faces point outside of the domain
	vectorField Uin = (-1.)*(patch().Sf()/patch().magSf()) * flowSpeed_;
	
	// Writing files
	std::ofstream MyFile("filename.txt");
	// Write to the file
	MyFile << "Time:" << std::endl;
	
	float current_time = this->db().time().value();
	
	if(current_time > 5.8){
		int div = current_time/5.8;
		current_time = current_time - div*5.8;
	}
	
	current_time += 10.0;
	
	
	MyFile << current_time << std::endl;
	MyFile << "y_coord:" << std::endl;
	
    // go over each face and add the BL profile for faces close to the wall
	forAll(patch().Cf(), faceI)
	{
		//scalar yOverDelta ( (1.-mag(centrepoint_ - patch().Cf()[faceI])/R_)/deltaByR_ );
		//torch::Tensor a = torch::tensor({1., 1.}, torch::requires_grad());
		
		if (spanwise_.compare("x") == 0) {
			float write_y = (patch().Cf()[faceI][0] - locationY_)/D_;
			if(write_y > 5.9){
				write_y = 5.9;
			}
			if(write_y < -6.0){
				write_y = -6.0;
			}
			MyFile << write_y << " " << locationStreamwise_ << std::endl;
		}
		else if (spanwise_.compare("y") == 0){
			float write_y = (patch().Cf()[faceI][1] - locationY_)/D_;
			if(write_y > 5.9){
				write_y = 5.9;
			}
			if(write_y < -6.0){
				write_y = -6.0;
			}
			MyFile << write_y << " " << locationStreamwise_ << std::endl;
		}
		else {
			float write_y = (patch().Cf()[faceI][2] - locationY_)/D_;
			if(write_y > 5.9){
				write_y = 5.9;
			}
			if(write_y < -6.0){
				write_y = -6.0;
			}
			MyFile << write_y << " " << locationStreamwise_ << std::endl;
		}
		
		
		//MyFile << write_y << " " << locationStreamwise_ << std::endl;
		
	}	
	
	system("./wake_BC");
	
	// Close the file
	MyFile.close();
	
	std::string myText;
	std::ifstream MyReadFile("boundaryData.txt");

	int j = 0;
	while (getline(MyReadFile, myText)) {
		
		int res = myText.compare("");

		if (res != 0) {

			std::size_t pos = myText.find(" ");
			std::string Ustr = myText.substr(0, pos);
			std::string Vstr = myText.substr(pos + 1);

			float U = std::stof(Ustr);
			float V = std::stof(Vstr);
			
			if (streamwise_.compare("x") == 0) {
				Uin[j][0] = U*flowSpeed_;
			}
			else if (streamwise_.compare("y") == 0) {
				Uin[j][1] = U*flowSpeed_;
			}
			else {
				Uin[j][3] = U*flowSpeed_;
			}
			
			if (spanwise_.compare("x") == 0) {
				Uin[j][0] = V*flowSpeed_;
			}
			else if (spanwise_.compare("y") == 0) {
				Uin[j][1] = V*flowSpeed_;
			}
			else {
				Uin[j][3] = V*flowSpeed_;
			}
			
			//Uin[j][0] = U*flowSpeed_;
			//Uin[j][1] = V*flowSpeed_;
			
			j++;
		}
	}

	MyReadFile.close();

	// set the value_ of this patch to the newly computed flow speed
    this->operator==(Uin);

    // call the base class method to make sure all the other bits and pieces get updated
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::wakeBC100ReFvPatchVectorField::write(Ostream& os) const
{
	
	fvPatchVectorField::write(os);
    os.writeKeyword("flowSpeed") << flowSpeed_ << token::END_STATEMENT << nl;
    os.writeKeyword("locationY") << locationY_ << token::END_STATEMENT << nl;
	os.writeKeyword("locationStreamwise") << locationStreamwise_ << token::END_STATEMENT << nl;
    os.writeKeyword("D") << D_ << token::END_STATEMENT << nl;
    os.writeKeyword("streamwise") << streamwise_ << token::END_STATEMENT << nl;
	os.writeKeyword("spanwise") << spanwise_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        wakeBC100ReFvPatchVectorField
    );
}

// ************************************************************************* //
