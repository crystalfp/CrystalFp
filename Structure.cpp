#include <cstring>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <climits>
#include "Structure.h"

using namespace cfp_internal;


Structure::Structure()
{
	// Per structure description
	mStepId = -1;	
	mNumAtoms = 0;		
	mUnitCell[15] = 0.0F; // means no unit cell
	mEnergyPerAtom = 0.0F;	
	mHasEnergy = false;

	// Fingerprint
	mFingerprintSectionLen = 0;
	mFingerprintNumSections = 0;
}


Structure::Structure(int aStepId, int aNumAtoms, const float* aCoords, const unsigned int* aZ, const float* aUnitCell, float aEnergyPerAtom, bool aHasEnergy)
{
	// Load the base data
	mStepId        = aStepId;	
	mNumAtoms      = aNumAtoms;		
	mEnergyPerAtom = aEnergyPerAtom;	
	mHasEnergy     = aHasEnergy;

	// If no atoms or no unit cell do not load it
	if(aUnitCell == 0)
	{
		mUnitCell[15] = 0.0F; // means no unit cell
	}
	else
	{
		memcpy(mUnitCell, aUnitCell, 16*sizeof(float));
	}
	mCoordinates.assign(aCoords, aCoords+3*aNumAtoms);
	mAtomZ.assign(aZ, aZ+aNumAtoms);

	// Fingerprint
	mFingerprintSectionLen = 0;
	mFingerprintNumSections = 0;
}


Structure::~Structure()
{
	mCoordinates.clear();	
	mAtomZ.clear();			
	mFingerprint.clear();
	mWeights.clear();
	mInteratomicDistances.clear();
}


// Copy constructor and assignment
Structure::Structure(const Structure& p)
{
	// Copy per structure data
	mStepId = p.mStepId;
	mNumAtoms = p.mNumAtoms;
	mEnergyPerAtom = p.mEnergyPerAtom;
	mHasEnergy = p.mHasEnergy;
	mCoordinates = p.mCoordinates;
	mAtomZ = p.mAtomZ;
	memcpy(mUnitCell, p.mUnitCell, 16*sizeof(float));

	// Copy the fingerprint data
	mFingerprint = p.mFingerprint;
	mFingerprintSectionLen = p.mFingerprintSectionLen;
	mFingerprintNumSections = p.mFingerprintNumSections;
	mWeights = p.mWeights;

	// Copy the interatomic distances
	mInteratomicDistances = p.mInteratomicDistances;
}


Structure& Structure::operator=(const Structure& p)
{
	// Make sure not same object
	if(this != &p)
	{
		// Copy per structure data
		mStepId = p.mStepId;
		mNumAtoms = p.mNumAtoms;
		mEnergyPerAtom = p.mEnergyPerAtom;
		mHasEnergy = p.mHasEnergy;
		mCoordinates = p.mCoordinates;
		mAtomZ = p.mAtomZ;
		memcpy(mUnitCell, p.mUnitCell, 16*sizeof(float));

		// Copy the fingerprint data
		mFingerprint = p.mFingerprint;
		mFingerprintSectionLen = p.mFingerprintSectionLen;
		mFingerprintNumSections = p.mFingerprintNumSections;
		mWeights = p.mWeights;

		// Copy the interatomic distances
		mInteratomicDistances = p.mInteratomicDistances;
	}

	// Return ref for multiple assignment
	return *this;
}


float Structure::getMaxDiagonalLength(void) const
{
	// First diagonal: a+b+c
	float dx = mUnitCell[0] + mUnitCell[4] + mUnitCell[8];
	float dy = mUnitCell[1] + mUnitCell[5] + mUnitCell[9];
	float dz = mUnitCell[2] + mUnitCell[6] + mUnitCell[10];
	float len = dx*dx+dy*dy+dz*dz;
	float max_basis_len = len;

	// Second diagonal: a-b+c
	dx = mUnitCell[0] - mUnitCell[4] + mUnitCell[8];
	dy = mUnitCell[1] - mUnitCell[5] + mUnitCell[9];
	dz = mUnitCell[2] - mUnitCell[6] + mUnitCell[10];
	len = dx*dx+dy*dy+dz*dz;
	if(len > max_basis_len) max_basis_len = len;

	// Third diagonal: a-b-c
	dx = mUnitCell[0] - mUnitCell[4] - mUnitCell[8];
	dy = mUnitCell[1] - mUnitCell[5] - mUnitCell[9];
	dz = mUnitCell[2] - mUnitCell[6] - mUnitCell[10];
	len = dx*dx+dy*dy+dz*dz;
	if(len > max_basis_len) max_basis_len = len;

	// Forth diagonal: a+b-c
	dx = mUnitCell[0] + mUnitCell[4] - mUnitCell[8];
	dy = mUnitCell[1] + mUnitCell[5] - mUnitCell[9];
	dz = mUnitCell[2] + mUnitCell[6] - mUnitCell[10];
	len = dx*dx+dy*dy+dz*dz;
	if(len > max_basis_len) max_basis_len = len;

	// Return the maximum value
	return sqrt(max_basis_len);
}

#include <iostream>

void Structure::serialize(std::ofstream& aStream) const
{
	aStream.write((char *)&mStepId, sizeof(int));
	aStream.write((char *)&mNumAtoms, sizeof(unsigned int));
	if(mNumAtoms) aStream.write((char *)&mCoordinates[0], sizeof(float)*3*mNumAtoms);
	if(mNumAtoms) aStream.write((char *)&mAtomZ[0], sizeof(unsigned int)*mNumAtoms);
	aStream.write((char *)mUnitCell, sizeof(float)*16);
	aStream.write((char *)&mEnergyPerAtom, sizeof(float));
	unsigned int x = (mHasEnergy) ? 1 : 0;
	aStream.write((char *)&x, sizeof(unsigned int));
	aStream.write((char *)&mFingerprintSectionLen, sizeof(unsigned int));
	aStream.write((char *)&mFingerprintNumSections, sizeof(unsigned int));
	if(mFingerprintSectionLen > 0 && mFingerprintNumSections > 0)
		aStream.write((char *)&mFingerprint[0], sizeof(float)*mFingerprintSectionLen*mFingerprintNumSections);
	x = mWeights.size();
	aStream.write((char *)&x, sizeof(unsigned int));
	if(x) aStream.write((char *)&mWeights[0], sizeof(float)*x);
	x = mInteratomicDistances.size();
	aStream.write((char *)&x, sizeof(unsigned int));
	unsigned int i;
	for(i=0; i < x; ++i)
	{
		unsigned int y = mInteratomicDistances[i].size();
		aStream.write((char *)&y, sizeof(unsigned int));
		if(y) aStream.write((char *)&mInteratomicDistances[i][0], sizeof(float)*y);
	}
}

void Structure::unserialize(std::ifstream& aStream)
{
	aStream.read((char *)&mStepId, sizeof(int));
	aStream.read((char *)&mNumAtoms, sizeof(unsigned int));
	mCoordinates.resize(3*mNumAtoms);
	if(mNumAtoms) aStream.read((char *)&mCoordinates[0], sizeof(float)*3*mNumAtoms);
	mAtomZ.resize(mNumAtoms);
	if(mNumAtoms) aStream.read((char *)&mAtomZ[0], sizeof(unsigned int)*mNumAtoms);
	aStream.read((char *)mUnitCell, sizeof(float)*16);
	aStream.read((char *)&mEnergyPerAtom, sizeof(float));
	unsigned int x;
	aStream.read((char *)&x, sizeof(unsigned int));
	mHasEnergy = (x != 0);
	aStream.read((char *)&mFingerprintSectionLen, sizeof(unsigned int));
	aStream.read((char *)&mFingerprintNumSections, sizeof(unsigned int));
	mFingerprint.resize(mFingerprintSectionLen*mFingerprintNumSections);
	if(mFingerprintSectionLen > 0 && mFingerprintNumSections > 0)
		aStream.read((char *)&mFingerprint[0], sizeof(float)*mFingerprintSectionLen*mFingerprintNumSections);
	aStream.read((char *)&x, sizeof(unsigned int));
	mWeights.resize(x);
	if(x) aStream.read((char *)&mWeights[0], sizeof(float)*x);
	aStream.read((char *)&x, sizeof(unsigned int));
	mInteratomicDistances.reserve(x);
	for(unsigned int i=0; i < x; ++i)
	{
		std::vector<float> v;
		mInteratomicDistances.push_back(v);

		unsigned int y;
		aStream.read((char *)&y, sizeof(unsigned int));

		mInteratomicDistances[i].resize(y);
		if(y) aStream.read((char *)&mInteratomicDistances[i][0], sizeof(float)*y);
	}
}


