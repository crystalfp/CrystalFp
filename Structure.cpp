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
