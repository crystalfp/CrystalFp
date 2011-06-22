
#ifndef CRYSTALFP_STRUCTURE_H
#define CRYSTALFP_STRUCTURE_H

#include <vector>

/// @namespace cfp_internal
/// Private namespace to hide internal implementation of the library
///
namespace cfp_internal
{


class Structure
{
public:
	// Constructors and destructor
	Structure();
	Structure(int aStepId, int aNumAtoms, const float* aCoords, const unsigned int* aZ, const float* aUnitCell, float aEnergyPerAtom, bool aHasEnergy);
	~Structure();
	
	// Longer diagonal length
	float getMaxDiagonalLength(void) const;

	// Copy constructor and assignment
	Structure(const Structure& p);
	Structure& operator=(const Structure& p);

public:
	// Structure description data
	int							mStepId;					///< Step number (it is an identifier)
	unsigned int				mNumAtoms;					///< Number of atoms in the structure
	std::vector<float>			mCoordinates;				///< Coordinates of atoms in the x, y, z... order
	std::vector<unsigned int>	mAtomZ;						///< Atom types
	float						mUnitCell[16];				///< Unit cell (mUnitCell[15] == 0 means no unit cell present)
	float						mEnergyPerAtom;				///< Optional energy for atom in this structure. If not present it is 0.
	bool						mHasEnergy;					///< True if this structure has an associated energy value

	// Fingerprint
	std::vector<float>			mFingerprint;				///< The fingerprint values
	unsigned int				mFingerprintSectionLen;		///< The length of each section of the fingerprint
	unsigned int				mFingerprintNumSections;	///< The number of sections composing the fingerprint
	std::vector<float>			mWeights;					///< Weights for the distance computation in the per type histogram

	// Interatomic distances for quasi-entropy computation
	std::vector< std::vector<float> > mInteratomicDistances;
};

}
#endif

