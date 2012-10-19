
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <iostream>
#include <fstream>
#include <cfloat>
#include <cstring>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "CrystalFp.h"
#include "StructureList.h"
#include "FingerprintingMethod.h"
#include "DistanceMethod.h"
#include "GroupingMethod.h"
#include "DistanceMatrix.h"

static const unsigned int NO_METHOD_SELECTED = UINT_MAX;

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

using namespace cfp;
using namespace cfp_internal;

/// Private data for CrystalFp.
///
struct CrystalFp::CrystalFpImpl
{
	StructureList							mSL;						///< The list of data for each structure
	unsigned int							mVerbose;					///< Verbose level: 0 - no msg; 1 - moderate; 2 - full
	float									mCutoffDistance;			///< Max distance to which fingerprint should be computed
	bool									mIsNanocluster;				///< If the sets of structures are not crystals, but something as nanoclusters
	float									mDiffrBinSize;				///< Value of the bin size used in the current fingerprint computation
	float									mDiffrPeakSize;				///< Value of the gaussian smoothing used by the current fingerprint computation
	unsigned int							mForcedFpLen;				///< Length to which the read back from checkpoint fingerprint should be forced (0 means not set)
	std::string								mCheckpointDir;				///< Directory for storing/reading the checkpointed fingerprints
	unsigned int							mFingerprintMethodIdx;		///< Fingerprint computation method used
	FingerprintingMethodsList				mFingerprintMethodsList;	///< List of fingerprinting methods
	unsigned int							mDistanceMethodIdx;			///< Distance measure type used
	DistanceMethodsList						mDistanceMethodsList;		///< List of distance measure methods
	DistanceMatrix							mDistanceMatrix;			///< Distances between fingerprints
	unsigned int							mGroupingMethodIdx;			///< The classification method used
	GroupingMethodsList						mGroupingMethodsList;		///< List of grouping methods
	unsigned int							mK;							///< The common neighbors count used by some grouping methods
	float									mMaxDistanceForGrouping; 	///< Classify tolerance used for grouping
	std::vector< std::set<unsigned int> >	mGroupedResults;			///< Grouped results
	unsigned int 							mNumGroupedEntries;			///< This number of first entries in vector_of_groups contains more than one member
	unsigned int 							mNumSingleEntries;			///< The remaining entries with only one entry
};


CrystalFp::CrystalFp(unsigned int aVerboseLevel) : mPimpl(new CrystalFpImpl)
{
	mPimpl->mVerbose = aVerboseLevel;
	resetAll();
}


CrystalFp::~CrystalFp()
{
	resetAll();
	delete mPimpl;
}


void CrystalFp::resetAll(void)
{
	mPimpl->mSL.clear();
	mPimpl->mCutoffDistance			= 0;
	mPimpl->mIsNanocluster			= false;
	mPimpl->mDiffrBinSize			= 0.05F;
	mPimpl->mDiffrPeakSize			= 0.02F;
	mPimpl->mForcedFpLen			= 0;
	mPimpl->mCheckpointDir.clear();
	mPimpl->mFingerprintMethodIdx	= NO_METHOD_SELECTED;
	mPimpl->mDistanceMethodIdx		= NO_METHOD_SELECTED;
	mPimpl->mDistanceMatrix.clear();
	mPimpl->mGroupingMethodIdx		= NO_METHOD_SELECTED;
	mPimpl->mK						= 0;
	mPimpl->mMaxDistanceForGrouping	= 0.01F;
	mPimpl->mGroupedResults.clear();
	mPimpl->mNumGroupedEntries		= 0;
	mPimpl->mNumSingleEntries		= 0;
}


void CrystalFp::addStructure(int aStep, unsigned int aNumAtoms, const float *aCoords, const unsigned int *aZ, const float *aUnitCell, bool aHasEnergy, float aEnergy, bool aEnergyIsPerAtom)
{
	if(aHasEnergy && !aEnergyIsPerAtom) aEnergy /= (float)aNumAtoms;

	Structure s(aStep, aNumAtoms, aCoords, aZ, aUnitCell, aEnergy, aHasEnergy);

	mPimpl->mSL.addStructure(s);
	mPimpl->mSL.selectAll();
}


void CrystalFp::addStructureBatch(int aStep, unsigned int aNumAtoms, const float *aCoords, const unsigned int *aZ, const float *aUnitCell, bool aHasEnergy, float aEnergy, bool aEnergyIsPerAtom)
{
	if(aHasEnergy && !aEnergyIsPerAtom) aEnergy /= (float)aNumAtoms;

	Structure s(aStep, aNumAtoms, aCoords, aZ, aUnitCell, aEnergy, aHasEnergy);

	mPimpl->mSL.addStructure(s);
}


void CrystalFp::addStructureBatchFinish(void)
{
	mPimpl->mSL.selectAll();
}


const std::vector<std::string> CrystalFp::getFingerprintMethodsNames(void) const
{
	std::vector<std::string> names;

	FingerprintingMethodsList::const_iterator ifml;
	for(ifml=mPimpl->mFingerprintMethodsList.begin(); ifml != mPimpl->mFingerprintMethodsList.end(); ++ifml)
	{
		names.push_back((*ifml)->getName());
	}

	return names;
}


const std::vector<std::string> CrystalFp::getDistanceMethodsNames(void) const
{
	std::vector<std::string> names;

	DistanceMethodsList::const_iterator ifml;
	for(ifml=mPimpl->mDistanceMethodsList.begin(); ifml != mPimpl->mDistanceMethodsList.end(); ++ifml)
	{
		names.push_back((*ifml)->getName());
	}

	return names;
}


const std::vector<std::string> CrystalFp::getGroupingMethodsNames(void) const
{
	std::vector<std::string> names;

	GroupingMethodsList::const_iterator ifml;
	for(ifml=mPimpl->mGroupingMethodsList.begin(); ifml != mPimpl->mGroupingMethodsList.end(); ++ifml)
	{
		names.push_back((*ifml)->getName());
	}

	return names;
}


size_t CrystalFp::getNumActiveStructures(void) const
{
	return mPimpl->mSL.getSelectedCount();
}


size_t CrystalFp::getNumTotalStructures(void) const
{
	return mPimpl->mSL.getTotalCount();
}


bool CrystalFp::hasUnitCell(void) const
{
	if(mPimpl->mIsNanocluster) return false;

	StructureList::const_iterator it;
	for(it=mPimpl->mSL.begin(); it != mPimpl->mSL.end(); ++it)
	{
		if((*it)->second.mUnitCell[15] == 0) return false;
	}

	return true;
}


int CrystalFp::idxToStep(size_t aIdx) const
{
	Structure& s = mPimpl->mSL.getStructureByIndex(aIdx);

	return s.mStepId;
}


bool CrystalFp::hasEnergies(void) const
{
	StructureList::const_iterator it;
	for(it=mPimpl->mSL.begin(); it != mPimpl->mSL.end(); ++it)
	{
		if(!(*it)->second.mHasEnergy) return false;
	}

	return true;
}


float CrystalFp::getMinEnergy(void) const
{
	float curr_min_energy = FLT_MAX;

	StructureList::const_iterator it;
	for(it=mPimpl->mSL.begin(); it != mPimpl->mSL.end(); ++it)
	{
		if(!(*it)->second.mHasEnergy) return -FLT_MAX;

		float energy = (*it)->second.mEnergyPerAtom;
		if(energy < curr_min_energy) curr_min_energy = energy;
	}

	return curr_min_energy;
}


void CrystalFp::energyThreshold(float aEnergyThreshold)
{
	mPimpl->mSL.filterOnEnergy(aEnergyThreshold);
}


void CrystalFp::noEnergyThreshold(void)
{
	mPimpl->mSL.selectAll();
}


float CrystalFp::computeCutoffDistance(float aMargin) const
{
	// Compute the minimum radius to use for the fingerprinting
	float max_basis_len = 0.0F;
	StructureList::const_iterator it;
	for(it=mPimpl->mSL.begin(); it != mPimpl->mSL.end(); ++it)
	{
		float len = (*it)->second.getMaxDiagonalLength();
		if(len > max_basis_len) max_basis_len = len;
	}

	// Return the cutoff distance (adding a security margin)
	return (max_basis_len / 2.0F) * (1.F + aMargin);
}


void CrystalFp::setCutoffDistance(float aCutoff)
{
	mPimpl->mCutoffDistance = aCutoff;
}


void CrystalFp::setNanoclusterStructureType(void)
{
	mPimpl->mIsNanocluster = true;
}


void CrystalFp::setDiffrBinSize(float aBinSize)
{
	mPimpl->mDiffrBinSize = aBinSize;
}


float CrystalFp::getDiffrBinSize(void) const
{
	return mPimpl->mDiffrBinSize;
}


void CrystalFp::setDiffrPeakSize(float aPeakSize)
{
	mPimpl->mDiffrPeakSize = aPeakSize;
}


float CrystalFp::getDiffrPeakSize(void) const
{
	return mPimpl->mDiffrPeakSize;
}


void CrystalFp::forceFpLength(unsigned int aDim)
{
	mPimpl->mForcedFpLen = aDim;
}


void CrystalFp::setCheckpointDir(const char* aDir)
{
	if(!aDir || aDir[0] == '\0') mPimpl->mCheckpointDir.clear();
	else                         mPimpl->mCheckpointDir.assign(aDir);
}


void CrystalFp::loadCheckpoint(void)
{
	if(mPimpl->mCheckpointDir.empty())
	{
		std::cerr << "Checkpoint dir not set. Ignoring" << std::endl;
		return;
	}

	// Temporary buffer in which the fingerprint is read
	float* temp_fingerprint = 0;
	unsigned int temp_fingerprint_len = 0;

	// The checkpoint directory contains a file for each structure
	// The filename is FPsxxxxxxxxxx.dat
	// Where s is the sign of the step number: P +; N -
	// x are 10 digits zero padded
	//
	// For each structure
	StructureList::const_iterator it;
	for(it=mPimpl->mSL.begin(); it != mPimpl->mSL.end(); ++it)
	{
		std::string filename(mPimpl->mCheckpointDir);
		if((*it)->first > 0) filename.append("/FPP"); else filename.append("/FPN");

		char nb[15];
		sprintf(nb, "%010d.dat", abs((*it)->first));
		filename.append(nb);

		// Open the checkpoint file
		std::ifstream infile(filename.c_str(), std::ios_base::in | std::ios_base::binary);
		if(infile.good())
		{
			// The file contains the following values:
			//
			//		step number					(int)
			//		fingerprint total length	(unsigned int)
			//		fingerprint parts			(unsigned int)
			//		fingerprint					(float[total length])
			//
			int step;
			infile.read((char *)&step, sizeof(int));
			unsigned int fingerprint_total_length;
			infile.read((char *)&fingerprint_total_length, sizeof(unsigned int));
			unsigned int fingerprint_parts;
			infile.read((char *)&fingerprint_parts, sizeof(unsigned int));
			if(fingerprint_total_length > temp_fingerprint_len)
			{
				delete [] temp_fingerprint;
				temp_fingerprint = new float[fingerprint_total_length];
				temp_fingerprint_len = fingerprint_total_length;
			}
			infile.read((char *)temp_fingerprint, sizeof(float)*fingerprint_total_length);

			// Fill the Structure class eventually reducing the dimensionality of the fingerprint
			if(mPimpl->mForcedFpLen > 0 && mPimpl->mForcedFpLen < fingerprint_total_length)
			{
				unsigned int plen = mPimpl->mForcedFpLen / fingerprint_parts;
				(*it)->second.mFingerprint.resize(plen*fingerprint_parts);
				for(unsigned int i=0; i < fingerprint_parts; ++i)
				{
					for(unsigned int j=0; j < plen; ++j)
					{
						(*it)->second.mFingerprint[i*plen+j] = temp_fingerprint[i*(fingerprint_total_length/fingerprint_parts)+j];
					}
				}
				(*it)->second.mFingerprintNumSections = fingerprint_parts;
				(*it)->second.mFingerprintSectionLen = plen;
			}
			else
			{
				(*it)->second.mFingerprint.assign(temp_fingerprint, temp_fingerprint+fingerprint_total_length);
				(*it)->second.mFingerprintNumSections = fingerprint_parts;
				(*it)->second.mFingerprintSectionLen = fingerprint_total_length/fingerprint_parts;
			}

			// Read it
			infile.close();
		}
	}
	
	// Release the buffer
	delete [] temp_fingerprint;
}


void CrystalFp::setFingerprintMethod(unsigned int aFingerprintType)
{
	mPimpl->mFingerprintMethodIdx = aFingerprintType;
	if(mPimpl->mFingerprintMethodIdx == NO_METHOD_SELECTED) return;

	// On invalid type throw exception
	if(mPimpl->mFingerprintMethodIdx >= mPimpl->mFingerprintMethodsList.size()) throw CrystalFpFatal("Invalid Fingerprinting method number");
}


const std::string CrystalFp::getFingerprintMethod(void) const
{
	const FingerprintingMethod* method = mPimpl->mFingerprintMethodsList.getMethod(mPimpl->mFingerprintMethodIdx);
	return method->getName();
}


void CrystalFp::unitCellInverse(const float *aUnitCell, float aUnitCellInverse[3][3])
{
	// Calc the adjoint
	// Each cofactor is the determinant of the minor of the corresponding element
	aUnitCellInverse[0][0] = aUnitCell[5]*aUnitCell[10] - aUnitCell[9]*aUnitCell[6];
	aUnitCellInverse[1][0] = aUnitCell[9]*aUnitCell[2]  - aUnitCell[1]*aUnitCell[10];
	aUnitCellInverse[2][0] = aUnitCell[1]*aUnitCell[6]  - aUnitCell[5]*aUnitCell[2];
	aUnitCellInverse[0][1] = aUnitCell[8]*aUnitCell[6]  - aUnitCell[4]*aUnitCell[10];
	aUnitCellInverse[1][1] = aUnitCell[0]*aUnitCell[10] - aUnitCell[8]*aUnitCell[2];
	aUnitCellInverse[2][1] = aUnitCell[4]*aUnitCell[2]  - aUnitCell[0]*aUnitCell[6];
	aUnitCellInverse[0][2] = aUnitCell[4]*aUnitCell[9]  - aUnitCell[8]*aUnitCell[5];
	aUnitCellInverse[1][2] = aUnitCell[8]*aUnitCell[1]  - aUnitCell[0]*aUnitCell[9];
	aUnitCellInverse[2][2] = aUnitCell[0]*aUnitCell[5]  - aUnitCell[4]*aUnitCell[1];

	// Calc the determinant, with the diagonal
	float D = aUnitCell[0]*aUnitCellInverse[0][0] + aUnitCell[4]*aUnitCellInverse[1][0] + aUnitCell[8]*aUnitCellInverse[2][0];

	// Divide adjoint by determinant to obtain the inverse matrix
	aUnitCellInverse[0][0] /= D;
	aUnitCellInverse[1][0] /= D;
	aUnitCellInverse[2][0] /= D;
	aUnitCellInverse[0][1] /= D;
	aUnitCellInverse[1][1] /= D;
	aUnitCellInverse[2][1] /= D;
	aUnitCellInverse[0][2] /= D;
	aUnitCellInverse[1][2] /= D;
	aUnitCellInverse[2][2] /= D;
}


void CrystalFp::computeExpansion(const float* aUnitCell, unsigned int* aExpansion) const
{
	// Test points on a sphere of radius 1
	const int v_div = 8;
	const int h_div = 8;
	const int num_try_points = h_div*(v_div-1)+2;
	float try_points[num_try_points][3];

	try_points[0][0] = 1.0F;
	try_points[0][1] = 0.0F;
	try_points[0][2] = 0.0F;

	int j = 1;
	for(int v=1; v < v_div; ++v)
	{
		float first = (float)cos(M_PI/v_div*v);
		for(int h=0; h < h_div; ++h)
		{
			try_points[j][0] = first;
			try_points[j][1] = (float)cos(M_PI/h_div*h);
			try_points[j][2] = (float)sin(M_PI/h_div*h);
			++j;
		}
	}
	try_points[j][0] = -1.0F;
	try_points[j][1] =  0.0F;
	try_points[j][2] =  0.0F;

	// Invert the matrix to obtain the fractional coordinates from the cartesian ones
	float minv[3][3];
	unitCellInverse(aUnitCell, minv);

	// Compute the enlargement factors (ie. the number of unit cell copies on each side of the original one)
	unsigned int ex = 1;
	unsigned int ey = 1;
	unsigned int ez = 1;
	unsigned int n;

	// For each try point
	for(int i=0; i < num_try_points; ++i)
	{
		// Enlarge to account for truncation
		float x = aUnitCell[0] + aUnitCell[4] + aUnitCell[8]  + try_points[i][0]*(mPimpl->mCutoffDistance+4*mPimpl->mDiffrPeakSize);
		float y = aUnitCell[1] + aUnitCell[5] + aUnitCell[9]  + try_points[i][1]*(mPimpl->mCutoffDistance+4*mPimpl->mDiffrPeakSize);
		float z = aUnitCell[2] + aUnitCell[6] + aUnitCell[10] + try_points[i][2]*(mPimpl->mCutoffDistance+4*mPimpl->mDiffrPeakSize);

		float xf = x*minv[0][0]+y*minv[0][1]+z*minv[0][2];
		float yf = x*minv[1][0]+y*minv[1][1]+z*minv[1][2];
		float zf = x*minv[2][0]+y*minv[2][1]+z*minv[2][2];

		if(xf < 0.F)
		{
			xf = -xf; // this is not true from a frac coords point of view, but the result is the same
			n = (unsigned int)ceil(xf);
			if(n > ex) ex = n;
		}
		else if(xf > 2.F)
		{
			n = (unsigned int)ceil(xf) - 1;
			if(n > ex) ex = n;
		}
		if(yf < 0.F)
		{
			yf = -yf;
			n = (unsigned int)ceil(yf);
			if(n > ey) ey = n;
		}
		else if(yf > 2.F)
		{
			n = (unsigned int)ceil(yf) - 1;
			if(n > ey) ey = n;
		}
		if(zf < 0.F)
		{
			zf = -zf;
			n = (unsigned int)ceil(zf);
			if(n > ez) ez = n;
		}
		else if(zf > 2.F)
		{
			n = (unsigned int)ceil(zf) - 1;
			if(n > ez) ez = n;
		}
	}
	
	// Return results
	aExpansion[0] = ex;
	aExpansion[1] = ey;
	aExpansion[2] = ez;
}


void CrystalFp::computeFingerprints(void)
{
	// Select the method to use
	if(mPimpl->mFingerprintMethodIdx == NO_METHOD_SELECTED) return;
	FingerprintingMethod* method = mPimpl->mFingerprintMethodsList.getMethod(mPimpl->mFingerprintMethodIdx);

	// If not set, then compute the cutoff distance
	if(mPimpl->mCutoffDistance <= 0) mPimpl->mCutoffDistance = computeCutoffDistance();

	// For each structure
	StructureList::const_iterator it;
	for(it=mPimpl->mSL.begin(); it != mPimpl->mSL.end(); ++it)
	{
		// If already loaded from checkpoint, don't recompute
		if(!(*it)->second.mFingerprint.empty()) continue;

		// Compute the infinite slab for crystal structures
		unsigned int expansion[3];
		if(mPimpl->mIsNanocluster)
		{
			expansion[0] = 0;
			expansion[1] = 0;
			expansion[2] = 0;
		}
		else
		{
			computeExpansion((*it)->second.mUnitCell, expansion);
		}

		// Compute fingerprint
		method->setCutoffDistance(mPimpl->mCutoffDistance);
		if(mPimpl->mIsNanocluster) method->setNanoclusterStructureType();
		method->setDiffrBinSize(mPimpl->mDiffrBinSize);
		method->setDiffrPeakSize(mPimpl->mDiffrPeakSize);
		method->computeFingerprint((*it)->second, expansion);

		// If checkpointing requested, do it
		if(!mPimpl->mCheckpointDir.empty())
		{
			std::string filename(mPimpl->mCheckpointDir);
			if((*it)->first > 0) filename.append("/FPP"); else filename.append("/FPN");

			char nb[15];
			sprintf(nb, "%010d.dat", abs((*it)->first));
			filename.append(nb);

			// Create the binary file
			std::ofstream outfile(filename.c_str(), std::ios_base::binary | std::ios_base::trunc | std::ios_base::out);
			if(!outfile.good())
			{
				std::cerr << "Cannot create fingerprints checkpoint file <" << filename << "> Ignoring step: " << (*it)->first << " for checkpointing" << std::endl;
			}
			else
			{
				int step								= (*it)->first;
				unsigned int fingerprint_parts			= (*it)->second.mFingerprintNumSections;
				unsigned int fingerprint_total_length	= (*it)->second.mFingerprintSectionLen*fingerprint_parts;
				float*       fingerprint				= &((*it)->second.mFingerprint[0]);

				outfile.write((const char *)&step, sizeof(int));
				outfile.write((const char *)&fingerprint_total_length, sizeof(unsigned int));
				outfile.write((const char *)&fingerprint_parts, sizeof(unsigned int));
				outfile.write((const char *)fingerprint, sizeof(float)*fingerprint_total_length);
				outfile.close();
			}
		}
	}
}


bool CrystalFp::hasFingerprints(void) const
{
	StructureList::const_iterator it;
	for(it=mPimpl->mSL.begin(); it != mPimpl->mSL.end(); ++it)
	{
		if((*it)->second.mFingerprint.empty()) return false;
	}

	return true;
}


unsigned int CrystalFp::getFingerprintNumSections(void) const
{
	StructureList::const_iterator it = mPimpl->mSL.begin();
	return (*it)->second.mFingerprintNumSections;
}


unsigned int CrystalFp::getFingerprintSectionLen(void) const
{
	StructureList::const_iterator it = mPimpl->mSL.begin();
	return (*it)->second.mFingerprintSectionLen;
}


float CrystalFp::getCutoffDistance(void) const
{
	return mPimpl->mCutoffDistance;
}


const float* CrystalFp::getFingerprint(size_t aStructureIdx) const
{
	return &(mPimpl->mSL.getStructureByIndex(aStructureIdx).mFingerprint[0]);
}


bool CrystalFp::isDiffractionLike(void) const
{
	if(mPimpl->mFingerprintMethodIdx == NO_METHOD_SELECTED) return false;

	const FingerprintingMethod* method = mPimpl->mFingerprintMethodsList.getMethod(mPimpl->mFingerprintMethodIdx);
	return method->isDiffractionLike();
}


void CrystalFp::setDistanceMethod(unsigned int aMeasureType)
{
	mPimpl->mDistanceMethodIdx = aMeasureType;
	if(mPimpl->mDistanceMethodIdx == NO_METHOD_SELECTED) return;

	// On invalid type throw exception
	if(mPimpl->mDistanceMethodIdx >= mPimpl->mDistanceMethodsList.size()) throw CrystalFpFatal("Invalid Distance method number");
}


const std::string CrystalFp::getDistanceMethod(void) const
{
	const DistanceMethod* method = mPimpl->mDistanceMethodsList.getMethod(mPimpl->mDistanceMethodIdx);
	return method->getName();
}


void CrystalFp::computeDistanceMatrix(void)
{
	// Select the method to use
	if(mPimpl->mDistanceMethodIdx == NO_METHOD_SELECTED) return;
	DistanceMethod* method = mPimpl->mDistanceMethodsList.getMethod(mPimpl->mDistanceMethodIdx);

	// Check if fingerprints are available
	if(!hasFingerprints()) throw CrystalFpFatal("Cannot compute distances, missing fingerprints");

	// Create temporary distances triangular matrix
	unsigned int ns = static_cast<unsigned int>(mPimpl->mSL.getSelectedCount());
	unsigned int num_elem = (ns*ns-ns)/2;
	std::vector<float> dist_matrix;
	dist_matrix.resize(num_elem);

	// Create index table
	unsigned int i, j;
	std::vector<unsigned int> row_col_table;
	for(i=0; i < (ns-1); ++i)
	{
		for(j=i+1; j < ns; ++j) 
		{
			row_col_table.push_back(i);
			row_col_table.push_back(j);
		}
	}

	// Compute distances
	#pragma omp parallel for default(none) shared(num_elem, row_col_table, method, dist_matrix)
	for(int ii=0; ii < (int)num_elem; ++ii)
	{
		unsigned int row = row_col_table[2*ii+0];
		unsigned int col = row_col_table[2*ii+1];

		Structure& s1 = mPimpl->mSL.getStructureByIndex(row);
		Structure& s2 = mPimpl->mSL.getStructureByIndex(col);

		float dist = method->computeDistance(s1, s2);

		dist_matrix[ii] = dist;
	}

	// Create output matrix
	mPimpl->mDistanceMatrix.setVector(dist_matrix, mPimpl->mSL.getSelectedCount());
}


float CrystalFp::getDistance(size_t aIdx1, size_t aIdx2) const
{
	return mPimpl->mDistanceMatrix(aIdx1, aIdx2);
}


float CrystalFp::getMaxDistance(void) const
{
	return mPimpl->mDistanceMatrix.max();
}


void CrystalFp::setGroupingMethod(unsigned int aGroupingMethod)
{
	mPimpl->mGroupingMethodIdx = aGroupingMethod;
	if(mPimpl->mGroupingMethodIdx == NO_METHOD_SELECTED) return;

	// On invalid type throw exception
	if(mPimpl->mGroupingMethodIdx >= mPimpl->mGroupingMethodsList.size()) throw CrystalFpFatal("Invalid grouping method number");
}


const std::string CrystalFp::getGroupingMethod(void) const
{
	const GroupingMethod* method = mPimpl->mGroupingMethodsList.getMethod(mPimpl->mGroupingMethodIdx);
	return method->getName();
}


void CrystalFp::groupResults(void)
{
	// Select the method to use
	if(mPimpl->mGroupingMethodIdx == NO_METHOD_SELECTED) return;
	GroupingMethod* method = mPimpl->mGroupingMethodsList.getMethod(mPimpl->mGroupingMethodIdx);

	// Distances should be available
	if(!hasDistanceMatrix()) throw CrystalFpFatal("Cannot group due to missing distances");

	// Set parameters
	method->setKvalue(mPimpl->mK);
	method->setMaxDistanceForGrouping(mPimpl->mMaxDistanceForGrouping);

	// Invoke the grouping method
	mPimpl->mGroupedResults.clear();
	method->doGrouping(mPimpl->mSL.getSelectedCount(), mPimpl->mDistanceMatrix, mPimpl->mGroupedResults);

	// Count single and multi entries
	mPimpl->mNumSingleEntries  = 0;
	mPimpl->mNumGroupedEntries = 0;

	std::vector< std::set<unsigned int> >::const_iterator ivs;
	for(ivs = mPimpl->mGroupedResults.begin(); ivs != mPimpl->mGroupedResults.end(); ++ivs)
	{
		if(ivs->size() > 1) ++mPimpl->mNumGroupedEntries;
		else                ++mPimpl->mNumSingleEntries;
	}
}


void CrystalFp::setMaxGroupingDistance(float aMaxDistance)
{
	mPimpl->mMaxDistanceForGrouping = aMaxDistance;
}



float CrystalFp::getMaxGroupingDistance(void) const
{
	return mPimpl->mMaxDistanceForGrouping;
}


void CrystalFp::setK(unsigned int aK)
{
	mPimpl->mK = aK;
}


bool CrystalFp::hasDistanceMatrix(void) const
{
	return !mPimpl->mDistanceMatrix.empty();
}


unsigned int CrystalFp::getNgroups(void) const
{
	return 	mPimpl->mNumGroupedEntries;
}


unsigned int CrystalFp::getNsingle(void) const
{
	return mPimpl->mNumSingleEntries;
}


bool CrystalFp::groupingNeedsK(void) const
{
	// Select the method to use
	if(mPimpl->mGroupingMethodIdx == NO_METHOD_SELECTED) return false;
	GroupingMethod* method = mPimpl->mGroupingMethodsList.getMethod(mPimpl->mGroupingMethodIdx);
	return method->needsK();
}


float CrystalFp::getTotalEnergy(size_t aIdx) const
{
	Structure& s = mPimpl->mSL.getStructureByIndex(aIdx);
	return s.mEnergyPerAtom*s.mNumAtoms;
}


float CrystalFp::getPerAtomEnergy(size_t aIdx) const
{
	Structure& s = mPimpl->mSL.getStructureByIndex(aIdx);
	return s.mEnergyPerAtom;
}


const float* CrystalFp::getUnitCell(size_t aIdx) const
{
	Structure& s = mPimpl->mSL.getStructureByIndex(aIdx);
	return s.mUnitCell;
}


const float* CrystalFp::getWeights(size_t aIdx) const
{
	Structure& s = mPimpl->mSL.getStructureByIndex(aIdx);
	return &s.mWeights[0];
}


unsigned int CrystalFp::getNatoms(size_t aIdx) const
{
	Structure& s = mPimpl->mSL.getStructureByIndex(aIdx);
	return s.mNumAtoms;
}


const unsigned int* CrystalFp::getAtomZ(size_t aIdx) const
{
	Structure& s = mPimpl->mSL.getStructureByIndex(aIdx);
	return &s.mAtomZ[0];
}


const float* CrystalFp::getCoords(size_t aIdx) const
{
	Structure& s = mPimpl->mSL.getStructureByIndex(aIdx);
	return &s.mCoordinates[0];
}


const std::vector< std::vector<float> >& CrystalFp::getInteratomicDistances(unsigned int aIdx) const
{
	Structure& s = mPimpl->mSL.getStructureByIndex(aIdx);
	return s.mInteratomicDistances;
}


unsigned int CrystalFp::reduceDuplicatesToRepresentative(std::vector<unsigned int>& aNewIndexList)
{
	unsigned int i;

	// Nothing to do
	if(getNgroups() < 1 || !hasFingerprints() || !hasDistanceMatrix())
	{
		for(i=0; i < getNumActiveStructures(); ++i) aNewIndexList.push_back(i);
		return 0;
	}

	// Build a list of included structures ordered as the input file
	std::vector<bool> included;
	included.resize(getNumActiveStructures(), true);

	std::vector< std::set<unsigned int> >::const_iterator ivs;
	std::set<unsigned int>::const_iterator is;
	if(hasEnergies())
	{
		// For each group
		for(ivs = mPimpl->mGroupedResults.begin(); ivs != mPimpl->mGroupedResults.end(); ++ivs)
		{
			// Groups with  only one element are included by default
			if(ivs->size() < 2) continue;

			// Compute energy min and max inside the group
			float        emin =  FLT_MAX;
			unsigned int imin = UINT_MAX;
			for(is = ivs->begin(); is != ivs->end(); ++is)
			{
				float e = getTotalEnergy(*is);
				if(e < emin) {emin = e; imin = *is;}
			}

			// Remove all except the minimum energy structure
			for(is = ivs->begin(); is != ivs->end(); ++is)
			{
				if(*is != imin) included[*is] = false;
			}
		}
	}
	else if(getNgroups() > 1)
	{
		std::vector< std::set<unsigned int> >::const_iterator ivs2;
		std::set<unsigned int>::const_iterator is2;
		for(ivs = mPimpl->mGroupedResults.begin(); ivs != mPimpl->mGroupedResults.end(); ++ivs)
		{
			if(ivs->size() < 2) continue;

			float        max_silhouette     = -2.;
			unsigned int max_silhouette_idx = UINT_MAX;

			for(is = ivs->begin(); is != ivs->end(); ++is)
			{
				// Mean distance inside the group
				float a = 0.0F;
				for(is2 = ivs->begin(); is2 != ivs->end(); ++is2)
				{
					if(is != is2) a += getDistance(*is, *is2);
				}
				a /= (ivs->size()-1);

				// For each other group
				float b = FLT_MAX;
				for(ivs2 = mPimpl->mGroupedResults.begin(); ivs2 != mPimpl->mGroupedResults.end(); ++ivs2)
				{
					if(ivs == ivs2) continue;

					float bt = 0.0F;
					for(is2 = ivs2->begin(); is2 != ivs2->end(); ++is2)
					{
						bt += mPimpl->mDistanceMatrix(*is, *is2);
					}
					bt /= ivs2->size();

					// Find the minimum distance
					if(bt < b) b = bt;
				}

				// Compute the silhouette coefficient and save it in the output matrix
				float sc = (b - a) / ((a > b) ? a : b);

				// Conserve the maximum silhouette point
				if(sc > max_silhouette)
				{
					max_silhouette = sc;
					max_silhouette_idx = *is;
				}
			}
			for(is = ivs->begin(); is != ivs->end(); ++is)
			{
				if(*is != max_silhouette_idx) included[*is] = false;
			}
		}
	}
	else // Only one group, select its first element
	{
		std::vector< std::set<unsigned int> >::const_iterator ivs;
		for(ivs = mPimpl->mGroupedResults.begin(); ivs != mPimpl->mGroupedResults.end(); ++ivs)
		{
			if(ivs->size() < 2) continue;

			std::set<unsigned int>::const_iterator is = ivs->begin();
			for(++is; is != ivs->end(); ++is)
			{
				included[*is] = false;
			}

			// There is only one group, so stop here
			break;
		}
	}

	// Mark as deselected the structures
	unsigned int ns = static_cast<unsigned int>(getNumActiveStructures());
	mPimpl->mSL.selectAllIncluded(included);
	unsigned int rns = static_cast<unsigned int>(getNumActiveStructures());

	// Recompute the distances
	mPimpl->mDistanceMatrix.resizeToIncluded(included);

	// Write the selection methood
	if(mPimpl->mVerbose > 0)
	{
		if(hasEnergies())
		{
			std::cerr << "Selection using min energy" << std::endl;
		}
		else if(getNgroups() > 1)
		{
			std::cerr << "Selection using max silhouette" << std::endl;
		}
		else
		{
			std::cerr << "Selection using first element" << std::endl;
		}
	}

	// Write the new list of indices
	for(i=0; i < ns; ++i) if(included[i]) {aNewIndexList.push_back(i);}

	return rns;
}


const std::vector< std::set<unsigned int> >& CrystalFp::getGroups(void) const
{
	return mPimpl->mGroupedResults;
}

// For the class serialize
#define FORMAT_VERSION			1
#define MAGIC_LEN				32
#define MAGIC_STRING			"CrystalFpSerialize"
#if defined(__LP64__) || defined(_LP64)
#define MAGIC_NUMBER			0x466C617473797243	// "CrystalF" as size_t
#else
#define MAGIC_NUMBER			0x73797243			// "Crys" as size_t
#endif

void CrystalFp::serialize(std::ofstream& aStream) const
{
	// Write the header line in the file (length will be 32 bytes null padded)
	char magic[MAGIC_LEN];
	memset(magic, 0, MAGIC_LEN);
	sprintf(magic, "%s V%d\n", MAGIC_STRING, FORMAT_VERSION);
	aStream.write(magic, MAGIC_LEN);

	// Write number of structures and the structures themself
	mPimpl->mSL.serialize(aStream);

	aStream.write((char *)&mPimpl->mCutoffDistance, sizeof(float));
	unsigned int x = (mPimpl->mIsNanocluster) ? 1 : 0;
	aStream.write((char *)&x, sizeof(unsigned int));
	aStream.write((char *)&mPimpl->mDiffrBinSize, sizeof(float));
	aStream.write((char *)&mPimpl->mDiffrPeakSize, sizeof(float));
	aStream.write((char *)&mPimpl->mForcedFpLen, sizeof(unsigned int));
	x = static_cast<unsigned int>(mPimpl->mCheckpointDir.size());
	aStream.write((char *)&x, sizeof(unsigned int));
	if(x) aStream.write(mPimpl->mCheckpointDir.c_str(), x);
	aStream.write((char *)&mPimpl->mFingerprintMethodIdx, sizeof(unsigned int));
	aStream.write((char *)&mPimpl->mDistanceMethodIdx, sizeof(unsigned int));
	aStream.write((char *)&mPimpl->mGroupingMethodIdx, sizeof(unsigned int));
	aStream.write((char *)&mPimpl->mK, sizeof(unsigned int));
	aStream.write((char *)&mPimpl->mMaxDistanceForGrouping, sizeof(float));
	aStream.write((char *)&mPimpl->mNumGroupedEntries, sizeof(unsigned int));
	aStream.write((char *)&mPimpl->mNumSingleEntries, sizeof(unsigned int));
	mPimpl->mDistanceMatrix.serialize(aStream);
	x = static_cast<unsigned int>(mPimpl->mGroupedResults.size());
	aStream.write((char *)&x, sizeof(unsigned int));
	for(unsigned int i=0; i < x; ++i)
	{
		unsigned int y = static_cast<unsigned int>(mPimpl->mGroupedResults.size());
		aStream.write((char *)&y, sizeof(unsigned int));
		std::set<unsigned int>::const_iterator is;
		for(is=mPimpl->mGroupedResults[i].begin(); is != mPimpl->mGroupedResults[i].end(); ++is)
		{
			unsigned int z = *is;
			aStream.write((char *)&z, sizeof(unsigned int));
		}
	}
}


void CrystalFp::unserialize(std::ifstream& aStream, bool aAppend, int aStepOffset)
{
	// Verify the magic number
	size_t magic_number;
	aStream.read((char *)&magic_number, sizeof(size_t));
	if(magic_number != MAGIC_NUMBER)
	{
		aStream.close();
		throw CrystalFpFatal("Invalid magic number for save file to be reloaded");
	}

	// Verify version
	char magic[MAGIC_LEN];
	aStream.read(magic, MAGIC_LEN-sizeof(size_t));
	int file_version = atoi(magic+sizeof(MAGIC_STRING)-sizeof(size_t)+1);
	if(file_version != FORMAT_VERSION)
	{
		aStream.close();
		std::cerr << "Invalid version for save file to be reloaded" << std::endl;
		std::cerr << "Version is: " << file_version << " instead of: " << FORMAT_VERSION << std::endl;
		throw CrystalFpFatal();
	}

	// Read the set of Structures
	mPimpl->mSL.unserialize(aStream, aAppend, aStepOffset);

	unsigned int x;
	aStream.read((char *)&mPimpl->mCutoffDistance, sizeof(float));
	aStream.read((char *)&x, sizeof(unsigned int));
	mPimpl->mIsNanocluster = (x != 0);
	aStream.read((char *)&mPimpl->mDiffrBinSize, sizeof(float));
	aStream.read((char *)&mPimpl->mDiffrPeakSize, sizeof(float));
	aStream.read((char *)&mPimpl->mForcedFpLen, sizeof(unsigned int));
	aStream.read((char *)&x, sizeof(unsigned int));
	if(x)
	{
		char *str = new char[x];
		aStream.read(str, x);
		mPimpl->mCheckpointDir.assign(str, x);
		delete [] str;
	}
	else
	{
		mPimpl->mCheckpointDir.clear();
	}
	aStream.read((char *)&mPimpl->mFingerprintMethodIdx, sizeof(unsigned int));
	aStream.read((char *)&mPimpl->mDistanceMethodIdx, sizeof(unsigned int));
	aStream.read((char *)&mPimpl->mGroupingMethodIdx, sizeof(unsigned int));
	aStream.read((char *)&mPimpl->mK, sizeof(unsigned int));
	aStream.read((char *)&mPimpl->mMaxDistanceForGrouping, sizeof(float));
	aStream.read((char *)&mPimpl->mNumGroupedEntries, sizeof(unsigned int));
	aStream.read((char *)&mPimpl->mNumSingleEntries, sizeof(unsigned int));
	mPimpl->mDistanceMatrix.unserialize(aStream);
	aStream.read((char *)&x, sizeof(unsigned int));
	mPimpl->mGroupedResults.clear();
	for(unsigned int i=0; i < x; ++i)
	{
		std::set<unsigned int> s;
		s.clear();
		unsigned int y;
		aStream.read((char *)&y, sizeof(unsigned int));
		for(unsigned int j=0; j < y; ++j)
		{
			unsigned int z;
			aStream.read((char *)&z, sizeof(unsigned int));
			s.insert(z);
		}
		mPimpl->mGroupedResults.push_back(s);
	}
}


void CrystalFp::dump(void) const
{
	std::cerr << "Num. structures: " << mPimpl->mSL.getSelectedCount() << std::endl;
}


