#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "CrystalFpScatterplot.h"

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif


//using namespace cfp_internal;
using namespace cfp;


/// The private data for the CrystalFpAnalysis class.
///
struct CrystalFpScatterplot::CrystalFpScatterplotImpl
{
	float				mMass;
	float				mStiffness;
	float				mDamping;
	float				mPerturbScale;
	float				mPerturbDecay;

	unsigned int		mNumPoints;

	std::vector<float>	mPositions;			///< Array of particle positions
	std::vector<float>	mForces;			///< Forces on the particles
	std::vector<float>	mVelocities;		///< Particle velocities
	std::vector<float>	mStress;			///< Particle stress

	std::vector<float>	mSavedPositions;	///< Saved particle positions (to implement multistart)
	float				mMaxStress;			///< Stress value for the saved position

	DiagnosticType		mDiagnostic;		///< The requested scatterplot diagnostic type
	unsigned int		mNumBins;			///< Number of bins for binned diagnostic
	float				mWobbleScale;		///< A small random perturbation to the binned result points to avoid visual artefacts.
};


CrystalFpScatterplot::CrystalFpScatterplot() : mPimpl(new CrystalFpScatterplotImpl())
{
	mCrystalFp				= 0;

	mPimpl->mMass			= 10.0F;
	mPimpl->mStiffness		=  1.0F;
	mPimpl->mDamping		=  0.7F;
	mPimpl->mPerturbScale	=  0.2F;
	mPimpl->mPerturbDecay	=  0.1F;

	mPimpl->mNumPoints		= 0;
	mPimpl->mMaxStress		= FLT_MAX;
	mPimpl->mDiagnostic		= DIAG_DO_NOTHING;
	mPimpl->mNumBins		= 100;
	mPimpl->mWobbleScale	= 0.001F;
}

CrystalFpScatterplot::~CrystalFpScatterplot()
{
	delete mPimpl;
}


void CrystalFpScatterplot::setNamedParam(const std::string& aName, const std::string& aValue)
{
	if(!strncasecmp(aName.c_str(), "mass", 1))
	{
		mPimpl->mMass = (float)atof(aValue.c_str());
	}
	else if(!strncasecmp(aName.c_str(), "stiffness", 1))
	{
		mPimpl->mStiffness = (float)atof(aValue.c_str());
	}
	else if(!strncasecmp(aName.c_str(), "damping", 1))
	{
		mPimpl->mDamping = (float)atof(aValue.c_str());
	}
	else if(!strncasecmp(aName.c_str(), "perturb", 1))
	{
		mPimpl->mPerturbScale = (float)atof(aValue.c_str());
	}
	else if(!strncasecmp(aName.c_str(), "retry", 1))
	{
		int retry = atoi(aValue.c_str());
		if(retry <= 0) retry = 1;

		mPimpl->mPerturbDecay = (retry > 1) ? (float)pow(0.001, 1/(retry-1)) : 1.0F;
	}
	else if(!strncasecmp(aName.c_str(), "bins", 1))
	{
		int bins = atoi(aValue.c_str());
		if(bins < 1) bins = 1;
		mPimpl->mNumBins = bins;
	}
	else if(!strncasecmp(aName.c_str(), "wobble", 1))
	{
		mPimpl->mWobbleScale = (float)atof(aValue.c_str());
	}
	else
	{
		throw CrystalFpFatal("Invalid named param for CrystalFpScatterplot");
	}
}

unsigned int CrystalFpScatterplot::initScatterplot(const CrystalFp *aCfp)
{
	mCrystalFp = aCfp;

	// To compute the scatterplot, the distances should be available
	if(!mCrystalFp->hasDistanceMatrix()) throw CrystalFpFatal("No distances available");
	
	// Initialize the random number generator
	srand((unsigned)time(NULL));

	// Initialize the arrays
	unsigned int ns = mCrystalFp->getNumActiveStructures();
	mPimpl->mNumPoints = ns;
	mPimpl->mForces.resize(2*ns);

	// Initialize random positions
	unsigned int i;
	mPimpl->mPositions.resize(2*ns);
	for(i=0; i < 2*ns; ++i) mPimpl->mPositions[i] = 2.0F*((float)rand()/(float)RAND_MAX - 0.5F);

	// Initialize velocities
	mPimpl->mVelocities.assign(2*ns, 0.0F);

	// Initialize stress (in case colored by stress is asked before having done a stepScatterplot())
	mPimpl->mStress.assign(ns, 0.0F);

	// Initialize the array for restarts
	mPimpl->mSavedPositions.clear();
	mPimpl->mMaxStress = FLT_MAX;

	// Return the number of scatterpoints
	return ns;
}


void CrystalFpScatterplot::getPoints(float* aCoords) const
{
    if(mPimpl->mNumPoints) memcpy(aCoords, &(mPimpl->mPositions[0]), mPimpl->mNumPoints*2*sizeof(float));
}


void CrystalFpScatterplot::colorByGroup(unsigned int aNumPoints, float *aResult) const
{
	int group;
	const std::vector< std::set<unsigned int> > all_groups = mCrystalFp->getGroups();
	std::vector< std::set<unsigned int> >::const_iterator ivs;
	std::set<unsigned int>::const_iterator is;

	for(ivs = all_groups.begin(),group=0; ivs != all_groups.end(); ++ivs, ++group)
	{
		if(ivs->size() == 1)
		{
			is = ivs->begin();
			if(*is < aNumPoints) aResult[*is] = 0;
		}
		else
		{
			for(is = ivs->begin(); is != ivs->end(); ++is)
			{
				if(*is < aNumPoints) aResult[*is] = group+10.0F;
			}
		}
	}
}


void CrystalFpScatterplot::getValues(float* aValues, ValueType aValueType) const
{
	// If no points, do nothing
	if(!mPimpl->mNumPoints) return;

	// Set the point values
	unsigned int i;
	switch(aValueType)
	{
	case VAL_TOTAL_ENERGY:
		for(i=0; i < mPimpl->mNumPoints; ++i) aValues[i] = mCrystalFp->getTotalEnergy(i);
		break;

	case VAL_PER_ATOM_ENERGY:
		for(i=0; i < mPimpl->mNumPoints; ++i) aValues[i] = mCrystalFp->getPerAtomEnergy(i);
		break;

	case VAL_STRESS:
		memcpy(aValues, &(mPimpl->mStress[0]), mPimpl->mNumPoints*sizeof(float));
		break;

	case VAL_GROUP:
		colorByGroup(mPimpl->mNumPoints, aValues);
		break;

	case VAL_STEP:
		for(i=0; i < mPimpl->mNumPoints; ++i) aValues[i] = (float)i;
		break;

	default:
		throw CrystalFpFatal("Invalid points value type");
	}
}


float CrystalFpScatterplot::stepScatterplot(float aTimestep)
{
	unsigned int ns = mPimpl->mNumPoints;
	if(!ns) return 0.0F;

	// Initialize forces and stress
	mPimpl->mForces.assign(ns*2, 0.0F);
	mPimpl->mStress.assign(ns, 0.0F);

	// Compute forces
	unsigned int i, j;
	for(i=0; i < ns-1; ++i)
	{
		for(j=i+1; j < ns; ++j)
		{
			// Vector i -> j
			float dx = mPimpl->mPositions[2*j+0] - mPimpl->mPositions[2*i+0];
			float dy = mPimpl->mPositions[2*j+1] - mPimpl->mPositions[2*i+1];

			// Normalize
			float len = sqrt(dx*dx+dy*dy);
			if(len < 1e-10F) len = 1e-10F;
			dx /= len;
			dy /= len;

			// Compute force
			float target_len = mCrystalFp->getDistance(i, j);

			//float force = mPimpl->mStiffness*(len - target_len);
			float force = (len > target_len) ? 10*mPimpl->mStiffness*(len - target_len) : mPimpl->mStiffness*(len - target_len);

			// Update the forces
			mPimpl->mForces[2*i+0] += force*dx;
			mPimpl->mForces[2*i+1] += force*dy;
			mPimpl->mForces[2*j+0] -= force*dx;
			mPimpl->mForces[2*j+1] -= force*dy;

			// Update the stress values
			if(force > 0)
			{
				mPimpl->mStress[i] += force;
				mPimpl->mStress[j] += force;
			}
			else
			{
				mPimpl->mStress[i] -= force;
				mPimpl->mStress[j] -= force;
			}
		}
	}

	// Compute new velocities and positions
	for(i=0; i < ns; ++i)
	{
		mPimpl->mVelocities[2*i+0] = mPimpl->mDamping*mPimpl->mVelocities[2*i+0] + mPimpl->mForces[2*i+0]*aTimestep/mPimpl->mMass;
		mPimpl->mVelocities[2*i+1] = mPimpl->mDamping*mPimpl->mVelocities[2*i+1] + mPimpl->mForces[2*i+1]*aTimestep/mPimpl->mMass;

		mPimpl->mPositions[2*i+0] += mPimpl->mVelocities[2*i+0]*aTimestep;
		mPimpl->mPositions[2*i+1] += mPimpl->mVelocities[2*i+1]*aTimestep;
	}

	// Return the cinetic energy of the set of points. Used to decide when to stop iterating.
	float cinetic_energy = 0;
	for(i=0; i < ns; ++i)
	{
		cinetic_energy += mPimpl->mVelocities[2*i+0]*mPimpl->mVelocities[2*i+0]+mPimpl->mVelocities[2*i+1]*mPimpl->mVelocities[2*i+1];
	}

	return cinetic_energy/2.0F;
}


void CrystalFpScatterplot::perturbPositions(void)
{
	// If first time
	if(mPimpl->mSavedPositions.empty())
	{
		mPimpl->mSavedPositions = mPimpl->mPositions;
		mPimpl->mMaxStress = computeCost();
	}
	else
	{
		float curr_cost = computeCost();

		if(curr_cost < mPimpl->mMaxStress)
		{
			// Better solution found
			mPimpl->mSavedPositions = mPimpl->mPositions;
			mPimpl->mMaxStress = curr_cost;
		}
		else
		{
			// Reload the current best solution
			mPimpl->mPositions = mPimpl->mSavedPositions;
		}
	}

	// Randomly move the points
	unsigned int ns = mPimpl->mNumPoints;
	for(unsigned int i=0; i < ns; ++i)
	{
		float dx = 2.0F*((float)rand()/(float)RAND_MAX - 0.5F);
		float dy = 2.0F*((float)rand()/(float)RAND_MAX - 0.5F);
		float len = sqrt(dx*dx+dy*dy);

		mPimpl->mPositions[2*i+0] += mPimpl->mPerturbScale*mPimpl->mStress[i]*dx/len;
		mPimpl->mPositions[2*i+1] += mPimpl->mPerturbScale*mPimpl->mStress[i]*dy/len;
	}

	// Slowly reduce the size of the perturbation
	mPimpl->mPerturbScale *= mPimpl->mPerturbDecay;
}

	
float CrystalFpScatterplot::computeCost(void) const
{
	// Median value
	std::vector<float> tmp = mPimpl->mStress;

	std::sort(tmp.begin(), tmp.end());
	unsigned int npoints = tmp.size();
	float median = (npoints % 2) ? tmp[npoints/2] : (tmp[npoints/2]+tmp[npoints/2-1])/2.0F;

	return median;
}


unsigned int CrystalFpScatterplot::initDiagnostic(DiagnosticType aDiagnostic)
{
	mPimpl->mDiagnostic = aDiagnostic;

	switch(aDiagnostic)
	{
	case DIAG_DISTANCES:
		return (mPimpl->mNumPoints * (mPimpl->mNumPoints-1))/2;

	case DIAG_BINNED_DISTANCES:
		return mPimpl->mNumBins*mPimpl->mNumBins;

	case DIAG_DO_NOTHING:
		return 0;

	default:
		throw CrystalFpFatal("Invalid scatterplot diagnostic type");
	}
}

struct Point
{
	Point(float aX, float aY, unsigned int aCnt) {mX=aX; mY=aY; mCnt=aCnt;}
	Point() {mX=0; mY=0; mCnt=0;}
	Point(const Point& aPoint) {mX=aPoint.mX; mY=aPoint.mY; mCnt=aPoint.mCnt;}

	float mX;
	float mY;
	unsigned int mCnt;
};
struct AscendingCntSort
{
     bool operator()(const Point& aPt1, const Point& aPt2)
     {
		 return aPt1.mCnt < aPt2.mCnt;
     }
};


void CrystalFpScatterplot::getDiagnosticValues(float* aCoords, float* aValues) const
{
	// Ignore 
	if(mPimpl->mDiagnostic == DIAG_DO_NOTHING) return;

	unsigned int i, j, k;
	float max_real_dist = mCrystalFp->getMaxDistance();
	unsigned int ns = mPimpl->mNumPoints;
	unsigned int nb = mPimpl->mNumBins;
	std::vector<Point> pt;
	std::vector<Point>::const_iterator ipt;

	// Find mmaximum projected distance
	float max_proj_dist = 0.0F;
	for(i=0; i < ns-1; ++i)
	{
		for(j=i+1; j < ns; ++j)
		{
			float dx = mPimpl->mPositions[2*j+0] - mPimpl->mPositions[2*i+0];
			float dy = mPimpl->mPositions[2*j+1] - mPimpl->mPositions[2*i+1];
			float len = sqrt(dx*dx+dy*dy);
			if(len > max_proj_dist) max_proj_dist = len;
		}
	}

	switch(mPimpl->mDiagnostic)
	{
	case DIAG_DISTANCES:
		for(i=k=0; i < ns-1; ++i)
		{
			for(j=i+1; j < ns; ++j)
			{
				// Projected distance
				float dx = mPimpl->mPositions[2*j+0] - mPimpl->mPositions[2*i+0];
				float dy = mPimpl->mPositions[2*j+1] - mPimpl->mPositions[2*i+1];
				float proj_dist_norm = sqrt(dx*dx+dy*dy)/max_proj_dist;

				// Real distance
				float real_dist_norm = mCrystalFp->getDistance(i, j)/max_real_dist;

				// Difference
				float diff = real_dist_norm - proj_dist_norm;
				if(diff < 0) diff = -diff;

				// Output
				aCoords[k++] = real_dist_norm;	// x: real distance
				aCoords[k++] = proj_dist_norm;	// y: projected distance
				aValues[k-2] = diff/1.414214F;	// distance from the diagonal
			}
		}
		break;

	case DIAG_BINNED_DISTANCES:
		// Output coordinates
		for(i=k=0; i < nb; ++i)
		{
			for(j=0; j <  nb; ++j)
			{
				//aCoords[k++] = (0.5F+j)/(float)nb;
				//aCoords[k++] = (0.5F+i)/(float)nb;
				if(mPimpl->mWobbleScale < 1e-5)
					pt.push_back(Point((0.5F+j)/(float)nb, (0.5F+i)/(float)nb, 0));
				else
					pt.push_back(Point((0.5F+j)/(float)nb + mPimpl->mWobbleScale*2.0F*((float)rand()/(float)RAND_MAX - 0.5F),
					                   (0.5F+i)/(float)nb + mPimpl->mWobbleScale*2.0F*((float)rand()/(float)RAND_MAX - 0.5F), 0));
			}
		}

		// Initialize counts
		//memset(aValues, 0, nb*nb*sizeof(float));

		// Bin the values
		for(i=0; i < ns-1; ++i)
		{
			for(j=i+1; j < ns; ++j)
			{
				// Projected distance
				float dx = mPimpl->mPositions[2*j+0] - mPimpl->mPositions[2*i+0];
				float dy = mPimpl->mPositions[2*j+1] - mPimpl->mPositions[2*i+1];
				float proj_dist_norm = sqrt(dx*dx+dy*dy)/max_proj_dist;

				// Real distance
				float real_dist_norm = mCrystalFp->getDistance(i, j)/max_real_dist;

				// Bin along x
				unsigned int bx = (unsigned int)(real_dist_norm * nb + 0.5F);
				if(bx >= nb) bx = nb - 1;

				// Bin along y
				unsigned int by = (unsigned int)(proj_dist_norm * nb + 0.5F);
				if(by >= nb) by = nb - 1;

				// Output
				//aValues[by*nb+bx] += 1.0F;
				pt[by*nb+bx].mCnt += 1;
			}
		}

		// Sort the array of points
		std::sort(pt.begin(), pt.end(), AscendingCntSort());

		// Output the values
		for(ipt=pt.begin(), k=0; ipt != pt.end(); ++ipt, ++k)
		{
			aCoords[2*k+0] = ipt->mX;
			aCoords[2*k+1] = ipt->mY;
			aValues[k] = (float)ipt->mCnt;
		}
		break;

	default:
		throw CrystalFpFatal("Invalid scatterplot diagnostic type");
	}
}


void CrystalFpScatterplot::dumpParams(void) const
{	
	std::cerr << std::endl;
	std::cerr << "Mass:           " << std::setw(12) << mPimpl->mMass << std::endl;
	std::cerr << "Stiffness:      " << std::setw(12) << mPimpl->mStiffness << std::endl;
	std::cerr << "Damping:        " << std::setw(12) << mPimpl->mDamping << std::endl;
	std::cerr << "Perturb scale:  " << std::setw(12) << mPimpl->mPerturbScale << std::endl;
	std::cerr << "Perturb decay:  " << std::setw(12) << mPimpl->mPerturbDecay << std::endl;

	std::cerr << "Num. points:    " << std::setw(12) << mPimpl->mNumPoints << std::endl;
	std::cerr << "Max stress:     " << std::setw(12) << mPimpl->mMaxStress << std::endl;
	std::cerr << "Diagnostic:     " << std::setw(12) << mPimpl->mDiagnostic << std::endl;
	std::cerr << "Num. bins:      " << std::setw(12) << mPimpl->mNumBins << std::endl;
	std::cerr << "Wobble scale:   " << std::setw(12) << mPimpl->mWobbleScale << std::endl;
}


