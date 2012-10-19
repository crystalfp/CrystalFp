
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cfloat>
#include <climits>
#include "DistanceMethod.h"
#include "DistanceMatrix.h"

#ifdef _MSC_VER
#include <cmath>

namespace cfp_internal
{

double cbrt(double x)
{
	// Check for simple cases:
	if(x == 0.)
		return 0.;
	else if(x == 1.)
		return 1.;
	else if(x == -1.)
		return -1.;
	else
	{
		double y;               // Guess
		double g = fabs(x);		// Do this guess on a positive number
		double l = g * 1E-14;	// The limit for optimal guess
		// the multiplication with x (its magnitude) should
		// ensure no infinite loops, at the cost
		// of some precision on high numbers.

		// Make initial guess:
		if(g < 1)
			y = x;
		else if(g < 10)
			y = x / 3;
		else if(g < 20)
			y = x / 6;
		else if(g < 50)
			y = x / 10;
		else if(g < 100)
			y = x / 20;
		else if(g < 1000)
			y = x / 50;
		else if(g < 5000)
			y = x / 100;
		else if(g < 10000)
			y = x / 500;
		else if(g < 50000)
			y = x / 1000;
		else if(g < 100000)
			y = x / 50000;
		else
			y = x / 100000;

		// Improve guess immediately:
		y = ((x / (y * y)) + 2. * y) / 3.;      // Newton's approx. for new guess
		double d = fabs(y * y * y - x);         // Calculate difference (Last difference of y^3 and x)
		while(l < d)
		{
			y = ((x / (y * y)) + 2. * y) / 3.;  // Newton's approx. for new guess
			d = fabs(y * y * y - x);            // Calculate difference
		}

		return y;
	}
}

}
#endif

using namespace cfp_internal;

float CosineDistance::computeDistance(const Structure& aStructure1, const Structure& aStructure2) const
{
	unsigned int i;
	unsigned int nsect1 = aStructure1.mFingerprintNumSections;
	unsigned int sectlen1 = aStructure1.mFingerprintSectionLen;

	// Only one section
	if(nsect1 == 1)
	{
		double distance = 0.;
		double a_norm   = 0.;
		double b_norm   = 0.;
		for(i=0; i < sectlen1; ++i)
		{
			distance += aStructure1.mFingerprint[i] * aStructure2.mFingerprint[i];
			a_norm   += aStructure1.mFingerprint[i] * aStructure1.mFingerprint[i];
			b_norm   += aStructure2.mFingerprint[i] * aStructure2.mFingerprint[i];
		}
		distance /= sqrt(a_norm*b_norm);
		distance = (1. - distance)/2.;

		return (float)distance;
	}

	// More than one section, created by pseudo diffraction per atom
	if(!aStructure1.mWeights.empty())
	{
		double distance  = 0.0;
		double a_norm    = 0.0;
		double b_norm    = 0.0;
		for(unsigned int sect=0; sect < nsect1; ++sect)
		{
			float w1 = aStructure1.mWeights[sect];
			float w2 = aStructure2.mWeights[sect];

			for(i=0; i < sectlen1; ++i)
			{
				distance += aStructure1.mFingerprint[i+sectlen1*sect] * aStructure2.mFingerprint[i+sectlen1*sect] * w1 * w2;
				a_norm   += aStructure1.mFingerprint[i+sectlen1*sect] * aStructure1.mFingerprint[i+sectlen1*sect] * w1 * w1;
				b_norm   += aStructure2.mFingerprint[i+sectlen1*sect] * aStructure2.mFingerprint[i+sectlen1*sect] * w2 * w2;
			}
		}

		distance /= sqrt(a_norm*b_norm);
		return (float)((1. - distance)/2.);
	}

	// Distances per atom
	unsigned int nsect2 = aStructure2.mFingerprintNumSections;

	// Mark atoms already paired
	bool* a_atom_used = new bool[nsect1];
	bool* b_atom_used = new bool[nsect2];
	for(i=0; i < nsect1; ++i) a_atom_used[i] = false;
	for(i=0; i < nsect2; ++i) b_atom_used[i] = false;

	// For each section find the most similar section in the other structure
	double distance = 0.;
	for(i=0; i < nsect1; ++i)
	{
		// If already paired, skip it
		if(a_atom_used[i]) continue;

		double       curr_distance = DBL_MAX;
		unsigned int curr_min_idx  = UINT_MAX;
		unsigned int z1 = aStructure1.mAtomZ[i];

		for(unsigned int j=0; j < nsect2; ++j)
		{
			// If already paired, skip it
			if(b_atom_used[j]) continue;

			unsigned int z2 = aStructure2.mAtomZ[j];

			// If different type, skip it
			if(z1 != z2) continue;

			// Compute the distance
			double one_distance = 0.0;
			double a_norm       = 0.0;
			double b_norm       = 0.0;
			for(unsigned int k=0; k < sectlen1; ++k)
			{
				one_distance += aStructure1.mFingerprint[i*sectlen1+k] * aStructure2.mFingerprint[j*sectlen1+k];
				a_norm       += aStructure1.mFingerprint[i*sectlen1+k] * aStructure1.mFingerprint[i*sectlen1+k];
				b_norm       += aStructure2.mFingerprint[j*sectlen1+k] * aStructure2.mFingerprint[j*sectlen1+k];
			}
			one_distance /= sqrt(a_norm*b_norm);
			one_distance = (1. - one_distance)/2.;

			// Find the most similar atom in the other structure
			if(one_distance < curr_distance)
			{
				curr_distance = one_distance;
				curr_min_idx = j;
			}
		}

		// If a pairing has been found
		if(curr_min_idx < UINT_MAX)
		{
			distance += curr_distance;
			a_atom_used[i] = true;
			b_atom_used[curr_min_idx] = true;
		}
	}
	delete [] a_atom_used;
	delete [] b_atom_used;

	return (float)distance;
}


float EuclideanDistance::computeDistance(const Structure& aStructure1, const Structure& aStructure2) const
{
	unsigned int i;
	unsigned int nsect1 = aStructure1.mFingerprintNumSections;
	unsigned int sectlen1 = aStructure1.mFingerprintSectionLen;

	if(nsect1 == 1)
	{
		double distance = 0.;
		for(i=0; i < sectlen1; ++i)
		{
			double d = aStructure1.mFingerprint[i] - aStructure2.mFingerprint[i];
			distance += d*d;
		}

		return (float)sqrt(distance);
	}

	// More than one section, created by pseudo diffraction per atom
	if(!aStructure1.mWeights.empty())
	{
		double distance = 0.;
		for(unsigned int sect=0; sect < nsect1; ++sect)
		{
			double sdistance = 0.;
			for(i=0; i < sectlen1; ++i)
			{
				double d = aStructure1.mFingerprint[i+sectlen1*sect] - aStructure2.mFingerprint[i+sectlen1*sect];
				sdistance += d*d;
			}
			distance += sqrt(sdistance*aStructure1.mWeights[sect]*aStructure2.mWeights[sect]);
		}

		return (float)distance;
	}

	// Distances per atom
	unsigned int nsect2 = aStructure2.mFingerprintNumSections;
	unsigned int sectlen2 = aStructure2.mFingerprintSectionLen;

	// Mark atoms already paired
	bool* a_atom_used = new bool[nsect1];
	bool* b_atom_used = new bool[nsect2];
	for(i=0; i < nsect1; ++i) a_atom_used[i] = false;
	for(i=0; i < nsect2; ++i) b_atom_used[i] = false;

	// For each section find the most similar section in the other structure
	double distance = 0.;
	for(i=0; i < nsect1; ++i)
	{
		// If already paired, skip it
		if(a_atom_used[i]) continue;

		double       curr_distance = DBL_MAX;
		unsigned int curr_min_idx  = UINT_MAX;
		unsigned int z1 = aStructure1.mAtomZ[i];

		for(unsigned int j=0; j < nsect2; ++j)
		{
			// If already paired, skip it
			if(b_atom_used[j]) continue;

			unsigned int z2 = aStructure2.mAtomZ[j];

			// If different type, skip it
			if(z1 != z2) continue;

			// Compute the distance
			double one_distance = 0.0;
			for(unsigned int k=0; k < sectlen1; ++k)
			{
				double d = aStructure1.mFingerprint[i+sectlen1*k] - aStructure2.mFingerprint[j+sectlen2*k];
				one_distance += d*d;
			}
			one_distance = sqrt(one_distance);

			// Find the most similar atom in the other structure
			if(one_distance < curr_distance)
			{
				curr_distance = one_distance;
				curr_min_idx = j;
			}
		}

		// If a pairing has been found
		if(curr_min_idx < UINT_MAX)
		{
			distance += curr_distance;
			a_atom_used[i] = true;
			b_atom_used[curr_min_idx] = true;
		}
	}
	delete [] a_atom_used;
	delete [] b_atom_used;

	return (float)distance;
}

float MinkowskiDistance::computeDistance(const Structure& aStructure1, const Structure& aStructure2) const
{
	unsigned int i;
	unsigned int nsect1 = aStructure1.mFingerprintNumSections;
	unsigned int sectlen1 = aStructure1.mFingerprintSectionLen;

	if(nsect1 == 1)
	{
		double distance = 0.;
		for(i=0; i < sectlen1; ++i)
		{
			double d = aStructure1.mFingerprint[i] - aStructure2.mFingerprint[i];
			if(d < 0.) d = -d;
			distance += cbrt(d);
		}

		return (float)(distance*distance*distance);
	}

	// More than one section, created by pseudo diffraction per atom
	if(!aStructure1.mWeights.empty())
	{
		double distance = 0.;
		for(unsigned int sect=0; sect < nsect1; ++sect)
		{
			double sdistance = 0.;
			for(i=0; i < sectlen1; ++i)
			{
				double d = aStructure1.mFingerprint[i+sectlen1*sect] - aStructure2.mFingerprint[i+sectlen1*sect];
				if(d < 0.) d = -d;
				sdistance += cbrt(d);
			}
			distance += sdistance*sdistance*sdistance*sqrt(aStructure1.mWeights[sect]*aStructure2.mWeights[sect]);
		}

		return (float)distance;
	}

	// Distances per atom
	unsigned int nsect2 = aStructure2.mFingerprintNumSections;
	unsigned int sectlen2 = aStructure2.mFingerprintSectionLen;

	// Mark atoms already paired
	bool* a_atom_used = new bool[nsect1];
	bool* b_atom_used = new bool[nsect2];
	for(i=0; i < nsect1; ++i) a_atom_used[i] = false;
	for(i=0; i < nsect2; ++i) b_atom_used[i] = false;

	// For each section find the most similar section in the other structure
	double distance = 0.;
	for(i=0; i < nsect1; ++i)
	{
		// If already paired, skip it
		if(a_atom_used[i]) continue;

		double       curr_distance = DBL_MAX;
		unsigned int curr_min_idx  = UINT_MAX;
		unsigned int z1 = aStructure1.mAtomZ[i];

		for(unsigned int j=0; j < nsect2; ++j)
		{
			// If already paired, skip it
			if(b_atom_used[j]) continue;

			unsigned int z2 = aStructure2.mAtomZ[j];

			// If different type, skip it
			if(z1 != z2) continue;

			// Compute the distance
			double one_distance = 0.0;
			for(unsigned int k=0; k < sectlen1; ++k)
			{
				double d = aStructure1.mFingerprint[i+sectlen1*k] - aStructure2.mFingerprint[j+sectlen2*k];
				if(d < 0.) d = -d;
				one_distance += cbrt(d);
			}
			one_distance = one_distance*one_distance*one_distance;

			// Find the most similar atom in the other structure
			if(one_distance < curr_distance)
			{
				curr_distance = one_distance;
				curr_min_idx = j;
			}
		}

		// If a pairing has been found
		if(curr_min_idx < UINT_MAX)
		{
			distance += curr_distance;
			a_atom_used[i] = true;
			b_atom_used[curr_min_idx] = true;
		}
	}
	delete [] a_atom_used;
	delete [] b_atom_used;

	return (float)distance;
}

