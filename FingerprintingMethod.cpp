
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "FingerprintingMethod.h"
#include "SmoothPeak.h"

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

using namespace cfp_internal;

#if 0
void NormalizedDiffractionHistogram::computeFingerprint(Structure& aStructure, const unsigned int* aExpansion)
{
	// Compute fingerprint sizes
	const unsigned int MAX_Z = 111;
	int num_atoms[MAX_Z+1];
	unsigned int i, j;
	for(i=1; i <= MAX_Z; ++i) num_atoms[i] = 0;
	std::vector<unsigned int>::const_iterator iz;
	for(iz=aStructure.mAtomZ.begin(); iz != aStructure.mAtomZ.end(); ++iz) ++num_atoms[*iz];
	unsigned int num_species = 0;
	for(i=1; i <= MAX_Z; ++i) if(num_atoms[i]) ++num_species;
	int atoms_idx[MAX_Z+1];
	int pos = 0;
	for(i=1; i <= MAX_Z; ++i) if(num_atoms[i]) atoms_idx[i] = pos++;

	unsigned int num_sections = (num_species*(num_species+1))/2;
	unsigned int nbins = (unsigned int)((mCutoff/mBinSize)+0.5F);
	unsigned int fp_len = nbins*num_sections;
	float delta = mCutoff/nbins;
    
	float* unit_cell = 0;
	float cell_volume = 0;
	if(!mIsNanocluster)
	{
		unit_cell = aStructure.mUnitCell;
		cell_volume = unit_cell[0]*unit_cell[5]*unit_cell[10] + unit_cell[1]*unit_cell[6]*unit_cell[8] +
					  unit_cell[2]*unit_cell[4]*unit_cell[9]  - unit_cell[2]*unit_cell[5]*unit_cell[8] -
					  unit_cell[1]*unit_cell[4]*unit_cell[10] - unit_cell[0]*unit_cell[6]*unit_cell[9];
	}

	// Create the table to decode the fused loop
	int ex =  static_cast<int>(aExpansion[0]);
	int ey =  static_cast<int>(aExpansion[1]);
	int ez =  static_cast<int>(aExpansion[2]);
	int n = 0;
	int ii;
	int di, dj, dk;
	int iorig;
	unsigned int imax = (2*aExpansion[0]+1)*(2*aExpansion[1]+1)*(2*aExpansion[2]+1)*3;
	int *d = new int[imax];
	for(di= -ex; di <= ex; ++di)
	{
		for(dj= -ey; dj <= ey; ++dj)
		{
			for(dk = -ez; dk <= ez; ++dk)
			{
				d[n++] = di;
				d[n++] = dj;
				d[n++] = dk;

				if(di == 0 && dj == 0 && dk == 0) iorig = n-3;
			}
		}
	}
	
	// Create the infinite slab
	int natoms = aStructure.mNumAtoms;
	const float* coords = &(aStructure.mCoordinates[0]);
	const unsigned int* atom_z = &(aStructure.mAtomZ[0]);

	// Prepare the accumulator for the parallel version
#ifdef _OPENMP
	int nt = omp_get_max_threads();
#else
	int nt = 1;
#endif

	// Every thread update a line of this temporary array
	float *fpp = new float[fp_len*nt];
	memset(fpp, 0, fp_len*nt*sizeof(float));

	#pragma omp parallel for private(ii) default(none) shared(iorig, natoms, imax, atom_z, fpp, d, coords, num_atoms, cell_volume, num_species, unit_cell, delta, nbins, fp_len, atoms_idx)
	for(ii=0; ii < (int)imax; ii += 3)
	{
#ifdef _OPENMP
		int ct = omp_get_thread_num();
#else
		int ct = 0;
#endif
		// Copy the atoms in the unit cell replicas
		for(int a=0; a < natoms; ++a)
		{
			float x = coords[3*a+0];
			float y = coords[3*a+1];
			float z = coords[3*a+2];

			int Zi = atom_z[a];
			int Ni = num_atoms[Zi];
			int Pi = atoms_idx[Zi];

			if(ii == iorig)
			{
				for(int b=a+1; b < natoms; ++b)
				{
					float dx = coords[3*b+0] - x;
					float dy = coords[3*b+1] - y;
					float dz = coords[3*b+2] - z;

					float dist_squared = dx*dx+dy*dy+dz*dz;
					int Zj = atom_z[b];
					int Nj = num_atoms[Zj];
					int Pj = atoms_idx[Zj];

					// Compute the peak value Fing
					float fing = (float)((mIsNanocluster) ? 1./(Nj*Ni*mBinSize) : 1./(4.*M_PI*dist_squared*(Nj/cell_volume)*2.*Ni*mBinSize));
			//float Fing = (float)((is_cluster) ? (Zi*Zj)/(Nuc*diffraction_bin_size) : (Zi*Zj)/(4.*M_PI*Rij*Rij*Nuc/Vuc*diffraction_bin_size));

					// The components AA, BB, etc. should be counted twice
					if(Zi == Zj) fing *= 2.0F;

					// Compute the section index for this part
					int idx_section = (Pj >= Pi) ? Pi*num_species-(Pi*(Pi+1))/2+Pj : Pj*num_species-(Pj*(Pj+1))/2+Pi;

					// Smooth the peak and accumulate
					smoothPeak(fing, sqrt(dist_squared), delta, nbins, fpp + (ct*fp_len) + idx_section*nbins, mPeakSize);
				}
			}
			else
			{
				int di = d[ii+0];
				int dj = d[ii+1];
				int dk = d[ii+2];

				float ox = x + di*unit_cell[0] + dj*unit_cell[4] + dk*unit_cell[8];
				float oy = y + di*unit_cell[1] + dj*unit_cell[5] + dk*unit_cell[9];
				float oz = z + di*unit_cell[2] + dj*unit_cell[6] + dk*unit_cell[10];

				for(int b=0; b < natoms; ++b)
				{
					float dx = coords[3*b+0] - ox;
					float dy = coords[3*b+1] - oy;
					float dz = coords[3*b+2] - oz;

					float dist_squared = dx*dx+dy*dy+dz*dz;
					int Zj = atom_z[b];
					int Nj = num_atoms[Zj];
					int Pj = atoms_idx[Zj];

					// Compute the peak value Fing
					float fing = (float)((mIsNanocluster) ? 1./(Nj*Ni*mBinSize) : 1./(4.*M_PI*dist_squared*(Nj/cell_volume)*2.*Ni*mBinSize));

					// The components AA, BB, etc. should be counted twice
					if(Zi == Zj) fing *= 2.0F;

					// Compute the section index for this part
					int idx_section = (Pj >= Pi) ? Pi*num_species-(Pi*(Pi+1))/2+Pj : Pj*num_species-(Pj*(Pj+1))/2+Pi;

					// Smooth the peak and accumulate
					smoothPeak(fing, sqrt(dist_squared), delta, nbins, fpp + (ct*fp_len) + idx_section*nbins, mPeakSize);
				}
			}
		}
	}

	// Normalize and accumulate the per-thread fingerprints
	aStructure.mFingerprintNumSections = num_sections;
	aStructure.mFingerprintSectionLen = nbins;
	aStructure.mFingerprint.assign(fp_len, -1.0F);

	#pragma omp parallel for private(n, ii) shared(fpp, fp_len, nt)
	for(ii=0; ii < (int)fp_len; ++ii)
	{
		for(n=0; n < nt; ++n)
		{
			aStructure.mFingerprint[ii] += fpp[n*fp_len+ii];
		}
	}

	delete [] fpp;
	delete [] d;

	// Compute weights
	for(i=1; i <= MAX_Z; ++i)
	{
		if(num_atoms[i] == 0) continue;
		for(j=i; j <= MAX_Z; ++j)
		{
			if(num_atoms[j] == 0) continue;

			float w = (float)(num_atoms[i]*num_atoms[j]);
			aStructure.mWeights.push_back(w);
		}
	}

	// Normalize the weights and store them
	float w = 0.0F;
	for(i=0; i < num_sections; ++i) w += aStructure.mWeights[i];
	for(i=0; i < num_sections; ++i) aStructure.mWeights[i] /= w;
}
#endif


void PerElementRdfHistogram::computeFingerprint(Structure& aStructure, const unsigned int* aExpansion)
{
	std::cerr << "Start" << std::endl;
	// Compute fingerprint sizes
	const unsigned int MAX_Z = 111;
	int num_atoms[MAX_Z+1];
	unsigned int i, j;
	//for(i=1; i <= MAX_Z; ++i) num_atoms[i] = 0;
	memset(num_atoms, 0, (MAX_Z+1)*sizeof(int));
	std::vector<unsigned int>::const_iterator iz;
	for(iz=aStructure.mAtomZ.begin(); iz != aStructure.mAtomZ.end(); ++iz) ++num_atoms[*iz];
	unsigned int num_species = 0;
	for(i=1; i <= MAX_Z; ++i) if(num_atoms[i]) ++num_species;
	int atoms_idx[MAX_Z+1];
	int pos = 0;
	for(i=1; i <= MAX_Z; ++i) if(num_atoms[i]) atoms_idx[i] = pos++;

	unsigned int num_sections = (num_species*(num_species+1))/2;
	unsigned int nbins = (unsigned int)((mCutoff/mBinSize)+0.5F);
	size_t fp_len = nbins*num_sections;
	float delta = mCutoff/nbins;
    
	float* unit_cell = 0;
	float cell_volume = 0;
	if(!mIsNanocluster)
	{
		unit_cell = aStructure.mUnitCell;
		cell_volume = unit_cell[0]*unit_cell[5]*unit_cell[10] + unit_cell[1]*unit_cell[6]*unit_cell[8] +
					  unit_cell[2]*unit_cell[4]*unit_cell[9]  - unit_cell[2]*unit_cell[5]*unit_cell[8] -
					  unit_cell[1]*unit_cell[4]*unit_cell[10] - unit_cell[0]*unit_cell[6]*unit_cell[9];
	}

	// Create the table to decode the fused loop
	int ex =  static_cast<int>(aExpansion[0]);
	int ey =  static_cast<int>(aExpansion[1]);
	int ez =  static_cast<int>(aExpansion[2]);
	int n = 0;
	int ii;
	int di, dj, dk;
	int iorig;
	size_t imax = (2*aExpansion[0]+1)*(2*aExpansion[1]+1)*(2*aExpansion[2]+1)*3;
	int *d = new int[imax];
	for(di= -ex; di <= ex; ++di)
	{
		for(dj= -ey; dj <= ey; ++dj)
		{
			for(dk = -ez; dk <= ez; ++dk)
			{
				d[n++] = di;
				d[n++] = dj;
				d[n++] = dk;

				if(di == 0 && dj == 0 && dk == 0) iorig = n-3;
			}
		}
	}
	
	// Create the infinite slab
	int natoms = aStructure.mNumAtoms;
	const float* coords = &(aStructure.mCoordinates[0]);
	const unsigned int* atom_z = &(aStructure.mAtomZ[0]);

	// Prepare the interatomic distances array
	aStructure.mInteratomicDistances.clear();
	for(ii=0; ii < natoms; ++ii) aStructure.mInteratomicDistances.push_back(std::vector<float>());
	std::vector< std::vector< std::pair<int, float> > >dp;
	for(ii=0; ii < (int)imax/3; ++ii) dp.push_back(std::vector< std::pair<int, float> >());

	// Prepare the accumulator for the parallel version
#ifdef _OPENMP
	int nt = omp_get_max_threads();
#else
	int nt = 1;
#endif
		std::cerr << "Start2" << std::endl;

	// Every thread update a line of this temporary array
	float *fpp = new float[fp_len*nt];
	memset(fpp, 0, fp_len*nt*sizeof(float));
		std::cerr << "Start3" << std::endl;

	#pragma omp parallel for private(ii) default(none) shared(iorig, natoms, imax, atom_z, fpp, d, coords, num_atoms, cell_volume, num_species, unit_cell, delta, nbins, fp_len, atoms_idx, dp)
	for(ii=0; ii < (int)imax; ii += 3)
	{
#ifdef _OPENMP
		int ct = omp_get_thread_num();
#else
		int ct = 0;
#endif
		// Copy the atoms in the unit cell replicas
		for(int a=0; a < natoms; ++a)
		{
			float x = coords[3*a+0];
			float y = coords[3*a+1];
			float z = coords[3*a+2];

			int Zi = atom_z[a];
			int Ni = num_atoms[Zi];
			int Pi = atoms_idx[Zi];

			if(ii == iorig)
			{
				for(int b=a+1; b < natoms; ++b)
				{
					float dx = coords[3*b+0] - x;
					float dy = coords[3*b+1] - y;
					float dz = coords[3*b+2] - z;

					float dist_squared = dx*dx+dy*dy+dz*dz;
					int Zj = atom_z[b];
					int Nj = num_atoms[Zj];
					int Pj = atoms_idx[Zj];

					// Compute the peak value Fing
					float fing = (float)((mIsNanocluster) ? 1./(Nj*Ni*mBinSize) : 1./(4.*M_PI*dist_squared*(Nj/cell_volume)*2.*Ni*mBinSize));

					// The components AA, BB, etc. should be counted twice
					if(Zi == Zj) fing *= 2.0F;

					// Compute the section index for this part
					int idx_section = (Pj >= Pi) ? Pi*num_species-(Pi*(Pi+1))/2+Pj : Pj*num_species-(Pj*(Pj+1))/2+Pi;

					// Smooth the peak and accumulate
					float dist = (float)sqrt(dist_squared);
					smoothPeak(fing, dist, delta, nbins, fpp + (ct*fp_len) + idx_section*nbins, mPeakSize);

					// Save the interatomic distance
					dp[ii/3].push_back(std::pair<int, float>(a, dist));
				}
			}
			else
			{
				int di = d[ii+0];
				int dj = d[ii+1];
				int dk = d[ii+2];

				float ox = x + di*unit_cell[0] + dj*unit_cell[4] + dk*unit_cell[8];
				float oy = y + di*unit_cell[1] + dj*unit_cell[5] + dk*unit_cell[9];
				float oz = z + di*unit_cell[2] + dj*unit_cell[6] + dk*unit_cell[10];

				for(int b=0; b < natoms; ++b)
				{
					float dx = coords[3*b+0] - ox;
					float dy = coords[3*b+1] - oy;
					float dz = coords[3*b+2] - oz;

					float dist_squared = dx*dx+dy*dy+dz*dz;
					int Zj = atom_z[b];
					int Nj = num_atoms[Zj];
					int Pj = atoms_idx[Zj];

					// Compute the peak value Fing
					float fing = (float)((mIsNanocluster) ? 1./(Nj*Ni*mBinSize) : 1./(4.*M_PI*dist_squared*(Nj/cell_volume)*2.*Ni*mBinSize));

					// The components AA, BB, etc. should be counted twice
					if(Zi == Zj) fing *= 2.0F;

					// Compute the section index for this part
					int idx_section = (Pj >= Pi) ? Pi*num_species-(Pi*(Pi+1))/2+Pj : Pj*num_species-(Pj*(Pj+1))/2+Pi;

					// Smooth the peak and accumulate
					float dist = (float)sqrt(dist_squared);
					smoothPeak(fing, dist, delta, nbins, fpp + (ct*fp_len) + idx_section*nbins, mPeakSize);

					// Save the interatomic distance
					dp[ii/3].push_back(std::pair<int, float>(a, dist));
				}
			}
		}
	}
		std::cerr << "Start4" << std::endl;

	// Normalize and accumulate the per-thread fingerprints
	aStructure.mFingerprintNumSections = num_sections;
	aStructure.mFingerprintSectionLen = nbins;
	aStructure.mFingerprint.assign(fp_len, -1.0F);

	#pragma omp parallel for private(n, ii) shared(fpp, fp_len, nt)
	for(ii=0; ii < (int)fp_len; ++ii)
	{
		for(n=0; n < nt; ++n)
		{
			aStructure.mFingerprint[ii] += fpp[n*fp_len+ii];
		}
	}

	delete [] fpp;
	delete [] d;

	// Compute weights
	for(i=1; i <= MAX_Z; ++i)
	{
		if(num_atoms[i] == 0) continue;
		for(j=i; j <= MAX_Z; ++j)
		{
			if(num_atoms[j] == 0) continue;

			float w = (float)(num_atoms[i]*num_atoms[j]);
			aStructure.mWeights.push_back(w);
		}
	}

	// Normalize the weights and store them
	float w = 0.0F;
	for(i=0; i < num_sections; ++i) w += aStructure.mWeights[i];
	for(i=0; i < num_sections; ++i) aStructure.mWeights[i] /= w;

	// Store the interatomic distances
	for(ii=0; ii < (int)imax/3; ++ii)
	{
		std::vector< std::pair<int, float> >::const_iterator idp;
		for(idp=dp[ii].begin(); idp != dp[ii].end(); ++idp)
		{
			aStructure.mInteratomicDistances[idp->first].push_back(idp->second);
		}
	}
		std::cerr << "End" << std::endl;
}
