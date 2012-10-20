
#include <cmath>
#include <cfloat>
#include <set>
#include "AnalysisMethod.h"
#include "SmoothPeak.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279
#endif

using namespace cfp_internal;

#ifdef _MSC_VER
namespace cfp_internal
{
extern double cbrt(double x);
}
#define cbrtf (float)cbrt
#else
#include <cmath>
#endif

void MethodGetIdx::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	for(unsigned int idx=0; idx < aCfp->getNumActiveStructures(); ++idx)
	{
		aValue[idx] = (float)idx;
	}
}


void MethodGetStep::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	for(unsigned int idx=0; idx < aCfp->getNumActiveStructures(); ++idx)
	{
		aValue[idx] = (float)aCfp->idxToStep(idx);
	}
}


void MethodGetEnergy::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	for(unsigned int idx=0; idx < aCfp->getNumActiveStructures(); ++idx)
	{
		aValue[idx] = aCfp->getPerAtomEnergy(idx);
	}
}


void MethodGetCellVolume::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	for(unsigned int idx=0; idx < aCfp->getNumActiveStructures(); ++idx)
	{
		const float* unit_cell = aCfp->getUnitCell(idx);
		float cell_volume =  unit_cell[0]*unit_cell[5]*unit_cell[10] + unit_cell[1]*unit_cell[6]*unit_cell[8] +
						     unit_cell[2]*unit_cell[4]*unit_cell[9]  - unit_cell[2]*unit_cell[5]*unit_cell[8] -
						     unit_cell[1]*unit_cell[4]*unit_cell[10] - unit_cell[0]*unit_cell[6]*unit_cell[9];
							 
		aValue[idx] = cell_volume/aCfp->getNatoms(idx);
	}
}


void MethodGetDeltaEnergy::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	unsigned int idx = 0;
	for(unsigned int i1=0; i1 < aCfp->getNumActiveStructures()-1; ++i1)
	{
		for(unsigned int i2=i1+1; i2 < aCfp->getNumActiveStructures(); ++i2)
		{
			aValue[idx++] = fabsf(aCfp->getPerAtomEnergy(i1)-aCfp->getPerAtomEnergy(i2));
		}
	}
}


void MethodGetDistances::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	size_t idx = 0;
	for(size_t i1=0; i1 < aCfp->getNumActiveStructures()-1; ++i1)
	{
		for(size_t i2=i1+1; i2 < aCfp->getNumActiveStructures(); ++i2)
		{
			aValue[idx++] = aCfp->getDistance(i1, i2);
		}
	}
}


void MethodGetEnergyFromMin::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	unsigned int idx;
	float min_energy = FLT_MAX;
	for(idx=0; idx < aCfp->getNumActiveStructures(); ++idx)
	{
		float e = aCfp->getPerAtomEnergy(idx);
		if(e < min_energy) min_energy = e;
	}
	for(idx=0; idx < aCfp->getNumActiveStructures(); ++idx)
	{
		aValue[idx] = aCfp->getPerAtomEnergy(idx) - min_energy;
	}
}


void MethodGetDistFromMin::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	size_t idx;
	float min_energy = FLT_MAX;
	size_t min_idx   = 0;
	for(idx=0; idx < aCfp->getNumActiveStructures(); ++idx)
	{
		float e = aCfp->getPerAtomEnergy(idx);
		if(e < min_energy) {min_energy = e; min_idx = idx;}
	}
	for(idx=0; idx < aCfp->getNumActiveStructures(); ++idx)
	{
		aValue[idx] = aCfp->getDistance(idx, min_idx);
	}
}

void MethodGetPointDepth::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	//std::vector<float> v;
	unsigned int k;
	size_t i, j;

	unsigned int len = aCfp->getFingerprintNumSections() * aCfp->getFingerprintSectionLen();
	float *vect = new float[len];
	float *e = new float[len];
	size_t ns = aCfp->getNumActiveStructures();
#ifdef SPHERICAL_DEPTH
	for(i=0; i < ns; ++i)
	{
		const float* fpi = aCfp->getFingerprint(i);

		float sqlen = 0;
		for(k=0; k < len; ++k) sqlen += fpi[k]*fpi[k];
		for(k=0; k < len; ++k) e[k] = 0;

		for(j=0; j < ns; ++j)
		{
			if(i == j) continue;

			const float* fpj = aCfp->getFingerprint(j);

			float dot = 0;
			for(k=0; k < len; ++k) dot += fpi[k]*fpj[k];
			dot /= sqlen;

			for(k=0; k < len; ++k) vect[k] = fpj[k] - dot*fpi[k];

			float vectlen = 0;
			for(k=0; k < len; ++k) vectlen += vect[k]*vect[k];
			vectlen = (float)sqrt(vectlen);

			for(k=0; k < len; ++k) e[k] += vect[k]/vectlen;
		}

		float elen = 0;
		for(k=0; k < len; ++k) elen += e[k]*e[k];
		elen = (float)sqrt(elen);

		float x = elen/ns-1.0F/ns;
		if(x < 0) x = 0.0F;
		float depth = 1.0F - x;
		aValue[i] = depth;
	}
#else
	for(i=0; i < ns; ++i)
	{
		const float* fpi = aCfp->getFingerprint(i);

		for(k=0; k < len; ++k) e[k] = 0;

		for(j=0; j < ns; ++j)
		{
			if(i == j) continue;

			const float* fpj = aCfp->getFingerprint(j);

			for(k=0; k < len; ++k) vect[k] = fpi[k] - fpj[k];
			float vectlen = 0;
			for(k=0; k < len; ++k) vectlen += vect[k]*vect[k];
			vectlen = (float)sqrt(vectlen);

			for(k=0; k < len; ++k) e[k] += vect[k]/vectlen;
		}

		float elen = 0;
		for(k=0; k < len; ++k) elen += e[k]*e[k];
		elen = (float)sqrt(elen);

		float x = elen/ns-1.0F/ns;
		if(x < 0) x = 0.0F;
		float depth = 1.0F - x;
		aValue[i] = depth;
	}
#endif
	delete [] vect;
	delete [] e;
}


void MethodGetOrderF2::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	unsigned int i, j;

	unsigned int ns = aCfp->getFingerprintNumSections();
	if(ns == 1)
	{
		// Get the fingerprint length
		unsigned int fplen = aCfp->getFingerprintSectionLen();

		for(i=0; i < aCfp->getNumActiveStructures(); ++i)
		{
			// Extract fingerprint values
			const float *yv = aCfp->getFingerprint(i);

			// Compute V^1/3
			const float* unit_cell = aCfp->getUnitCell(i);
			float Vuc = unit_cell[0]*unit_cell[5]*unit_cell[10] + unit_cell[1]*unit_cell[6]*unit_cell[8] +
						unit_cell[2]*unit_cell[4]*unit_cell[9]  - unit_cell[2]*unit_cell[5]*unit_cell[8] -
						unit_cell[1]*unit_cell[4]*unit_cell[10] - unit_cell[0]*unit_cell[6]*unit_cell[9];
			float R0 = cbrtf(Vuc);

			// Compute the degree of order
			float op = 0.0F;
			for(j=0; j < fplen; ++j) op += yv[j]*yv[j];
			aValue[i] = aCfp->getDiffrBinSize()*op/R0;
		}
	}
	else
	{
		// Get the fingerprint length
		unsigned int fplen = aCfp->getFingerprintSectionLen();

		for(i=0; i < aCfp->getNumActiveStructures(); ++i)
		{
			float degree_of_order = 0.0F;
			const float* Wij = aCfp->getWeights(i);

			// Compute V^1/3
			const float* unit_cell = aCfp->getUnitCell(i);
			float Vuc = unit_cell[0]*unit_cell[5]*unit_cell[10] + unit_cell[1]*unit_cell[6]*unit_cell[8] +
						unit_cell[2]*unit_cell[4]*unit_cell[9]  - unit_cell[2]*unit_cell[5]*unit_cell[8] -
						unit_cell[1]*unit_cell[4]*unit_cell[10] - unit_cell[0]*unit_cell[6]*unit_cell[9];
			float R0 = cbrtf(Vuc);

			for(j=0; j < ns; ++j)
			{
				// Extract fingerprint values
				const float *yv = aCfp->getFingerprint(i)+j*aCfp->getFingerprintSectionLen();

				// Compute the degree of order
				float op = 0.0F;
				for(unsigned int k=0; k < fplen; ++k) op += yv[k]*yv[k];
				degree_of_order += aCfp->getDiffrBinSize()*op/R0*Wij[j];
			}
			aValue[i] = degree_of_order;
		}
	}
}

void MethodGetOrderF2R2::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	unsigned int i, j;

	unsigned int ns = aCfp->getFingerprintNumSections();
	if(ns == 1)
	{
		// Get the fingerprint length
		unsigned int fplen = aCfp->getFingerprintSectionLen();

		for(i=0; i < aCfp->getNumActiveStructures(); ++i)
		{
			// Compute Nuc/Vuc
			float Nuc = (float)aCfp->getNatoms(i);
			const float* unit_cell = aCfp->getUnitCell(i);
			float Vuc = unit_cell[0]*unit_cell[5]*unit_cell[10] + unit_cell[1]*unit_cell[6]*unit_cell[8] +
						unit_cell[2]*unit_cell[4]*unit_cell[9]  - unit_cell[2]*unit_cell[5]*unit_cell[8] -
						unit_cell[1]*unit_cell[4]*unit_cell[10] - unit_cell[0]*unit_cell[6]*unit_cell[9];

			// Extract fingerprint values
			const float *yv = aCfp->getFingerprint(i);

			// Compute the degree of order
			float op = 0.0F;
			for(j=0; j < fplen; ++j) op += yv[j]*yv[j]*(j*aCfp->getDiffrBinSize())*(j*aCfp->getDiffrBinSize());
			float degree_of_order = aCfp->getDiffrBinSize()*op*Nuc/Vuc;
			aValue[i] = degree_of_order;
		}
	}
	else
	{
		// Get the fingerprint length
		unsigned int fplen = aCfp->getFingerprintSectionLen();

		for(i=0; i < aCfp->getNumActiveStructures(); ++i)
		{
			// Compute Nuc/Vuc
			float Nuc = (float)aCfp->getNatoms(i);
			const float* unit_cell = aCfp->getUnitCell(i);
			float Vuc = unit_cell[0]*unit_cell[5]*unit_cell[10] + unit_cell[1]*unit_cell[6]*unit_cell[8] +
						unit_cell[2]*unit_cell[4]*unit_cell[9]  - unit_cell[2]*unit_cell[5]*unit_cell[8] -
						unit_cell[1]*unit_cell[4]*unit_cell[10] - unit_cell[0]*unit_cell[6]*unit_cell[9];

			float degree_of_order = 0.0F;
			const float* Wij = aCfp->getWeights(i);

			for(j=0; j < ns; ++j)
			{
				// Extract fingerprint values
				const float *yv = aCfp->getFingerprint(i)+j*aCfp->getFingerprintSectionLen();

				// Compute the degree of order
				float op = 0.0F;
				for(unsigned int k=0; k < fplen; ++k) op += yv[k]*yv[k]*(k*aCfp->getDiffrBinSize())*(k*aCfp->getDiffrBinSize());
				degree_of_order += aCfp->getDiffrBinSize()*op*Wij[j];
			}

			degree_of_order *= Nuc/Vuc;
			aValue[i] = degree_of_order;
		}
	}
}



void MethodGetQuasiEntropy::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	for(unsigned int idx=0; idx < aCfp->getNumActiveStructures(); ++idx)
	{
		aValue[idx] = computeQuasiEntropy(aCfp, idx);
	}
}

float MethodGetQuasiEntropy::computeQuasiEntropy(const cfp::CrystalFp* aCfp, unsigned int aIdx) const
{
	unsigned int i;
	unsigned int ai, aj;

	// Prepare list of different Z values
	std::set<unsigned int> different_z;
	const unsigned int* z = aCfp->getAtomZ(aIdx);
	unsigned int natoms = aCfp->getNatoms(aIdx);
	for(i=0; i < natoms; ++i) different_z.insert(z[i]);

	// Compute the volume of the unit cell
	const float *unit_cell = aCfp->getUnitCell(aIdx);
	float Vuc = unit_cell[0]*unit_cell[5]*unit_cell[10] + unit_cell[1]*unit_cell[6]*unit_cell[8] +
				unit_cell[2]*unit_cell[4]*unit_cell[9]  - unit_cell[2]*unit_cell[5]*unit_cell[8] -
				unit_cell[1]*unit_cell[4]*unit_cell[10] - unit_cell[0]*unit_cell[6]*unit_cell[9];

	// Gather data for the single atom fingerprint computation
	float diffraction_bin_size = aCfp->getDiffrBinSize();
	float diffraction_peak_size = aCfp->getDiffrPeakSize();
	float max_distance = aCfp->getCutoffDistance();
	unsigned int nbins = (int)(max_distance / diffraction_bin_size + 0.5F);
	float delta = max_distance/nbins;

	// Allocate the one-atom fingerprints and initialize to -1 (to normalize the result)
	float *fp = new float[nbins*natoms];
	for(i=0; i < nbins*natoms; ++i) fp[i] = -1.0F;

	// Compute all single-atom fingerprints
	const std::vector< std::vector<float> >& idist = aCfp->getInteratomicDistances(aIdx);
	for(ai=0; ai < natoms; ++ai)
	{
		for(i=0; i < idist[ai].size(); ++i)
		{
			float Rij = idist[ai][i];
			float Fing = (float)(1./(4.*M_PI*Rij*Rij*(natoms/Vuc)*diffraction_bin_size));
			smoothPeak(Fing, Rij, delta, nbins, fp+ai*nbins, diffraction_peak_size);
		}
	}

#if 0
	printf("\nComputeOneQuasiEntropy %d\n", idx);
	for(i=0; i < nbins; ++i) printf(" %g", fp[i]);
	printf("\n");
#endif

	// For each atomic specie
	float quasi_entropy = 0.0F;
	std::set<unsigned int>::const_iterator idz;
	for(idz=different_z.begin(); idz != different_z.end(); ++idz)
	{
		// Number of atoms of this specie
		float Na = 0.0F;
		for(i=0; i < natoms; ++i) if(z[i] == *idz) ++Na;

		// Compute collective diversity
		double coll_div  = 0.;
		double ncoll_div = 0.;
		for(ai=0; ai < natoms-1; ++ai)
		{
			if(z[ai] != *idz) continue;

			for(aj=ai+1; aj < natoms; ++aj)
			{
				if(z[aj] != *idz) continue;

				// Compute the distance between the two one-atom fingerprints
				const float *afp = fp + ai*nbins;
				const float *bfp = fp + aj*nbins;
				double Daiaj     = 0.;
				double a_norm    = 0.;
				double b_norm    = 0.;

				for(i=0; i < nbins; ++i)
				{
					Daiaj  += afp[i] * bfp[i];
					a_norm += afp[i] * afp[i];
					b_norm += bfp[i] * bfp[i];
				}
				Daiaj /= sqrt(a_norm*b_norm);
				Daiaj = (1. - Daiaj)/2.;

				// Accumulate collective diversity
				coll_div  += (1.-Daiaj)*log(1.-Daiaj);
				ncoll_div += 1.;
			}
		}

		// Sum over A (only if collective diversity can be computed)
		if(ncoll_div != 0.) quasi_entropy += -Na/(float)natoms*(float)(coll_div/ncoll_div);
	}

	delete [] fp;
	return quasi_entropy;
}

void MethodGetMinDimension::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int /*aIdx*/) const
{
	unsigned int seclen = aCfp->getFingerprintSectionLen();
	unsigned int numsec = aCfp->getFingerprintNumSections();

	for(unsigned int idx=0; idx < aCfp->getNumActiveStructures(); ++idx)
	{
		const float* fp = aCfp->getFingerprint(idx);
		unsigned int i;
		for(i=0; i < seclen; ++i)
		{
			bool is_at_min = true;
			for(unsigned int part=0; part < numsec; ++part)
			{
				if(fp[part*seclen+i] > -1.0F+1e-7) {is_at_min = false; break;}
			}
			if(!is_at_min) break;
		}

		aValue[idx] = (float)((seclen-i)*numsec);
	}
}


const std::vector<std::string> MethodGetFingerprint::getLabels(const cfp::CrystalFp* aCfp) const
{
	std::vector<std::string> lbl;

	if(aCfp->isDiffractionLike())
	{
		lbl.push_back("Distance");
		lbl.push_back("Intensity");
	}
	else
	{
		lbl.push_back("Index");
		lbl.push_back("Distance");
	}

	return lbl;
}


void MethodGetFingerprint::getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx) const
{
	unsigned int i;
	unsigned int len;

	switch(aIdx)
	{
	case 1:
		if(mPartSelected >= aCfp->getFingerprintNumSections())
		{
			len = aCfp->getFingerprintSectionLen()*aCfp->getFingerprintNumSections();
			memcpy(aValue, aCfp->getFingerprint(mStructureIdx), len*sizeof(float));
		}
		else
		{
			len = aCfp->getFingerprintSectionLen();
			memcpy(aValue, aCfp->getFingerprint(mStructureIdx)+mPartSelected*aCfp->getFingerprintSectionLen(), len*sizeof(float));
		}
		break;

	case 0:
		len = aCfp->getFingerprintSectionLen()*((mPartSelected >= aCfp->getFingerprintNumSections()) ? aCfp->getFingerprintNumSections() : 1);

		if(aCfp->isDiffractionLike())
		{
			for(i=0; i < len; ++i) aValue[i] = i * aCfp->getDiffrPeakSize();
		}
		else
		{
			for(i=0; i < len; ++i) aValue[i] = (float)i;
		}
		break;
	}
}


