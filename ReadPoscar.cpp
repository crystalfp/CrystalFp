
#include <cstdlib>
#include <cstdio>
#include "CrystalFp.h"
#include "ReadPoscar.h"

static bool loadData(int aStep, FILE *aPoscarFp, FILE *aEnergyFp, bool aEnergyIsPerAtom, std::vector<unsigned int>& aAtomZ, CrystalFp& aCfp)
{
	char linebuf[82];
	unsigned int k;
	unsigned int j;
	float fx;
	float fy;
	float fz;
	float scale;
	float uc[16];

	// Ignore the title
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;

	// Read the scale factor
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;
	sscanf(linebuf, "%f", &scale);

	// Read the unit cell base vectors
	for(j=0; j < 3; ++j)
	{
		fgets(linebuf, sizeof linebuf, aPoscarFp);
		sscanf(linebuf, "%f %f %f", &uc[4*j+0], &uc[4*j+1], &uc[4*j+2]);
		uc[4*j+0] *= scale;
		uc[4*j+1] *= scale;
		uc[4*j+2] *= scale;
		uc[4*j+3]  = 0.0F;
	}
	uc[12+0]  = 0.0F;
	uc[12+1]  = 0.0F;
	uc[12+2]  = 0.0F;
	uc[12+3]  = 1.0F;

	// Read the array of atoms numbers
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;
	long int res;
	char *endptr;
	const char *next = linebuf;
	unsigned int natoms = 0;
	unsigned int ntypes = 0;
	for(;;)
	{
		res = strtol(next, &endptr, 10);
		if(endptr == next) break;

		natoms += res;
		++ntypes;

		next = endptr;
	}

	int *natomtype = new int[ntypes];
	next = linebuf;
	for(k=0;;++k)
	{
		res = strtol(next, &endptr, 10);
		if(endptr == next) break;
		natomtype[k] = res;

		next = endptr;
	}

	// Skip the "Selective Dynamics" line
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) {delete [] natomtype; return false;}
	if(linebuf[0] == 'S' || linebuf[0] == 's')
	{
		if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) {delete [] natomtype; return false;}
	}

	// Set the kind of coordinates (cartesian or fractional)
	bool cartesian = false;
	if(linebuf[0] == 'C' || linebuf[0] == 'c' || linebuf[0] == 'K' || linebuf[0] == 'k')
	{
		cartesian = true;
	}
	else if(linebuf[0] != 'D' && linebuf[0] != 'd')
	{
		delete [] natomtype;
		return false;
	}

	// Sanity check
	if(natoms < 1)
	{
		delete [] natomtype;
		return false;
	}

	// Allocate the various array for AddStructure
	unsigned int *z = new unsigned int[natoms];
	float *coords = new float[3*natoms];

	// Read the atom records
	unsigned int currbnd = natomtype[0];
	for(j=k=0; j < natoms; ++j)
	{
		fgets(linebuf, sizeof linebuf, aPoscarFp);
		if(sscanf(linebuf, "%f %f %f", &fx, &fy, &fz) < 3)
		{
			delete [] natomtype;
			delete [] z;
			delete [] coords;
			return false;
		}

		// Change atom type
		if(j >= currbnd)
		{
			++k;
			currbnd += natomtype[k];
		}

		// Add the atom converting coordinates from fractional
		z[j] = (aAtomZ.size() > k) ? aAtomZ[k] : k+1;
		if(cartesian)
		{
			coords[3*j+0] = fx;
			coords[3*j+1] = fy;
			coords[3*j+2] = fz;
		}
		else
		{
			coords[3*j+0] = fx*uc[0]+fy*uc[4]+fz*uc[8];
			coords[3*j+1] = fx*uc[1]+fy*uc[5]+fz*uc[9];
			coords[3*j+2] = fx*uc[2]+fy*uc[6]+fz*uc[10];
		}
	}

	// Read the energy, if present and add the structure
	if(aEnergyFp)
	{
		float energy = 0.0F;
		fgets(linebuf, sizeof linebuf, aEnergyFp);
		sscanf(linebuf, "%f", &energy);

		aCfp.addStructureBatch(aStep, natoms, coords, z, uc, true, energy, aEnergyIsPerAtom);
	}
	else
	{
		aCfp.addStructureBatch(aStep, natoms, coords, z, uc, false, 0.0F, false);
	}

	delete [] natomtype;
	delete [] z;
	delete [] coords;

	return true;
}

static bool skipStep(FILE *aPoscarFp, FILE *aEnergyFp)
{
	char linebuf[82];
	unsigned int j;

	// Read and ignore the title
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;

	// Read and ignore the scale factor
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;

	// Read and ignore the unit cell base vectors
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;

	// Read the array of atoms numbers
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;
	long int res;
	char *endptr;
	const char *next = linebuf;
	unsigned int natoms = 0;
	for(;;)
	{
		res = strtol(next, &endptr, 10);
		if(endptr == next) break;

		natoms += res;

		next = endptr;
	}

	// Skip the "Direct" line
	if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;

	// Read the atom records
	for(j=0; j < natoms; ++j)
	{
		if(fgets(linebuf, sizeof linebuf, aPoscarFp) == NULL) return false;
	}

	// Read the energy
	if(aEnergyFp)
	{
		if(fgets(linebuf, sizeof linebuf, aEnergyFp) == NULL) return false;
	}

	return true;
}


void readPoscarAndEnergies(const char* aFilenamePoscar,
							    const char* aFilenameEnergies,
							    bool aEnergyIsPerAtom,
							    unsigned int aStartStep,
							    unsigned int aEndStep,
							    std::vector<unsigned int>& aAtomZ,
							    CrystalFp& aCfp)
{
	if(!aFilenamePoscar)
	{
		throw CrystalFpFatal("No files given");
	}

	// Read POSCAR
	FILE* fpp = fopen(aFilenamePoscar, "r");
	if(!fpp)
	{
		fprintf(stderr, "Cannot open POSCAR file %s\n", aFilenamePoscar);
		throw CrystalFpFatal();
	}

	FILE *fpe = 0;
	if(aFilenameEnergies)
	{
		fpe = fopen(aFilenameEnergies, "r");
		if(!fpe)
		{
			fprintf(stderr, "Cannot open energy file %s\n", aFilenameEnergies);
			fclose(fpp);
			throw CrystalFpFatal();
		}
	}

	// Skip the inital steps
	unsigned int n = 1;
	for(; n < aStartStep; ++n) if(!skipStep(fpp, fpe)) {fclose(fpp); if(fpe) fclose(fpe); return;}

	// Read requested steps
	while((aEndStep < 1 || n <= aEndStep) && loadData(n, fpp, fpe, aEnergyIsPerAtom, aAtomZ, aCfp)) ++n;

	// Finish the batch loading of structure
	aCfp.addStructureBatchFinish();

	// Close and return
	fclose(fpp);
	if(fpe) fclose(fpe);
}
