
#ifndef READPOSCAR_H
#define READPOSCAR_H

#include <vector>
using namespace cfp;

extern void readPoscarAndEnergies(const char* aFilenamePoscar,
									   const char* aFilenameEnergies,
									   bool aEnergyIsPerAtom,
									   unsigned int aStartStep,
									   unsigned int aEndStep,
									   std::vector<unsigned int>& aAtomZ,
									   CrystalFp& aCfp);

#endif

