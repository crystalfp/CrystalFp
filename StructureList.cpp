
#include <map>
#include <iostream>
#include "Structure.h"
#include "StructureList.h"
#include "CrystalFpExceptions.h"

using namespace cfp_internal;

StructureList::StructureList()
{
}


StructureList::~StructureList()
{
	clear();
}


void StructureList::addStructure(Structure& aElement)
{
	mStructures.insert(std::pair<int,Structure>(aElement.mStepId, aElement));
}


void StructureList::clear(void)
{
	mStructures.clear();
	mStructuresDirectory.clear();
}


void StructureList::filterOnEnergy(float aEnergyThreshold)
{
	std::map<int,Structure>::iterator im;

	mStructuresDirectory.clear();

	for(im=mStructures.begin(); im != mStructures.end(); ++im)
	{
		if(!im->second.mHasEnergy || im->second.mEnergyPerAtom <= aEnergyThreshold) mStructuresDirectory.push_back(im);
	}
}


void StructureList::selectAll(void)
{
	std::map<int,Structure>::iterator im;

	mStructuresDirectory.clear();

	for(im=mStructures.begin(); im != mStructures.end(); ++im) mStructuresDirectory.push_back(im);
}


void StructureList::selectNone(void)
{
	mStructuresDirectory.clear();
}


Structure& StructureList::getStructureBySequence(int aSequence)
{
	std::map<int,Structure>::iterator im;

	if((im = mStructures.find(aSequence)) == mStructures.end()) throw cfp::CrystalFpFatal("Sequence not found");
	return im->second;
}


Structure& StructureList::getStructureByIndex(unsigned int aIdx)
{
	if(aIdx >= mStructuresDirectory.size()) throw cfp::CrystalFpFatal("Index out of range");
	return mStructuresDirectory[aIdx]->second;
}


void StructureList::selectAllIncluded(const std::vector<bool>& aIncluded)
{
	std::vector<std::map<int, Structure>::iterator> temp_structures_directory;
	for(unsigned int idx=0; idx < mStructuresDirectory.size(); ++idx)
	{
		if(aIncluded[idx]) temp_structures_directory.push_back(mStructuresDirectory[idx]);
	}
	mStructuresDirectory = temp_structures_directory;
}
