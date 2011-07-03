
#include <map>
#include <iostream>
#include <fstream>
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


void StructureList::serialize(std::ofstream& aStream) const
{
	// Write number of structures
	unsigned int ns = mStructuresDirectory.size();
	aStream.write((char *)&ns, sizeof(unsigned int));

	// Write the included (or selected) structures
	unsigned int idx;
	for(idx=0; idx < ns; ++idx)
	{
		mStructuresDirectory[idx]->second.serialize(aStream);
	}
}


void StructureList::unserialize(std::ifstream& aStream, bool aAppend, int aStepOffset)
{
	// Read the number of structures
	unsigned int ns;
	aStream.read((char *)&ns, sizeof(unsigned int));

	// Recreate the structures list (empting it if no append required)
	if(!aAppend) mStructures.clear();

	unsigned int idx;
	for(idx=0; idx < ns; ++idx)
	{
		Structure s;
		s.unserialize(aStream);

		if(aAppend) s.mStepId += aStepOffset;
		mStructures.insert(std::pair<int,Structure>(s.mStepId, s));
	}
	mStructuresDirectory.clear();
	std::map<int,Structure>::iterator im;
	for(im=mStructures.begin(); im != mStructures.end(); ++im) mStructuresDirectory.push_back(im);
}

