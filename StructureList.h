
#ifndef CRYSTALFP_STRUCTURELIST_H
#define CRYSTALFP_STRUCTURELIST_H

#ifdef _MSC_VER
#pragma warning(disable:4786)
#pragma warning(disable:4290)
#endif

#include <map>
#include <vector>
#include "Structure.h"

namespace cfp_internal
{

class StructureList
{
public:
	StructureList();
	~StructureList();

	void addStructure(Structure& aElement);

	void clear(void);

	typedef std::vector<std::map<int, Structure>::iterator>::iterator iterator;
	typedef std::vector<std::map<int, Structure>::iterator>::const_iterator const_iterator;

	iterator begin(void) {return mStructuresDirectory.begin();}
	iterator end(void)   {return mStructuresDirectory.end();}
	const_iterator begin(void) const {return mStructuresDirectory.begin();}
	const_iterator end(void) const   {return mStructuresDirectory.end();}

	void filterOnEnergy(float aEnergyThreshold);
	void selectAll(void);
	void selectNone(void);
	void selectAllIncluded(const std::vector<bool>& aIncluded);

	unsigned int getTotalCount(void) const {return mStructures.size();}
	unsigned int getSelectedCount(void) const {return mStructuresDirectory.size();}
	unsigned int size(void) const {return mStructuresDirectory.size();}

	Structure& getStructureBySequence(int aSequence);
	Structure& getStructureByIndex(unsigned int aIdx);

	void serialize(std::ofstream& aStream) const;
	void unserialize(std::ifstream& aStream, bool aAppend=false, int aStepOffset=10000);

private:
	std::map<int, Structure> mStructures;
	std::vector<std::map<int, Structure>::iterator> mStructuresDirectory;
};

}
#endif

