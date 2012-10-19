
#ifndef GROUPINGMETHOD_H
#define GROUPINGMETHOD_H

#include <vector>
#include <set>
#include <string>
#include <list>
#include "MethodsList.h"
#include "Structure.h"
#include "DistanceMatrix.h"

namespace cfp_internal
{

class GroupingMethod
{
public:
	GroupingMethod()
	{
		mK = 0;
		mMaxDistanceForGrouping = 0.01F;
	}
	std::string getName(void) const { return mName; }
	virtual ~GroupingMethod() {}
	void setKvalue(unsigned int aK) { mK = aK;}
	void setMaxDistanceForGrouping(float aMaxDistanceForGrouping) { mMaxDistanceForGrouping = aMaxDistanceForGrouping;}
	virtual void doGrouping(size_t aNumStructures, const DistanceMatrix& aDistances, std::vector< std::set<unsigned int> >& aResult) =0;
	virtual bool needsK(void) const {return false;}

protected:
	std::string		mName;
	unsigned int	mK;
	float			mMaxDistanceForGrouping;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class PseudoSNNGrouping : public GroupingMethod
{
public:
	PseudoSNNGrouping() {mName = "Pseudo SNN";}
	virtual void doGrouping(size_t aNumStructures, const DistanceMatrix& aDistances, std::vector< std::set<unsigned int> >& aResult);
	void computeConnectionMatrix(size_t aNumStructures, const DistanceMatrix& aDistances, unsigned int* aConnection) const;
	void doDepthFirstVisit(unsigned int aIdx, bool *aAssigned, std::set<unsigned int>& aGroups, unsigned int* aConnection, size_t aNumStructures) const;
	virtual bool needsK(void) const {return true;}
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class HierarchicalGrouping : public GroupingMethod
{
public:
	virtual void doGrouping(size_t aNumStructures, const DistanceMatrix& aDistances, std::vector< std::set<unsigned int> >& aResult);

protected:
	struct Node
	{
		Node() {}
		Node(unsigned int i) { idx.push_back(i); }

		std::vector<unsigned int> idx;

		void Merge(const Node& n)
		{
			std::vector<unsigned int>::const_iterator ii;
			for(ii=n.idx.begin(); ii != n.idx.end(); ++ii) idx.push_back(*ii);
		}
		void Merge(const std::list<Node>::iterator& n)
		{
			std::vector<unsigned int>::const_iterator ii;
			for(ii=n->idx.begin(); ii != n->idx.end(); ++ii) idx.push_back(*ii);
		}
	};

public:
	float ClusterDistanceSingleLinkage(std::list<HierarchicalGrouping::Node>::iterator& aNi,
									   std::list<HierarchicalGrouping::Node>::iterator& aNj,
									   const DistanceMatrix& aDistances) const;
	float ClusterDistanceCompleteLinkage(std::list<HierarchicalGrouping::Node>::iterator& aNi,
										 std::list<HierarchicalGrouping::Node>::iterator& aNj,
										 const DistanceMatrix& aDistances) const;

protected:
	float (HierarchicalGrouping::*ClusterDistance)(std::list<HierarchicalGrouping::Node>::iterator& aNi,
												   std::list<HierarchicalGrouping::Node>::iterator& aNj,
												   const DistanceMatrix& aDistances) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class HierarchicalSingleLinkageGrouping : public HierarchicalGrouping
{
public:
	HierarchicalSingleLinkageGrouping()
	{
		mName = "Hierarchical grouping (single linkage)";
		ClusterDistance = &HierarchicalGrouping::ClusterDistanceSingleLinkage;
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class HierarchicalCompleteLinkageGrouping : public HierarchicalGrouping
{
public:
	HierarchicalCompleteLinkageGrouping()
	{
		mName = "Hierarchical grouping (complete linkage)";
		ClusterDistance = &HierarchicalGrouping::ClusterDistanceCompleteLinkage;
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class GroupingMethodsList : public MethodsList<GroupingMethod>
{
public:
	GroupingMethodsList()
	{
		mMethodsList.push_back(new PseudoSNNGrouping);
		mMethodsList.push_back(new HierarchicalSingleLinkageGrouping);
		mMethodsList.push_back(new HierarchicalCompleteLinkageGrouping);
	}
	~GroupingMethodsList() {clear();}
};

}
#endif

