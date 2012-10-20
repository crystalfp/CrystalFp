
//#define TRACE_GROUPING

#ifdef TRACE_GROUPING
#include <iostream>
#include <iomanip>
#endif
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstring>
#include "GroupingMethod.h"

using namespace cfp_internal;

void PseudoSNNGrouping::computeConnectionMatrix(size_t aNumStructures, const DistanceMatrix& aDistances, unsigned int* aConnection) const
{
	// Build the connection matrix
	for(unsigned int r=0; r < aNumStructures; ++r)
	{
		for(unsigned int c=r+1; c < aNumStructures; ++c)
		{
			float dist = aDistances(r, c);
			aConnection[c*aNumStructures+r] = aConnection[r*aNumStructures+c] = (dist < mMaxDistanceForGrouping) ? 1 : 0;
		}
	}
	for(unsigned int k=0; k < aNumStructures; ++k) aConnection[k*(aNumStructures+1)] = 0;
}


void PseudoSNNGrouping::doDepthFirstVisit(unsigned int aIdx, bool *aAssigned, std::set<unsigned int>& aGroups, unsigned int* aConnection, size_t aNumStructures) const
{
	// Assign the node to the connected component
	aAssigned[aIdx] = true;
	aGroups.insert(aIdx);

	// For each node to which it is connected
	for(unsigned int j=0; j < aNumStructures; ++j)
	{
		if(aAssigned[j]) continue;
		if(aConnection[aIdx*aNumStructures+j] != 0) doDepthFirstVisit(j, aAssigned, aGroups, aConnection, aNumStructures);
	}
}

void PseudoSNNGrouping::doGrouping(size_t aNumStructures, const DistanceMatrix& aDistances, std::vector< std::set<unsigned int> >& aResult)
{
	// Compute the connection matrix
	unsigned int* connection = new unsigned int[aNumStructures*aNumStructures];
	computeConnectionMatrix(aNumStructures, aDistances, connection);

	// Compute the number of shared NN
	unsigned int idx1, idx2, j;
	for(idx1=0; idx1 < aNumStructures-1; ++idx1)
	{
		for(idx2=idx1+1; idx2 < aNumStructures; ++idx2)
		{
			// Do nothing if nodes not connected
			if(connection[idx1*aNumStructures+idx2] == 0) continue;

			// For all the connections to node1
			for(j=0; j < aNumStructures; ++j)
			{
				// if j is shared between node-1 and node-2
				if(j != idx1 && j != idx2 && connection[idx1*aNumStructures+j] != 0 && connection[idx2*aNumStructures+j] != 0)
				{
					++connection[idx1*aNumStructures+idx2];
					++connection[idx2*aNumStructures+idx1];
				}
			}
		}
	}

	// Now remove connections with only one connection except for standalone pairs
	for(idx1=0; idx1 < aNumStructures-1; ++idx1)
	{
		for(idx2=idx1+1; idx2 < aNumStructures; ++idx2)
		{
			// Do nothing if nodes have none or more than one SNN
			if(connection[idx1*aNumStructures+idx2] != 1) continue;

			// Count the number of connections from node1
			unsigned int nconn1 = 0;
			for(j=0; j < aNumStructures; ++j) if(connection[idx1*aNumStructures+j] != 0) ++nconn1;

			// Count the number of connections from node2
			unsigned int nconn2 = 0;
			for(j=0; j < aNumStructures; ++j) if(connection[idx2*aNumStructures+j] != 0) ++nconn2;

			// If it is a bridge, break it
			if(nconn1 > 1 && nconn2 > 1) connection[idx1*aNumStructures+idx2] = 0;
		}
	}

	// If requested remove nodes with less than K shared nearest neighbors
	if(mK > 0)
	{
		for(idx1=0; idx1 < aNumStructures-1; ++idx1)
		{
			for(idx2=idx1+1; idx2 < aNumStructures; ++idx2)
			{
				// Do nothing if nodes have none or more than one SNN
				if(connection[idx1*aNumStructures+idx2] == 0) continue;
				if(connection[idx1*aNumStructures+idx2] < (1+mK)) connection[idx1*aNumStructures+idx2] = 0;
			}
		}
	}

	// Make the matrix symmetrical again
	for(idx1=0; idx1 < aNumStructures-1; ++idx1)
	{
		for(idx2=idx1+1; idx2 < aNumStructures; ++idx2)
		{
			connection[idx2*aNumStructures+idx1] = connection[idx1*aNumStructures+idx2];
		}
	}

	// Array to keep track of nodes already assigned to a group
	bool *assigned = new bool[aNumStructures];
	for(j=0; j < aNumStructures; ++j) assigned[j] = false;

	// For each node do a DFS to extract the nodes
	for(idx1=0; idx1 < aNumStructures; ++idx1)
	{
		// Skip if already assigned
		if(assigned[idx1]) continue;

		// Start a new group
		std::set<unsigned int> s;

		// Find the connected component
		doDepthFirstVisit(idx1, assigned, s, connection, aNumStructures);

		// Insert the new group
		aResult.push_back(s);
	}

	// Release and return
	delete [] assigned;
	delete [] connection;
}



float HierarchicalGrouping::ClusterDistanceSingleLinkage(std::list<HierarchicalGrouping::Node>::iterator& aNi,
														 std::list<HierarchicalGrouping::Node>::iterator& aNj,
														 const DistanceMatrix& aDistances) const
{
	size_t leni = aNi->idx.size();
	size_t lenj = aNj->idx.size();
	if(leni == 1 && lenj == 1)
	{
		return aDistances(aNi->idx[0], aNj->idx[0]);
	}
	else
	{
		float dist = FLT_MAX;
		for(size_t i=0; i < leni; ++i)
		{
			for(size_t j=0; j < lenj; ++j)
			{
				float d = aDistances(aNi->idx[i], aNj->idx[j]);

				if(d < dist) dist = d;
			}
		}
		return dist;
	}
}

float HierarchicalGrouping::ClusterDistanceCompleteLinkage(std::list<HierarchicalGrouping::Node>::iterator& aNi,
														   std::list<HierarchicalGrouping::Node>::iterator& aNj,
														   const DistanceMatrix& aDistances) const
{
	size_t leni = aNi->idx.size();
	size_t lenj = aNj->idx.size();
	if(leni == 1 && lenj == 1)
	{
		return aDistances(aNi->idx[0], aNj->idx[0]);
	}
	else
	{
		float dist = -FLT_MAX;
		for(size_t i=0; i < leni; ++i)
		{
			for(size_t j=0; j < lenj; ++j)
			{
				float d = aDistances(aNi->idx[i], aNj->idx[j]);
				if(d > dist) dist = d;
			}
		}
		return dist;
	}
}


void HierarchicalGrouping::doGrouping(size_t aNumStructures, const DistanceMatrix& aDistances, std::vector< std::set<unsigned int> >& aResult)
{
	std::list<Node>::const_iterator n;

	// Initialize root (to point to all)
	std::list<Node> root;
	for(unsigned int i=0; i < aNumStructures; ++i) root.push_back(Node(i));
#ifdef TRACE_GROUPING
	float curr_min = 0.0F;
#endif

	// Iterate till the distance becomes greather than the given threshold
	for(; root.size() > 1;)
	{
#ifdef TRACE_GROUPING
		// Corresponding threshold value
		std::cerr << std::fixed << std::setprecision(4) << curr_min << "  ";

		// Debug: output initial grouping
		for(n=root.begin(); n != root.end(); ++n)
		{
			std::cerr << "(" ;
			std::vector<unsigned int>::const_iterator ii;
			for(ii=n->idx.begin(); ii != n->idx.end(); ++ii)
			{
				if(ii != n->idx.begin()) std::cerr << " " ;
				std::cerr << *ii;
			}
			std::cerr << ")" ;
		}
		std::cerr << std::endl;
#endif

		// Find minimum dist and minimum pairs
		std::list<Node>::iterator ni;
		std::list<Node>::iterator nj;
		float min_dist = FLT_MAX;
		std::list<Node>::iterator min_i;
		std::list<Node>::iterator min_j;
		for(ni=root.begin(); ni != root.end(); ++ni)
		{
			for(nj=ni, ++nj; nj != root.end(); ++nj)
			{
				float dist = (this->*ClusterDistance)(ni, nj, aDistances);
				if(dist < min_dist)
				{
					min_dist = dist;
					min_i    = ni;
					min_j    = nj;
				}
			}
		}

		// Exit if the threshold has been reached
		if(min_dist > mMaxDistanceForGrouping) break;

#ifdef TRACE_GROUPING
		// Update min distance
		curr_min = min_dist;
#endif

		// Update the group list. Merge node j at the end of node i
		min_i->Merge(min_j);

		// Remove merged node
		root.erase(min_j);
	}

#ifdef TRACE_GROUPING
	// If the last step is a single cluster, display it
	if(root.size() == 1)
	{
		// Corresponding threshold value
		std::cerr << std::fixed << std::setprecision(4) << curr_min << "  ";

		// Debug: output initial grouping
		for(n=root.begin(); n != root.end(); ++n)
		{
			std::cerr << "(" ;
			std::vector<unsigned int>::const_iterator ii;
			for(ii=n->idx.begin(); ii != n->idx.end(); ++ii)
			{
				if(ii != n->idx.begin()) std::cerr << " " ;
				std::cerr << *ii;
			}
			std::cerr << ")" ;
		}
		std::cerr << std::endl;
	}
#endif

	// Load the configuration in the final structure
	for(n=root.begin(); n != root.end(); ++n)
	{
		std::vector<unsigned int>::const_iterator ii;
		std::set<unsigned int> s;
		for(ii=n->idx.begin(); ii != n->idx.end(); ++ii) s.insert(*ii) ;
		aResult.push_back(s);
	}
}


