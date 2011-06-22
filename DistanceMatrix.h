
#ifndef DISTANCEMATRIX_H
#define DISTANCEMATRIX_H

#include <vector>
#include <algorithm>
#include "CrystalFpExceptions.h"

class DistanceMatrix
{
public:
	DistanceMatrix(unsigned int aNumElements=0)
	{
		mNumElements = aNumElements;
		if(aNumElements == 0) return;
		mDistances.reserve((aNumElements*aNumElements-aNumElements)/2);
		mDistances.resize((aNumElements*aNumElements-aNumElements)/2);
	}

	~DistanceMatrix()
	{
		mDistances.clear();
	}
	
	// Copy constructor and assignment
	DistanceMatrix(const DistanceMatrix& aDistMat)
	{
		mNumElements = aDistMat.size();
		if(mNumElements == 0) return;
		mDistances.reserve((mNumElements*mNumElements-mNumElements)/2);
		mDistances = aDistMat.getVector();
	}

	DistanceMatrix& operator=(const DistanceMatrix& aDistMat)
	{
		// Make sure not same object
		if(this != &aDistMat)
		{
			mNumElements = aDistMat.size();
			if(mNumElements != 0)
			{
				mDistances.reserve((mNumElements*mNumElements-mNumElements)/2);
				mDistances = aDistMat.getVector();
			}
		}
		
		// Return ref for multiple assignment
		return *this;
	}

	void setVector(const std::vector<float>& aDistVect, unsigned int aNumElements)
	{
		mNumElements = aNumElements;
		if(mNumElements == 0) return;
		mDistances = aDistVect;
	}
	
	void resize(unsigned int aNumElements)
	{
		mNumElements = aNumElements;
		if(aNumElements == 0) return;
		mDistances.reserve((aNumElements*aNumElements-aNumElements)/2);
		mDistances.resize((aNumElements*aNumElements-aNumElements)/2);
	}
	
	const std::vector<float>& getVector(void) const
	{
		return mDistances;
	}

	unsigned int size(void) const
	{
		return mNumElements;
	}

	bool empty(void) const
	{
		return mNumElements == 0;
	}
	
	float min() const
	{
		return *std::min_element(mDistances.begin(), mDistances.end());
	}
	
	float max() const
	{
		return *std::max_element(mDistances.begin(), mDistances.end());
	}

	float& operator() (unsigned int aRow, unsigned int aCol)
	{
		if(aRow >= mNumElements) throw cfp::CrystalFpFatal("Invalid aRow in DistanceMatrix");
		if(aCol >= mNumElements) throw cfp::CrystalFpFatal("Invalid aCol in DistanceMatrix");

		if(aRow == aCol) {mDummy = 0.0F; return mDummy;}

		unsigned int idx;
		if(aRow < aCol)
		{
			idx = aRow*(2*mNumElements-aRow-1)/2 + aCol - aRow - 1;
		}
		else
		{
			idx = aCol*(2*mNumElements-aCol-1)/2 + aRow - aCol - 1;
		}

		return mDistances[idx];
	}

	float operator() (unsigned int aRow, unsigned int aCol) const
	{
		if(aRow >= mNumElements) throw cfp::CrystalFpFatal("Invalid aRow in DistanceMatrix");
		if(aCol >= mNumElements) throw cfp::CrystalFpFatal("Invalid aCol in DistanceMatrix");

		if(aRow == aCol) return 0.0F;

		unsigned int idx;
		if(aRow < aCol)
		{
			idx = aRow*(2*mNumElements-aRow-1)/2 + aCol - aRow - 1;
		}
		else
		{
			idx = aCol*(2*mNumElements-aCol-1)/2 + aRow - aCol - 1;
		}

		return mDistances[idx];
	}

	void resizeToIncluded(const std::vector<bool>& aIncluded)
	{
		unsigned int i, j;

		// Obtain the new size
		unsigned int new_size = 0;
		for(i=0; i < aIncluded.size(); ++i) if(aIncluded[i]) ++new_size;

		// Check
		if(new_size == mNumElements) return;
		if(new_size > mNumElements) throw cfp::CrystalFpFatal("Resize attempted to a bigger distance matrix");

		// Map old row/col to new ones
		std::vector<unsigned int> map;
		map.resize(mNumElements);
		for(i=j=0; i < aIncluded.size(); ++i) if(aIncluded[i]) map[i] = j++;

		// Temporary distance matrix
		std::vector<float> new_distances;
		new_distances.resize((new_size*new_size-new_size)/2);

		// Remove unselected row/columns
		for(unsigned int row=0; row < mNumElements-1; ++row)
		{
			if(!aIncluded[row]) continue;

			for(unsigned int col=row+1; col < mNumElements; ++col)
			{
				if(!aIncluded[col]) continue;

				unsigned int idx = row*(2*mNumElements-row-1)/2 + col - row - 1;

				unsigned int nrow = map[row];
				unsigned int ncol = map[col];
			
				unsigned int nidx = nrow*(2*new_size-nrow-1)/2 + ncol - nrow - 1;

				new_distances[nidx] = mDistances[idx];
			}
		}

		// Update
		mDistances   = new_distances;
		mNumElements = new_size;
	}

	void clear(void)
	{
		mDistances.clear();
		mNumElements = 0;
	}

private:
	std::vector<float>	mDistances;
	unsigned int		mNumElements;
	float				mDummy;
};

#endif
