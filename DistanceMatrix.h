
#ifndef DISTANCEMATRIX_H
#define DISTANCEMATRIX_H

#include <vector>
#include <algorithm>
#include <fstream>
#include "CrystalFpExceptions.h"

class DistanceMatrix
{
public:
	DistanceMatrix(size_t aNumElements=0) : mNumElements(aNumElements), mDummy(0.0F)
	{
		if(aNumElements == 0) return;
		mDistances.reserve((aNumElements*aNumElements-aNumElements)/2);
		mDistances.resize((aNumElements*aNumElements-aNumElements)/2);
	}

	~DistanceMatrix()
	{
		mDistances.clear();
	}
	
	// Copy constructor and assignment
	DistanceMatrix(const DistanceMatrix& aDistMat) : mDummy(0.0F)
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
			mDummy = 0.0F;
		}
		
		// Return ref for multiple assignment
		return *this;
	}

	void setVector(const std::vector<float>& aDistVect, size_t aNumElements)
	{
		mNumElements = aNumElements;
		if(mNumElements == 0) return;
		mDistances = aDistVect;
	}
	
	void resize(size_t aNumElements)
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

	size_t size(void) const
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

	float& operator() (size_t aRow, size_t aCol)
	{
		if(aRow >= mNumElements) throw cfp::CrystalFpFatal("Invalid aRow in DistanceMatrix");
		if(aCol >= mNumElements) throw cfp::CrystalFpFatal("Invalid aCol in DistanceMatrix");

		if(aRow == aCol) {mDummy = 0.0F; return mDummy;}

		size_t idx;
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

	float operator() (size_t aRow, size_t aCol) const
	{
		if(aRow >= mNumElements) throw cfp::CrystalFpFatal("Invalid aRow in DistanceMatrix");
		if(aCol >= mNumElements) throw cfp::CrystalFpFatal("Invalid aCol in DistanceMatrix");

		if(aRow == aCol) return 0.0F;

		size_t idx;
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
		size_t i;
		unsigned int j;

		// Obtain the new size
		size_t new_size = 0;
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
		for(size_t row=0; row < mNumElements-1; ++row)
		{
			if(!aIncluded[row]) continue;

			for(size_t col=row+1; col < mNumElements; ++col)
			{
				if(!aIncluded[col]) continue;

				size_t idx = row*(2*mNumElements-row-1)/2 + col - row - 1;

				unsigned int nrow = map[row];
				unsigned int ncol = map[col];
			
				size_t nidx = nrow*(2*new_size-nrow-1)/2 + ncol - nrow - 1;

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

	
	void serialize(std::ofstream& aStream) const
	{
		unsigned int x = static_cast<unsigned int>(mDistances.size());
		aStream.write((char *)&x, sizeof(unsigned int));
		if(x) aStream.write((char *)&mDistances[0], sizeof(float)*x);
		x = static_cast<unsigned int>(mNumElements);
		aStream.write((char *)&x, sizeof(unsigned int));
	}


	void unserialize(std::ifstream& aStream)
	{
		unsigned int x;
		aStream.read((char *)&x, sizeof(unsigned int));
		mDistances.resize(static_cast<size_t>(x));
		if(x) aStream.read((char *)&mDistances[0], sizeof(float)*x);
		aStream.read((char *)&x, sizeof(unsigned int));
		mNumElements = x;
	}

private:
	std::vector<float>	mDistances;
	size_t				mNumElements;
	float				mDummy;
};

#endif
