
#ifndef DISTANCEMETHOD_H
#define DISTANCEMETHOD_H

#include <string>
#include <vector>
#include "MethodsList.h"
#include "Structure.h"
#include "StructureList.h"

namespace cfp_internal
{

class DistanceMethod
{
public:
	DistanceMethod()
	{
	};
	std::string getName(void) const { return mName; }
	virtual ~DistanceMethod() {};
	virtual float computeDistance(const Structure& aStructure1, const Structure& aStructure2) =0;

protected:
	std::string mName;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class CosineDistance : public DistanceMethod
{
public:
	CosineDistance() {mName = "Cosine distance";}
	virtual float computeDistance(const Structure& aStructure1, const Structure& aStructure2);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class EuclideanDistance : public DistanceMethod
{
public:
	EuclideanDistance() {mName = "Euclidean distance";}
	virtual float computeDistance(const Structure& aStructure1, const Structure& aStructure2);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MinkowskiDistance : public DistanceMethod
{
public:
	MinkowskiDistance() {mName = "Minkowski (p=1/3) distance";}
	virtual float computeDistance(const Structure& aStructure1, const Structure& aStructure2);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class DistanceMethodsList : public MethodsList<DistanceMethod>
{
public:
	DistanceMethodsList()
	{
		mMethodsList.push_back(new CosineDistance);
		mMethodsList.push_back(new EuclideanDistance);
		mMethodsList.push_back(new MinkowskiDistance);
	}
	~DistanceMethodsList() {clear();}
};

}
#endif

