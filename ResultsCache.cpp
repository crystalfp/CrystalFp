
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <cstdio>
#include <cstring>
#include <iostream>
#include "ResultsCache.h"

bool ResultsCache::retrieveResult(unsigned int aCategory, unsigned int aMethodIdx, unsigned int aLength, int aId1, int aId2, float* aOut)
{
	// Compute the key
	sprintf(aKeyTmp, "%u|%u|%u|%d|%d", aCategory, aMethodIdx, aLength, aId1, aId2);
	aKey = aKeyTmp;

	// If it is in cache, copy the value to the output array
	std::map< std::string, std::vector<float> >::const_iterator im = mCache.find(aKey);
	if(im != mCache.end())
	{
		memcpy(aOut, &(im->second[0]), aLength*sizeof(float));
		return true;
	}
	return false;
}

void ResultsCache::addResult(unsigned int aCategory, unsigned int aMethodIdx, unsigned int aLength, int aId1, int aId2, const float* aOut)
{
	// Compute the key
	sprintf(aKeyTmp, "%u|%u|%u|%d|%d", aCategory, aMethodIdx, aLength, aId1, aId2);
	aKey = aKeyTmp;

	// Cache the value
	std::vector<float> x;
	x.assign(aOut, aOut+aLength);
	mCache[aKey] = x;
}


