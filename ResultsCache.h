
#ifndef RESULTS_CACHE_H
#define RESULTS_CACHE_H

#include <map>
#include <vector>
#include <string>

class ResultsCache
{
public:
	ResultsCache() {}

	~ResultsCache() {mCache.clear();}

	void clear(void) {mCache.clear();}

	void addResult(unsigned int aCategory, unsigned int aMethodIdx, unsigned int aLength, int aId1, int aId2, const float* aOut);

	bool retrieveResult(unsigned int aCategory, unsigned int aMethodIdx, unsigned int aLength, int aId1, int aId2, float* aOut);

private:
	char aKeyTmp[256];
	std::string aKey;
	std::map< std::string, std::vector<float> > mCache;

};

#endif

