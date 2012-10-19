
#include <iostream>
#include <cfloat>
#include <cstdlib>

#include "CrystalFpAnalysis.h"
#include "AnalysisMethod.h"
#include "CrystalFpExceptions.h"
#include "ResultsCache.h"

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif
static const unsigned int NO_METHOD_SELECTED = UINT_MAX;

using namespace cfp_internal;
using namespace cfp;

/// The private data for the CrystalFpAnalysis class.
///
struct CrystalFpAnalysis::CrystalFpAnalysisImpl
{
	AnalysisMethodsSimpleList	mSimpleMethodsList;
	AnalysisMethodsHistList		mHistMethodsList;
	AnalysisMethodsSpecialList	mSpecialMethodsList;
	unsigned int				mAnalysisTypeX;
	unsigned int				mAnalysisTypeY;
	unsigned int				mAnalysisCategory;
	AnalysisMethod*				mAnalysisMethodX;
	AnalysisMethod*				mAnalysisMethodY;
	size_t						mHistNumBins;
	unsigned int				mStructureIdx;
	ResultsCache				mCache;
	bool						mCacheResults;

};

CrystalFpAnalysis::CrystalFpAnalysis(const CrystalFp *aCfp) : mPimpl(new CrystalFpAnalysisImpl())
{
	mCrystalFp = aCfp;

	mPimpl->mAnalysisTypeX    = NO_METHOD_SELECTED;
	mPimpl->mAnalysisTypeY    = NO_METHOD_SELECTED;
	mPimpl->mAnalysisCategory = NO_METHOD_SELECTED;

	mPimpl->mHistNumBins      = 100;
	mPimpl->mStructureIdx     = 0;

	mPimpl->mAnalysisMethodX  = 0;
	mPimpl->mAnalysisMethodY  = 0;
	mPimpl->mCacheResults     = false;
}

CrystalFpAnalysis::~CrystalFpAnalysis()
{
	mPimpl->mCache.clear();
	delete mPimpl;
}

void CrystalFpAnalysis::resetCache(const CrystalFp *aCfp)
{
	if(aCfp) mCrystalFp = aCfp;
	
	mPimpl->mCache.clear();
}

void CrystalFpAnalysis::enableCaching(void)
{
	mPimpl->mCacheResults = true;
}


const std::vector<std::string> CrystalFpAnalysis::getAnalysisMethodsNames(unsigned int aCategory) const
{
	std::vector<std::string> names;

	switch(aCategory)
	{
	case CATEGORY_SIMPLE:
		for(AnalysisMethodsSimpleList::const_iterator ifml=mPimpl->mSimpleMethodsList.begin(); ifml != mPimpl->mSimpleMethodsList.end(); ++ifml)
		{
			names.push_back((*ifml)->getName());
		}
		break;

	case CATEGORY_HIST:
		for(AnalysisMethodsHistList::const_iterator ifml=mPimpl->mHistMethodsList.begin(); ifml != mPimpl->mHistMethodsList.end(); ++ifml)
		{
			names.push_back((*ifml)->getNameHist());
		}
		break;

	case CATEGORY_SPECIAL:
		for(AnalysisMethodsSpecialList::const_iterator ifml=mPimpl->mSpecialMethodsList.begin(); ifml != mPimpl->mSpecialMethodsList.end(); ++ifml)
		{
			names.push_back((*ifml)->getName());
		}
		break;

	case CATEGORY_ALL:
		names.push_back("All per structure analysis");
		break;
	}
	return names;
}

bool CrystalFpAnalysis::setAnalysisMethod(unsigned int aAnalysisTypeX, unsigned int aAnalysisTypeY)
{
	unsigned int n1 = static_cast<unsigned int>(mPimpl->mSimpleMethodsList.size());
	if(aAnalysisTypeX < n1)
	{
		if(aAnalysisTypeY >= n1) return false;
		mPimpl->mAnalysisTypeX = aAnalysisTypeX;
		mPimpl->mAnalysisTypeY = aAnalysisTypeY;
		mPimpl->mAnalysisCategory = CATEGORY_SIMPLE;
		mPimpl->mAnalysisMethodX = mPimpl->mSimpleMethodsList.getMethod(aAnalysisTypeX);
		mPimpl->mAnalysisMethodY = mPimpl->mSimpleMethodsList.getMethod(aAnalysisTypeY);

		if(!mPimpl->mAnalysisMethodX->isValid(mCrystalFp)) return false;
		if(!mPimpl->mAnalysisMethodY->isValid(mCrystalFp)) return false;
		if(mPimpl->mAnalysisMethodX->numElements(mCrystalFp) != mPimpl->mAnalysisMethodY->numElements(mCrystalFp)) return false;
		return true;
	}

	unsigned int n2 = static_cast<unsigned int>(mPimpl->mHistMethodsList.size());
	if(aAnalysisTypeX < (n1+n2))
	{
		mPimpl->mAnalysisTypeX = aAnalysisTypeX-n1;
		mPimpl->mAnalysisCategory = CATEGORY_HIST;
		mPimpl->mAnalysisMethodX = mPimpl->mHistMethodsList.getMethod(aAnalysisTypeX-n1);

		if(!mPimpl->mAnalysisMethodX->isValid(mCrystalFp))  return false;
		return true;
	}

	unsigned int n3 = static_cast<unsigned int>(mPimpl->mSpecialMethodsList.size());
	if(aAnalysisTypeX < (n1+n2+n3))
	{
		mPimpl->mAnalysisTypeX = aAnalysisTypeX-n1-n2;
		mPimpl->mAnalysisCategory = CATEGORY_SPECIAL;
		mPimpl->mAnalysisMethodX = mPimpl->mSpecialMethodsList.getMethod(aAnalysisTypeX-n1-n2);

		if(!mPimpl->mAnalysisMethodX->isValid(mCrystalFp))  return false;
		return true;
	}

	if(aAnalysisTypeX == (n1+n2+n3))
	{
		mPimpl->mAnalysisTypeX = 0;
		mPimpl->mAnalysisCategory = CATEGORY_ALL;
		return true;
	}
	
	 return false;
}


void CrystalFpAnalysis::setNamedParam(const std::string& aName, const std::string& aValue)
{
	if(!strncasecmp(aName.c_str(), "bin", 3))
	{
		mPimpl->mHistNumBins = (unsigned int)atoi(aValue.c_str());
	}
	else if(!strncasecmp(aName.c_str(), "part", 4))
	{
		if(mPimpl->mAnalysisCategory == CATEGORY_SPECIAL) mPimpl->mAnalysisMethodX->setFingerprintPart((unsigned int)atoi(aValue.c_str()));
	}
	else if(!strncasecmp(aName.c_str(), "idx", 4))
	{
		if(mPimpl->mAnalysisCategory == CATEGORY_SPECIAL)
		{
			mPimpl->mStructureIdx = (unsigned int)atoi(aValue.c_str());
			mPimpl->mAnalysisMethodX->setStructureIdx(mPimpl->mStructureIdx);
		}
	}
	else
	{
		throw CrystalFpFatal("Invalid named param for CrystalFpAnalysis");
	}
}

void CrystalFpAnalysis::setNamedParam(const std::string& aName, unsigned int aValue)
{
	if(!strncasecmp(aName.c_str(), "bin", 3))
	{
		mPimpl->mHistNumBins = aValue;
	}
	else if(!strncasecmp(aName.c_str(), "part", 4))
	{
		if(mPimpl->mAnalysisCategory == CATEGORY_SPECIAL) mPimpl->mAnalysisMethodX->setFingerprintPart(aValue);
	}
	else if(!strncasecmp(aName.c_str(), "idx", 4))
	{
		if(mPimpl->mAnalysisCategory == CATEGORY_SPECIAL)
		{
			mPimpl->mStructureIdx = aValue;
			mPimpl->mAnalysisMethodX->setStructureIdx(aValue);
		}
	}
}


std::vector<size_t> CrystalFpAnalysis::numValues(void)
{
	std::vector<size_t> nums;

	if(mPimpl->mAnalysisCategory == CATEGORY_SIMPLE)
	{
		nums.push_back(1);
		nums.push_back(mPimpl->mAnalysisMethodX->numElements(mCrystalFp));
		nums.push_back(mPimpl->mAnalysisMethodY->numElements(mCrystalFp));
	}
	else if(mPimpl->mAnalysisCategory == CATEGORY_HIST)
	{
		nums.push_back(1);
		nums.push_back(mPimpl->mHistNumBins);
		nums.push_back(mPimpl->mHistNumBins);
	}
	else if(mPimpl->mAnalysisCategory == CATEGORY_SPECIAL)
	{
		nums.push_back(1);
		for(unsigned int i=0; i < mPimpl->mAnalysisMethodX->numArrays(); ++i)
		{
			nums.push_back(mPimpl->mAnalysisMethodX->numElements(mCrystalFp, i));
		}
	}
	else if(mPimpl->mAnalysisCategory == CATEGORY_ALL)
	{
		nums.push_back(1);
		for(AnalysisMethodsSimpleList::const_iterator ifml=mPimpl->mSimpleMethodsList.begin(); ifml != mPimpl->mSimpleMethodsList.end(); ++ifml)
		{
			if((*ifml)->numElements(mCrystalFp) != mCrystalFp->getNumActiveStructures()) continue;
			nums.push_back(mCrystalFp->getNumActiveStructures());
		}
	}
	else
	{
		throw CrystalFpFatal("Analysis method not set");
	}
	
	return nums;
}

const std::vector<std::string> CrystalFpAnalysis::getLabels(void) const
{
	std::vector<std::string> labels;
	std::vector<std::string> tmp;

	switch(mPimpl->mAnalysisCategory)
	{
	case CATEGORY_SIMPLE:
		labels = mPimpl->mAnalysisMethodX->getLabels();
		tmp    = mPimpl->mAnalysisMethodY->getLabels();
		labels.insert(labels.end(), tmp.begin(), tmp.end());
		break;

	case CATEGORY_HIST:
		labels = mPimpl->mAnalysisMethodX->getLabels();
		labels.push_back("Frequency");
		break;

	case CATEGORY_SPECIAL:
		labels = mPimpl->mAnalysisMethodX->getLabels(mCrystalFp);
		break;

	case CATEGORY_ALL:
		for(AnalysisMethodsSimpleList::const_iterator ifml=mPimpl->mSimpleMethodsList.begin(); ifml != mPimpl->mSimpleMethodsList.end(); ++ifml)
		{
			if((*ifml)->numElements(mCrystalFp) != mCrystalFp->getNumActiveStructures()) continue;
			tmp = (*ifml)->getLabels();
			labels.insert(labels.end(), tmp.begin(), tmp.end());
		}
		break;
	}

	return labels;
}

void CrystalFpAnalysis::getValues(float **aValues) const
{
	unsigned int i;

	switch(mPimpl->mAnalysisCategory)
	{
	case CATEGORY_SIMPLE:
		if(mPimpl->mCacheResults)
		{
			size_t sz = mPimpl->mAnalysisMethodX->numElements(mCrystalFp);

			if(!mPimpl->mCache.retrieveResult(CATEGORY_SIMPLE, mPimpl->mAnalysisTypeX, sz, 0, 0, aValues[0]))
			{
				mPimpl->mAnalysisMethodX->getValues(mCrystalFp, aValues[0]);

				mPimpl->mCache.addResult(CATEGORY_SIMPLE, mPimpl->mAnalysisTypeX, sz, 0, 0, aValues[0]);
			}
			if(!mPimpl->mCache.retrieveResult(CATEGORY_SIMPLE, mPimpl->mAnalysisTypeY, sz, 0, 0, aValues[1]))
			{
				mPimpl->mAnalysisMethodY->getValues(mCrystalFp, aValues[1]);

				mPimpl->mCache.addResult(CATEGORY_SIMPLE, mPimpl->mAnalysisTypeY, sz, 0, 0, aValues[1]);
			}
		}
		else
		{
			mPimpl->mAnalysisMethodX->getValues(mCrystalFp, aValues[0]);
			mPimpl->mAnalysisMethodY->getValues(mCrystalFp, aValues[1]);
		}
		break;

	case CATEGORY_HIST:
		{
		size_t sz = mPimpl->mAnalysisMethodX->numElements(mCrystalFp);
		float* val = new float[sz];
		if(mPimpl->mCacheResults)
		{
			if(!mPimpl->mCache.retrieveResult(CATEGORY_HIST, mPimpl->mAnalysisTypeX, sz, 0, 0, val))
			{
				mPimpl->mAnalysisMethodX->getValues(mCrystalFp, val);

				mPimpl->mCache.addResult(CATEGORY_HIST, mPimpl->mAnalysisTypeX, sz, 0, 0, val);
			}
		}
		else
		{
			mPimpl->mAnalysisMethodX->getValues(mCrystalFp, val);
		}

		float minv =  FLT_MAX;
		float maxv = -FLT_MAX;
		for(i=0; i < sz; ++i)
		{
			if(val[i] > maxv) maxv = val[i];
			if(val[i] < minv) minv = val[i];
		}
		float delta = (maxv-minv)/mPimpl->mHistNumBins;
		aValues[0][0] = minv;
		for(i=1; i < mPimpl->mHistNumBins; ++i) aValues[0][i] = aValues[0][i-1]+delta;
		for(i=0; i < mPimpl->mHistNumBins; ++i) aValues[1][i] = 0;
		float tot = 0;
		for(i=0; i < sz; ++i)
		{
			unsigned int idx = (unsigned int)((val[i]-minv)*mPimpl->mHistNumBins/(maxv*(1.0000001)-minv));
			++aValues[1][idx];
			++tot;
		}
		for(i=0; i < mPimpl->mHistNumBins; ++i) aValues[1][i] /= tot;

		delete [] val;
		}
		break;

	case CATEGORY_SPECIAL:
		{
			size_t na = mPimpl->mAnalysisMethodX->numArrays();
			for(i=0; i < na; ++i)
			{
				if(mPimpl->mCacheResults)
				{
					size_t sz = mPimpl->mAnalysisMethodX->numElements(mCrystalFp, i);
					if(!mPimpl->mCache.retrieveResult(CATEGORY_SPECIAL, mPimpl->mAnalysisTypeX, sz, i, mPimpl->mStructureIdx, aValues[i]))
					{
						mPimpl->mAnalysisMethodX->getValues(mCrystalFp, aValues[i], i);

						mPimpl->mCache.addResult(CATEGORY_SPECIAL, mPimpl->mAnalysisTypeX, sz, i, mPimpl->mStructureIdx, aValues[i]);
					}
				}
				else
				{
					mPimpl->mAnalysisMethodX->getValues(mCrystalFp, aValues[i], i);
				}
			}
		}
		break;

	case CATEGORY_ALL:
		i = 0;
		for(AnalysisMethodsSimpleList::const_iterator ifml=mPimpl->mSimpleMethodsList.begin(); ifml != mPimpl->mSimpleMethodsList.end(); ++ifml)
		{
			// No caching. It is used only in batch
			if((*ifml)->numElements(mCrystalFp) != mCrystalFp->getNumActiveStructures()) continue;
			(*ifml)->getValues(mCrystalFp, aValues[i++]);
		}
		break;
	}
}


#if 0
bool CrystalFpAnalysis::SetAnalysisMethod(unsigned int category, unsigned int analysis_type_x, unsigned int analysis_type_y, const CrystalFpAnalysisParams& params)
{
	unsigned int nvx, nvy;
	std::vector<unsigned int> idx_out;
	unsigned int num_arrays;

	switch(category)
	{
	case CATEGORY_SIMPLE:
		pimpl->curr_analysis_type_x = analysis_type_x;
		pimpl->curr_analysis_type_y = analysis_type_y;

		pimpl->curr_method_x = pimpl->simple_methods_list->GetMethod(static_cast<size_t>(analysis_type_x));
		if(pimpl->curr_method_x == 0) return false;
		pimpl->curr_method_y = pimpl->simple_methods_list->GetMethod(static_cast<size_t>(analysis_type_y));
		if(pimpl->curr_method_y == 0) return false;

		pimpl->curr_method_x->SetupParameters(sl, params);
		pimpl->curr_method_y->SetupParameters(sl, params);
		if(!pimpl->curr_method_x->IsValid() || !pimpl->curr_method_y->IsValid()) return false;

		nvx = pimpl->curr_method_x->NumElements();
		nvy = pimpl->curr_method_y->NumElements();
		if(nvx != nvy) return false;
		pimpl->num_values = nvx;

		idx_out = pimpl->curr_method_x->GetValuesToOutput();
		pimpl->num_arrays = (unsigned int)idx_out.size();
		idx_out = pimpl->curr_method_y->GetValuesToOutput();
		pimpl->num_arrays += (unsigned int)idx_out.size();
		break;

	case CATEGORY_HIST:
		pimpl->curr_analysis_type_x = analysis_type_x;
		pimpl->curr_analysis_type_y = analysis_type_x;
		pimpl->curr_bin_size = params.num_x_bins;

		pimpl->curr_method_x = pimpl->simple_methods_list->GetMethod(static_cast<size_t>(analysis_type_x));
		if(pimpl->curr_method_x == 0) return false;

		pimpl->curr_method_x->SetupParameters(sl, params);
		if(!pimpl->curr_method_x->IsValid()) return false;

		pimpl->num_values = pimpl->curr_bin_size;
		pimpl->num_arrays = 2;
		break;

	case CATEGORY_SPECIAL:
		pimpl->curr_analysis_type_x = analysis_type_x;
		pimpl->curr_analysis_type_y = analysis_type_x;

		pimpl->curr_method_special = pimpl->special_methods_list->GetMethod(static_cast<size_t>(analysis_type_x));
		if(pimpl->curr_method_special == 0) return false;

		pimpl->curr_method_special->SetupParameters(sl, params);
		if(!pimpl->curr_method_special->IsValid()) return false;

		pimpl->num_values = pimpl->curr_method_special->NumElements();
		idx_out = pimpl->curr_method_special->GetValuesToOutput();
		pimpl->num_arrays = (unsigned int)idx_out.size();
		break;

	case CATEGORY_ALL:
		pimpl->num_values = (unsigned int)sl->GetNumActiveStructures();
		num_arrays = 0;
		for(MethodsList<AnalysisMethodSimple>::const_iterator ifml=pimpl->simple_methods_list->begin();
			ifml != pimpl->simple_methods_list->end();
			++ifml)
		{
			(*ifml)->SetupParameters(sl, params);

			if((*ifml)->NumElements() != pimpl->num_values) continue;
			num_arrays += (*ifml)->NumArrays();
		}
		pimpl->num_arrays = num_arrays;
		break;

	default:
		return false;
	}

	pimpl->curr_category = category;

	return true;
}


unsigned int CrystalFpAnalysis::GetNumValues(void) const
{
	return pimpl->num_values;
}


unsigned int CrystalFpAnalysis::GetNumArrays(void) const
{
	return pimpl->num_arrays;
}


unsigned int CrystalFpAnalysis::GetNumXValues(void) const
{
	std::vector<unsigned int> idx_out;

	switch(pimpl->curr_category)
	{
	case CATEGORY_SIMPLE:
		idx_out = pimpl->curr_method_x->GetValuesToOutput();
		return (unsigned int)idx_out.size();

	case CATEGORY_HIST:
	case CATEGORY_SPECIAL:
	case CATEGORY_ALL:
		return 1;
	}

	return 1;
}


const std::vector<std::string> CrystalFpAnalysis::GetLabels(void) const
{
	std::vector<std::string> labels;
	std::vector<std::string> tmp;

	switch(pimpl->curr_category)
	{
	case CATEGORY_SIMPLE:
		labels = pimpl->curr_method_x->GetLabels();
		tmp    = pimpl->curr_method_y->GetLabels();
		labels.insert(labels.end(), tmp.begin(), tmp.end());
		break;

	case CATEGORY_HIST:
		labels = pimpl->curr_method_x->GetLabels();
		labels.push_back("Frequency");
		break;

	case CATEGORY_SPECIAL:
		labels = pimpl->curr_method_special->GetLabels();
		break;

	case CATEGORY_ALL:
		for(MethodsList<AnalysisMethodSimple>::const_iterator ifml=pimpl->simple_methods_list->begin();
			ifml != pimpl->simple_methods_list->end();
			++ifml)
		{
			if((*ifml)->NumElements() != pimpl->num_values) continue;
			tmp = (*ifml)->GetLabels();
			labels.insert(labels.end(), tmp.begin(), tmp.end());
		}
		break;
	}

	return labels;
}


unsigned int CrystalFpAnalysis::GetValues(float **values) const
{
	unsigned int i, j, nx;
	std::vector<unsigned int> idx_out;
	std::vector<unsigned int>::const_iterator io;
	std::vector<unsigned int> id_offs;

	switch(pimpl->curr_category)
	{
	case CATEGORY_SIMPLE:
		idx_out = pimpl->curr_method_x->GetValuesToOutput();
		id_offs = pimpl->curr_method_x->GetCachingIdOffs();
		for(io=idx_out.begin(),i=0; io != idx_out.end(); ++io,++i)
		{
			if(values[i] == 0) continue;

			if(!pimpl->simple_values_cache.Retrieve(id_offs[*io] + pimpl->curr_analysis_type_x, pimpl->num_values, values[i]))
			{
				std::vector< std::vector<float> > vals = pimpl->curr_method_x->GetAllValues();
				std::vector<unsigned int>::const_iterator ioffs;
				for(ioffs=id_offs.begin(),j=0; ioffs != id_offs.end(); ++ioffs,++j)
				{
					pimpl->simple_values_cache.Add(*ioffs + pimpl->curr_analysis_type_x, vals[j]);
				}
				if(!pimpl->simple_values_cache.Retrieve(id_offs[*io] + pimpl->curr_analysis_type_x, pimpl->num_values, values[i]))
				{
					std::cerr << "Something wrong retrieving x " << *io << " for out value " << i << std::endl;
				}
			}
		}
		nx = i;
		idx_out = pimpl->curr_method_y->GetValuesToOutput();
		id_offs = pimpl->curr_method_y->GetCachingIdOffs();
		for(io=idx_out.begin(); io != idx_out.end(); ++io,++i)
		{
			if(values[i] == 0) continue;

			if(!pimpl->simple_values_cache.Retrieve(id_offs[*io] + pimpl->curr_analysis_type_y, pimpl->num_values, values[i]))
			{
				std::vector< std::vector<float> > vals = pimpl->curr_method_y->GetAllValues();
				std::vector<unsigned int>::const_iterator ioffs;
				for(ioffs=id_offs.begin(),j=0; ioffs != id_offs.end(); ++ioffs,++j)
				{
					pimpl->simple_values_cache.Add(*ioffs + pimpl->curr_analysis_type_y, vals[j]);
				}
				if(!pimpl->simple_values_cache.Retrieve(id_offs[*io] + pimpl->curr_analysis_type_y, pimpl->num_values, values[i]))
				{
					std::cerr << "Something wrong retrieving y " << *io << " for out value " << i << std::endl;
				}
			}
		}
		return nx;

	case CATEGORY_HIST:
		{
		idx_out = pimpl->curr_method_x->GetValuesToOutput();
		id_offs = pimpl->curr_method_x->GetCachingIdOffs();
		io = idx_out.begin();
		nx = pimpl->curr_method_x->NumElements();
		float *v = new float[nx];

		if(!pimpl->simple_values_cache.Retrieve(id_offs[*io] + pimpl->curr_analysis_type_x, nx, v))
		{
			std::vector< std::vector<float> > vals = pimpl->curr_method_x->GetAllValues();
			std::vector<unsigned int>::const_iterator ioffs;
			for(ioffs=id_offs.begin(),j=0; ioffs != id_offs.end(); ++ioffs,++j)
			{
				pimpl->simple_values_cache.Add(*ioffs + pimpl->curr_analysis_type_x, vals[j]);
			}
			if(!pimpl->simple_values_cache.Retrieve(id_offs[*io] + pimpl->curr_analysis_type_x, nx, v))
			{
				std::cerr << "Something wrong retrieving " << *io << " for histogramming" << std::endl;
			}
		}

		float minv = v[0];
		float maxv = v[0];

		for(i=1; i < nx; ++i)
		{
			if(v[i] < minv) minv = v[i];
			if(v[i] > maxv) maxv = v[i];
		}
		float delta = (maxv - minv)/pimpl->num_values;

		for(i=0; i < pimpl->num_values; ++i) values[0][i] = minv + i*delta;
		for(i=0; i < pimpl->num_values; ++i) values[1][i] = 0.0F;
		for(i=0; i < nx; ++i)
		{
			unsigned int nb = (unsigned int)((v[i]-minv)/delta);
			if(nb >= pimpl->num_values) nb = pimpl->num_values-1;
			values[1][nb] += 1.0F;
		}
		for(i=0; i < pimpl->num_values; ++i) values[1][i] /= (float)nx;

		delete [] v;
		}
		return 1;

	case CATEGORY_SPECIAL:
		idx_out = pimpl->curr_method_special->GetValuesToOutput();
		id_offs = pimpl->curr_method_special->GetCachingIdOffs();
		for(io=idx_out.begin(),i=0; io != idx_out.end(); ++io,++i)
		{
			if(values[i] == 0) continue;

			if(!pimpl->special_values_cache.Retrieve(id_offs[*io] + pimpl->curr_analysis_type_x, pimpl->num_values, values[i]))
			{
				std::vector< std::vector<float> > vals = pimpl->curr_method_special->GetAllValues();
				std::vector<unsigned int>::const_iterator ioffs;
				for(ioffs=id_offs.begin(),j=0; ioffs != id_offs.end(); ++ioffs,++j)
				{
					pimpl->special_values_cache.Add(*ioffs + pimpl->curr_analysis_type_x, vals[j]);
				}
				if(!pimpl->special_values_cache.Retrieve(id_offs[*io] + pimpl->curr_analysis_type_x, pimpl->num_values, values[i]))
				{
					std::cerr << "Something wrong retrieving special " << *io << " for out value " << i << std::endl;
				}
			}
		}
		return 1;

	case CATEGORY_ALL:
		j = 0;
		for(MethodsList<AnalysisMethodSimple>::const_iterator ifml=pimpl->simple_methods_list->begin();
			ifml != pimpl->simple_methods_list->end();
			++ifml)
		{
			if((*ifml)->NumElements() != pimpl->num_values) continue;
			std::vector< std::vector<float> > vals = (*ifml)->GetAllValues();
			for(std::vector< std::vector<float> >::const_iterator iv=vals.begin(); iv != vals.end(); ++iv)
			{
				memcpy(values[j], &((*iv)[0]), sizeof(float)*pimpl->num_values);
				++j;
			}
		}
		return 1;
	}

	return 0;
}


CrystalFpAnalysisParams::CrystalFpAnalysisParams()
{
	full_fingerprint = false;
	this_step = 0;
	this_part = 0;
	which_atom = 0;
	is_single_linkage = true;
	num_x_bins = 100;
}


CrystalFpAnalysisParams::CrystalFpAnalysisParams(const char* params)
{
	unsigned int i = 0;
	for(; params[i] && !isdigit(params[i]); ++i) {}
	if(params[i] == 0) return;
	full_fingerprint = atoi(params+i) != 0;
	for(; params[i] && isdigit(params[i]); ++i) {}
	for(; params[i] && !isdigit(params[i]); ++i) {}
	if(params[i] == 0) return;
	this_step = atoi(params+i);
	for(; params[i] && isdigit(params[i]); ++i) {}
	for(; params[i] && !isdigit(params[i]); ++i) {}
	if(params[i] == 0) return;
	this_part = atoi(params+i);
	for(; params[i] && isdigit(params[i]); ++i) {}
	for(; params[i] && !isdigit(params[i]); ++i) {}
	if(params[i] == 0) return;
	which_atom = atoi(params+i);
	for(; params[i] && isdigit(params[i]); ++i) {}
	for(; params[i] && !isdigit(params[i]); ++i) {}
	if(params[i] == 0) return;
	is_single_linkage = atoi(params+i) != 0;
	for(; params[i] && isdigit(params[i]); ++i) {}
	for(; params[i] && !isdigit(params[i]); ++i) {}
	if(params[i] == 0) return;
	num_x_bins = atoi(params+i);
}

#endif
