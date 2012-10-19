
#ifndef CRYSTALFP_METHODS_LIST_H
#define CRYSTALFP_METHODS_LIST_H

#include <vector>
#include "CrystalFpExceptions.h"

namespace cfp_internal
{

template<typename T>
class MethodsList
{
public:
	T *getMethod(unsigned int aMethodIdx) const
	{
		if(aMethodIdx >= mMethodsList.size()) throw cfp::CrystalFpFatal("Invalid index in getMethod");
		return mMethodsList[aMethodIdx];
	}

	typedef typename std::vector<T *>::iterator iterator;
	typedef typename std::vector<T *>::const_iterator const_iterator;

	iterator begin(void)			 {return mMethodsList.begin();}
	iterator end(void)				 {return mMethodsList.end();}
	const_iterator begin(void) const {return mMethodsList.begin();}
	const_iterator end(void) const   {return mMethodsList.end();}

	void clear(void)
	{
		typename std::vector<T *>::iterator im;
		for(im=mMethodsList.begin(); im != mMethodsList.end(); ++im) delete (*im);
		mMethodsList.clear();
	}
	size_t size(void) { return mMethodsList.size(); }


protected:
	std::vector<T *> mMethodsList;
};


}
#endif

