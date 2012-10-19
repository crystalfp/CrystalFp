
#ifndef FINGERPRINTINGMETHOD_H
#define FINGERPRINTINGMETHOD_H

#include <string>
#include "MethodsList.h"
#include "Structure.h"

namespace cfp_internal
{

class FingerprintingMethod
{
public:
	FingerprintingMethod(const char* aName) : mName(aName), mCutoff(10.0F), mBinSize(0.05F), mPeakSize(0.02F), mIsNanocluster(false)
	{
		// mName			= aName;
		// mCutoff			= 10.0F;
		// mBinSize		= 0.05F;
		// mPeakSize		= 0.02F;
		// mIsNanocluster	= false;
	};
	const std::string&  getName(void) const { return mName; }
	virtual      ~FingerprintingMethod() {};
	virtual void computeFingerprint(Structure& aStructure, const unsigned int* aExpansion) =0;
	void		 setCutoffDistance(float aCutoff) {mCutoff = aCutoff;}
	void		 setNanoclusterStructureType(void) {mIsNanocluster = true;}
	void		 setDiffrBinSize(float aBinSize) {mBinSize = aBinSize;}
	void		 setDiffrPeakSize(float aPeakSize) {mPeakSize = aPeakSize;}
	virtual bool isDiffractionLike(void) const {return false;}

protected:
	std::string mName;
	float		mCutoff;
	float		mBinSize;
	float		mPeakSize;
	bool		mIsNanocluster;
};

//class NormalizedDiffractionHistogram : public FingerprintingMethod
//{
//public:
//	NormalizedDiffractionHistogram() : FingerprintingMethod("Normalized diffraction") {}
//	virtual void computeFingerprint(Structure& aStructure, const unsigned int* aExpansion);
//	virtual bool isDiffractionLike(void) const {return true;}
//};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class PerElementRdfHistogram : public FingerprintingMethod
{
public:
	PerElementRdfHistogram() : FingerprintingMethod("Per element diffraction") {}
	virtual void computeFingerprint(Structure& aStructure, const unsigned int* aExpansion);
	virtual bool isDiffractionLike(void) const {return true;}
};

class FingerprintingMethodsList : public MethodsList<FingerprintingMethod>
{
public:
	FingerprintingMethodsList()
	{
		//mMethodsList.push_back(new NormalizedDiffractionHistogram);
		mMethodsList.push_back(new PerElementRdfHistogram);
	}
	~FingerprintingMethodsList() {clear();}
};

}
#endif

