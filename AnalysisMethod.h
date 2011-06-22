
#ifndef ANALYSISMETHOD_H
#define ANALYSISMETHOD_H

#include <vector>

#include <string>
#include <cstring>
#include <list>
#include <climits>
#include "MethodsList.h"
#include "CrystalFp.h"

namespace cfp_internal
{

class AnalysisMethod
{
public:
	AnalysisMethod(const char* aName, const char* aNameHist=0)
	{
		mName = aName;
		if(aNameHist) mNameHist = aNameHist;
		mPartSelected = UINT_MAX; // All fingerprint parts selected
		mStructureIdx = 0;
	}
	virtual ~AnalysisMethod() {}

	std::string								getName(void) const { return mName; }
	std::string								getNameHist(void) const { return mNameHist; }
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const =0;
	virtual unsigned int					numElements(const cfp::CrystalFp* aCfp) const  { return aCfp->getNumActiveStructures(); }
	virtual unsigned int					numElements(const cfp::CrystalFp* aCfp, unsigned int /*aIdx*/) const  { return aCfp->getNumActiveStructures(); }
	virtual unsigned int					numArrays(void) const {return 1;}
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> names; return names; }
	virtual const std::vector<std::string>	getLabels(const cfp::CrystalFp* /*aCfp*/) const { std::vector<std::string> names; return names; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const =0;

	void									setFingerprintPart(unsigned int aPartSelected) {mPartSelected = aPartSelected;}
	void									setStructureIdx(unsigned int aStructureIdx) {mStructureIdx = aStructureIdx;}

protected:
	std::string			mName;
	std::string			mNameHist;
	unsigned int		mPartSelected;
	unsigned int		mStructureIdx;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetIdx : public AnalysisMethod
{
public:
	MethodGetIdx() : AnalysisMethod("Index") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 0; }
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Index"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetStep : public AnalysisMethod
{
public:
	MethodGetStep() : AnalysisMethod("Step") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 0; }
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Step"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetEnergy : public AnalysisMethod
{
public:
	MethodGetEnergy() : AnalysisMethod("Energy per atom", "Per atom energy histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 0 && aCfp->hasEnergies(); }
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Energy per atom"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetCellVolume : public AnalysisMethod
{
public:
	MethodGetCellVolume() : AnalysisMethod("Cell volume per atom", "Per atom cell volume histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 0 && aCfp->hasUnitCell(); }
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Cell volume per atom"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetDeltaEnergy : public AnalysisMethod
{
public:
	MethodGetDeltaEnergy() : AnalysisMethod("Per atom energy difference", "Per atom energy difference histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 1 && aCfp->hasEnergies(); }
	virtual unsigned int					numElements(const cfp::CrystalFp* aCfp) const { unsigned int n = aCfp->getNumActiveStructures(); return (unsigned int)((n*(n-1))/2);}
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Per atom energy difference"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetDistances : public AnalysisMethod
{
public:
	MethodGetDistances() : AnalysisMethod("Distances", "Distances histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 1 && aCfp->hasDistanceMatrix(); }
	virtual unsigned int					numElements(const cfp::CrystalFp* aCfp) const { unsigned int n = aCfp->getNumActiveStructures(); return (unsigned int)((n*(n-1))/2);}
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Distance"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetEnergyFromMin : public AnalysisMethod
{
public:
	MethodGetEnergyFromMin() : AnalysisMethod("Energy per atom from min.", "Energy per atom from min. histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 0 && aCfp->hasEnergies(); }
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Energy per atom from min."); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetDistFromMin : public AnalysisMethod
{
public:
	MethodGetDistFromMin() : AnalysisMethod("Distances from min.", "Distances from min. histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 0 && aCfp->hasDistanceMatrix(); }
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Distance from min."); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetPointDepth : public AnalysisMethod
{
public:
	MethodGetPointDepth() : AnalysisMethod("Point depth", "Point depth histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 1 && aCfp->hasFingerprints(); }
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Depth"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetOrderF2 : public AnalysisMethod
{
public:
	MethodGetOrderF2() : AnalysisMethod("Order (F2)", "Order (F2) histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const {
		return aCfp->getNumActiveStructures() > 1 && aCfp->hasFingerprints() && aCfp->isDiffractionLike() && aCfp->hasUnitCell();}
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Degree of order (F2)"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetOrderF2R2 : public AnalysisMethod
{
public:
	MethodGetOrderF2R2() : AnalysisMethod("Order (F2R2)", "Order (F2R2) histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const {
		return aCfp->getNumActiveStructures() > 1 && aCfp->hasFingerprints() && aCfp->isDiffractionLike() && aCfp->hasUnitCell();}
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Degree of order (F2R2)"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetQuasiEntropy : public AnalysisMethod
{
public:
	MethodGetQuasiEntropy() : AnalysisMethod("Quasi-entropy", "Quasi-entropy histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const {
		return aCfp->getNumActiveStructures() > 1 && aCfp->hasFingerprints() && aCfp->isDiffractionLike() && aCfp->hasUnitCell();}
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Quasi-entropy"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;

private:
	float computeQuasiEntropy(const cfp::CrystalFp* aCfp, unsigned int aIdx) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MethodGetMinDimension : public AnalysisMethod
{
public:
	MethodGetMinDimension() : AnalysisMethod("Min. dimension", "Min. dimension histogram") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 0  && aCfp->hasFingerprints() && aCfp->isDiffractionLike(); }
	virtual const std::vector<std::string>	getLabels(void) const { std::vector<std::string> lbl; lbl.push_back("Min. dimension"); return lbl; }
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class AnalysisMethodsSimpleList : public MethodsList<AnalysisMethod>
{
public:
	AnalysisMethodsSimpleList()
	{
		mMethodsList.push_back(new MethodGetIdx);
		mMethodsList.push_back(new MethodGetStep);
		mMethodsList.push_back(new MethodGetEnergy);
		mMethodsList.push_back(new MethodGetCellVolume);
		mMethodsList.push_back(new MethodGetDeltaEnergy);
		mMethodsList.push_back(new MethodGetDistances);
		mMethodsList.push_back(new MethodGetEnergyFromMin);
		mMethodsList.push_back(new MethodGetDistFromMin);
		mMethodsList.push_back(new MethodGetPointDepth);
		mMethodsList.push_back(new MethodGetOrderF2);
		mMethodsList.push_back(new MethodGetOrderF2R2);
		mMethodsList.push_back(new MethodGetQuasiEntropy);
		mMethodsList.push_back(new MethodGetMinDimension);
	}
	~AnalysisMethodsSimpleList() {clear();}
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class AnalysisMethodsHistList : public MethodsList<AnalysisMethod>
{
public:
	AnalysisMethodsHistList()
	{
		mMethodsList.push_back(new MethodGetEnergy);
		mMethodsList.push_back(new MethodGetCellVolume);
		mMethodsList.push_back(new MethodGetDeltaEnergy);
		mMethodsList.push_back(new MethodGetDistances);
		mMethodsList.push_back(new MethodGetEnergyFromMin);
		mMethodsList.push_back(new MethodGetDistFromMin);
		mMethodsList.push_back(new MethodGetPointDepth);
		mMethodsList.push_back(new MethodGetOrderF2);
		mMethodsList.push_back(new MethodGetOrderF2R2);
		mMethodsList.push_back(new MethodGetQuasiEntropy);
		mMethodsList.push_back(new MethodGetMinDimension);
	}
	~AnalysisMethodsHistList() {clear();}
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MethodGetFingerprint : public AnalysisMethod
{
public:
	MethodGetFingerprint() : AnalysisMethod("One fingerprint") {}
	virtual bool							isValid(const cfp::CrystalFp* aCfp) const { return aCfp->getNumActiveStructures() > 0 && aCfp->hasFingerprints(); }
	virtual unsigned int					numElements(const cfp::CrystalFp* aCfp, unsigned int /*aIdx*/) const
												{ return aCfp->getFingerprintSectionLen()*((mPartSelected >= aCfp->getFingerprintNumSections()) ? aCfp->getFingerprintNumSections() : 1);}
	virtual unsigned int					numArrays(void) const {return 2;}
	virtual const std::vector<std::string>	getLabels(const cfp::CrystalFp* aCfp) const;
	virtual void							getValues(const cfp::CrystalFp* aCfp, float* aValue, unsigned int aIdx=0) const;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class AnalysisMethodsSpecialList : public MethodsList<AnalysisMethod>
{
public:
	AnalysisMethodsSpecialList()
	{
		mMethodsList.push_back(new MethodGetFingerprint);
	}
	~AnalysisMethodsSpecialList() {clear();}
};

}

#endif

