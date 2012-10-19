/// @file CrystalFpAnalysis.h
/// Access to various analysis on a given CrystalFp object.
///
#ifndef CRYSTALFPANALYSIS_H
#define CRYSTALFPANALYSIS_H

#include "CrystalFp.h"
namespace cfp
{

/// Public interface to the analysis algorithms for CrystalFp results.
///
///	This is the interface to the algoritms needed to analyze the CrystalFp results. The results are line charts or scatterplots.
///
///	@author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///	@date 2008-06-06 (initial version)
///	@date 2009-09-01 (version 1.0)
///	@version 1.0
///
class CrystalFpAnalysis
{
public:
	/// Constructor
	/// @param aCfp The object describing the structures to be analyzed
	///
	CrystalFpAnalysis(const CrystalFp *aCfp);

	/// Destructor
	///
	~CrystalFpAnalysis();

	/// Reset
	/// Clean cache and reload the structure list under analysis
	/// @param aCfp The object describing the structures to be analyzed. If null the original object is not changed
	///
	void resetCache(const CrystalFp *aCfp=0);

	/// Return the names of the implemented distance computation methods
	/// @param aCategory Index of the category of analysis to use (simple, histogram, special)
	/// @return The list of distance measure methods names for the given category
	///
	const std::vector<std::string> getAnalysisMethodsNames(unsigned int aCategory) const;

	/// Analysis categories
	///
	enum CategoryTypes {
		CATEGORY_SIMPLE,	///< One value on x one value on y
		CATEGORY_HIST,		///< Histogram of one value
		CATEGORY_SPECIAL,	///< Other analysis included the 2D ones
		CATEGORY_ALL		///< All simple analysis
	};

	/// Set the analysis to perform
	///
	/// @param aAnalysisTypeX The index in the global list of methods to be used or the one to be used for the X value
	/// @param aAnalysisTypeY The index in the global list of methods to be used for the Y value (if any)
	/// @return Return false if the method index or the combination of methods is invalid or if the data cannot support the requested analysis method
	///
	bool setAnalysisMethod(unsigned int aAnalysisTypeX, unsigned int aAnalysisTypeY);

	/// Retrieve the length of each returned value plus the number of x values
	///
	/// @return A list with for first element the number of x variables, then the lengths of each variables, starting with the x ones
	///
	/// @exception CrystalFpFatal If analysis method has not been set
	///
	std::vector<size_t> numValues(void);

	/// Get the labels for the analysis result arrays
	/// @return The vector of labels
	///
	const std::vector<std::string> getLabels(void) const;

	/// Return the analysis values into preallocated arrays
	///
	/// @param aValues An array of pointers that point to arrays in which the results will be stored
	///
	void getValues(float **aValues) const;

	/// Set the parameter identified by the name.
	/// Mainly used to set parameters from the command line
	///
	/// @param aName Name of the parameter to be set (currently: "bin", "part", "idx")
	/// @param aValue String representation of the value to be set (internally it is converted to the needed data type)
	///
	/// @exception CrystalFpFatal If the param name is invalid
	///
	void setNamedParam(const std::string& aName, const std::string& aValue);

	/// Set the parameter identified by the name.
	/// Mainly used to set parameters from the command line
	///
	/// @param aName Name of the parameter to be set (currently: "bin", "part", "idx")
	/// @param aValue String representation of the value to be set (internally it is converted to the needed data type)
	///
	/// @exception CrystalFpFatal If the param name is invalid
	///
	void setNamedParam(const std::string& aName, unsigned int aValue);

	/// Enables the results caching.
	/// After this call all the computed values are retained in a cache and reused if needed.
	///
	void enableCaching(void);


private:
	const CrystalFp* mCrystalFp;

private:
	struct CrystalFpAnalysisImpl;
	struct CrystalFpAnalysisImpl* mPimpl;
};

}
#endif

