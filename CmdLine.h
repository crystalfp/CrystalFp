
#ifndef CMDLINE_H
#define CMDLINE_H

#include <vector>
#include <map>
#include <string>
#include <climits>
#include "simpleopt/SimpleOpt.h"

/// Parse the command line flags.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-22 (initial version)
///     @version 1.0
///
///
class CmdLine
{
public:
	/// Constructor.
	///
	/// It parses the command line throwing exceptions on error
	///
	/// @param[in] aCnt Number of command line parameters (from main argc)
	/// @param[in] aVal Array of command line parameters (from main argv)
	/// @param[in] aUseDefaultsIfNoArguments If true and aCnt < 2 use hardcoded default values
	///
	CmdLine(int aCnt, char **aVal, bool aUseDefaultsIfNoArguments=false);

public:
	unsigned int				mVerboseLevel;				///< Verbose level (0: no verbose, 1: mildly verbose, 2: verbose)
	std::vector<unsigned int>	mAtomZ;						///< List of element types
	const char*					mPoscarFile;				///< POSCAR structure file name
	const char*					mEnergyFile;				///< Energy file name
	unsigned int				mStartStep;					///< First step to be loaded (steps starts from 1)
	unsigned int				mEndStep;					///< Last step to be loaded (zero means load all)
	bool						mEnergyIsPerAtom;			///< If the energy is per atom. Otherwise it is per structure
	float						mEnergyThreshold;			///< Energy threshold (ignored if mHasEnergyThreshold is false)
	bool						mHasEnergyThreshold;		///< If there is an energy threshold
	bool						mHasReverseEnergyThreshold;	///< If the energy threshold is from the minimum energy in the set
	float						mCutoffDistance;			///< Cutof distance for fingerprint building (zero means no cutoff)	
	bool						mIsNanocluster;				///< The structures are not crystals
	float						mDiffrBinSize;				///< For diffraction-like methods gives the bin width for the histogram
	float						mDiffrPeakSize;				///< For diffraction-like methods gives the gaussian smoothing peak width
	unsigned int				mForcedFpLen;				///< If not zero, forces the total length of the fingerprint to such value
	const char*					mCheckpointDir;				///< Path of a directory to save the fingerprints as soon as they are computed (if unset do not checkpoint)
	bool						mOverwriteChkptDir;			///< Instead of reloading existing checkpoint, overwrite it
	bool						mListFingerprintingMethods;	///< List fingerprint methods and exit
	unsigned int				mFingerprintingMethod;		///< Fingerprinting method to be used (index same order as the list)
	bool						mListDistanceMethods;		///< List distance measure methods and exit
	unsigned int				mDistanceMethod;			///< Distance measure method to be used (index same order as the list)
	bool						mListGroupingMethods;		///< List grouping methods and exit
	unsigned int				mGroupingMethod;			///< Grouping method to be used (index same order as the list)
	unsigned int				mK;							///< K value needed by some grouping algorithm
	float						mMaxDistanceForGrouping;	///< Grouping max distance threshold
	const char*					mSummaryFile;				///< If set output a short summary file amenable ro be read back by a machine
	const char*					mFldOutFp;					///< If set output fingerprints as AVS FLD format to this file
	const char*					mFldOutDist;				///< If set output distance matrix as AVS FLD format to this file
	const char*					mSortedDistFile;			///< If set output sorted distances to this file
	const char*					mMapFile;					///< If duplicated are removed this file maps the original step to the new one
	bool						mRemoveDuplicates;			///< If set the duplicated elements are reduced to only one representative per group
	bool						mListAnalysisMethods;
	unsigned int				mAnalysisMethod;			///< Analysis to be applied on the fingerprints results
	unsigned int				mAnalysisMethod2;			///< Analysis to be applied on the fingerprints results (second value if needed)
	const char*					mAnalysisFile;				///< Output file for analysis results
	std::map<std::string,std::string> mAnalysisParams;		///< Pairs of parameter-name parameter-value for the analysis algorithm
	bool						mCreateScatterplot;
	std::map<std::string,std::string> mScatterplotParams;	///< Pairs of parameter-name parameter-value for the analysis algorithm
	const char*					mScatterplotFile;			///< Output file for the scatterplot points
	const char*					mDiagnosticFile;			///< Output file for the diagnostic chart for the current scatterplot


	static const unsigned int NO_METHOD_SELECTED = UINT_MAX;

private:

	/// Return the text corresponding to an error code.
	///
	/// @param[in] aOptParser The command line parser object
	/// @return The human readable error message
	///
	const char *getLastErrorText(CSimpleOpt& aOptParser);

	/// Print the help about the parameters.
	/// @param[in] aParserOptions The table of options definitions
	///
	void showHelp(const CSimpleOpt::SOption *aParserOptions);
};


#include <stdexcept>


/// Fatal error in CrystalFp.
/// The message explain the reason
///
class CmdLineFatal : public std::runtime_error
{
public:
	CmdLineFatal(const char *msg="") : runtime_error(msg)
	{}
};

/// Early termination message.
/// It does not signal a fatal error, but something like printing help.
/// It is never throw by the library, but it is used by the driver program.
///
class CmdLineSuccess : public std::exception
{
public:
	CmdLineSuccess() : exception()
	{}
};


#endif

