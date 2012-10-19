
#include <iostream>
#include <cctype>
#include <cstdlib>
#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif
#include "simpleopt/SimpleOpt.h"
#include "CmdLine.h"
#include "AtomSymbols.h"

const char *CmdLine::getLastErrorText(CSimpleOpt& aOptParser)
{
    switch(aOptParser.LastError())
	{
    case SO_SUCCESS:            return "Success";
    case SO_OPT_INVALID:        return "Unrecognized option";
    case SO_OPT_MULTIPLE:       return "Option matched multiple strings";
    case SO_ARG_INVALID:        return "Option does not accept argument";
    case SO_ARG_INVALID_TYPE:   return "Invalid argument format";
    case SO_ARG_MISSING:        return "Required argument is missing";
    case SO_ARG_INVALID_DATA:   return "Invalid argument data";
    default:					return "Unknown error";
    }
}


void CmdLine::showHelp(const CSimpleOpt::SOption *aParserOptions)
{
	size_t i, j, cnt;

	// Count entries and create an indicator array
	for(cnt=0; aParserOptions[cnt].pszArg != NULL; ++cnt) {}
	bool *done = new bool[cnt];
	for(i=0; i < cnt; ++i) done[i] = false;

	// For each different option
	for(i=0; i < cnt; ++i)
	{
		if(done[i]) continue;
		done[i] = true;

		std::cerr << aParserOptions[i].pszArg;
		for(j=i+1; j < cnt; ++j)
		{
			if(done[j] || aParserOptions[j].nId != aParserOptions[i].nId) continue;
			done[j] = true;
			std::cerr << "  " << aParserOptions[j].pszArg;
		}

		// Translate the kind of argument
		const char* type = "";
		switch(aParserOptions[i].nArgType)
		{
		case SO_NONE:   
			type = "(no argument)";
			break;

		case SO_REQ_SEP:
		case SO_REQ_CMB:
			type = "(required argument)";
			break;

		case SO_OPT:
			type = "(optional argument)";
			break;

		case SO_MULTI:
			type = "(multiple arguments)";
			break;
		}

		std::cerr << " " << type << std::endl;
		std::cerr << "        " << aParserOptions[i].pszHelp << std::endl << std::endl;
	}

	delete [] done;
}



CmdLine::CmdLine(int aCnt, char **aVal, bool aUseDefaultsIfNoArguments)
{
	// To support debugging these values are set if no arguments on the command line
	if(aUseDefaultsIfNoArguments && aCnt < 2 && getenv("CFP_DEBUG"))
	{
#ifdef _MSC_VER
		mPoscarFile					= "C:/mv/cfp/NewCrystalFpLib/gaas-8at_new.log";
		mEnergyFile					= "C:/mv/cfp/NewCrystalFpLib/gaas-8at_new.log.energy";
		mCheckpointDir				= "C:/mv/cfp/NewCrystalFpLib/test";
#else
		mPoscarFile					= "/users/mvalle/cfp/NewCrystalFpLib/gaas-8at_new.log";
		mEnergyFile					= "/users/mvalle/cfp/NewCrystalFpLib/gaas-8at_new.log.energy";
		mCheckpointDir				= "/users/mvalle/cfp/NewCrystalFpLib/test";
#endif
		mVerboseLevel				= 2;
		mStartStep					= 3;
		mEndStep					= 20;
		mEnergyIsPerAtom			= false;
		mEnergyThreshold			= 0;	
		mHasEnergyThreshold			= false;
		mHasReverseEnergyThreshold	= false;
		mCutoffDistance				= 0;	
		mIsNanocluster				= false;
		mDiffrBinSize				= 0.05F;
		mDiffrPeakSize				= 0.02F;
		mForcedFpLen				= 0;	
		mOverwriteChkptDir			= true;
		mListFingerprintingMethods	= false;
		mFingerprintingMethod		= 0;
		mListDistanceMethods		= false;
		mDistanceMethod				= 0;
		mListGroupingMethods		= false;
		mGroupingMethod				= 0;
		mK							= 2;
		mMaxDistanceForGrouping		= 0.01F;
		mSummaryFile				= 0;
		mFldOutFp					= 0;
		mFldOutDist					= 0;
		mSortedDistFile				= 0;
		mMapFile					= 0;
		mRemoveDuplicates			= false;
		mListAnalysisMethods		= false;
		mAnalysisMethod				= CmdLine::NO_METHOD_SELECTED;
		mAnalysisMethod2			= CmdLine::NO_METHOD_SELECTED;
		mAnalysisFile				= 0;
		mCreateScatterplot			= false;
		mScatterplotFile			= 0;
		mDiagnosticFile				= 0;
		mSerializeFile				= 0;

		return;
	}

	// Initialize options values
	mVerboseLevel				= 0;
	mPoscarFile					= 0;
	mEnergyFile					= 0;
	mStartStep					= 1;
	mEndStep					= 0;
	mEnergyIsPerAtom			= true;
	mEnergyThreshold			= 0.0F;
	mHasEnergyThreshold			= false;
	mHasReverseEnergyThreshold	= false;
	mCutoffDistance				= 0;
	mIsNanocluster				= false;
	mDiffrBinSize				= 0.05F;
	mDiffrPeakSize				= 0.02F;
	mForcedFpLen				= 0;
	mCheckpointDir				= 0;
	mOverwriteChkptDir			= false;
	mListFingerprintingMethods	= false;
	mFingerprintingMethod		= CmdLine::NO_METHOD_SELECTED;
	mListDistanceMethods		= false;
	mDistanceMethod				= CmdLine::NO_METHOD_SELECTED;
	mListGroupingMethods		= false;
	mGroupingMethod				= CmdLine::NO_METHOD_SELECTED;
	mK							= 0;
	mMaxDistanceForGrouping		= 0.01F;
	mSummaryFile				= 0;
	mFldOutFp					= 0;
	mFldOutDist					= 0;
	mSortedDistFile				= 0;
	mMapFile					= 0;
	mRemoveDuplicates			= false;
	mListAnalysisMethods		= false;
	mAnalysisMethod				= CmdLine::NO_METHOD_SELECTED;
	mAnalysisMethod2			= CmdLine::NO_METHOD_SELECTED;
	mAnalysisFile				= 0;
	mCreateScatterplot			= false;
	mScatterplotFile			= 0;
	mDiagnosticFile				= 0;
	mSerializeFile				= 0;

	// Setup the command line parser
	enum {
		OPT_VERSION,
		OPT_VERBOSE,
		OPT_HELP,
		OPT_TYPES,
		OPT_MAX_STEP,
		OPT_MIN_STEP,
		OPT_ENERGY_PER_STR,
		OPT_THRESH,
		OPT_REV_THRESH,
		OPT_CUTOFF_DIST,
		OPT_NANO_CLUST,
		OPT_BIN_SIZE,
		OPT_PEAK_SIZE,
		OPT_FORCE_DIM,
		OPT_CHKPT_DIR,
		OPT_IGNORE_CHKPT,
		OPT_FP_METHOD,
		OPT_DIST_METHOD,
		OPT_GROUPING_METHOD,
		OPT_GROUPING_THRESH,
		OPT_K_VALUE,
		OPT_SUMMARY,
		OPT_FLD_FP,
		OPT_FLD_DIST,
		OPT_DIST_DISTRIB,
		OPT_REMOVE_DUPL,
		OPT_ANALYSIS,
		OPT_ANALYSIS_FILE,
		OPT_ANALYSIS_PARAM,
		OPT_SCATTER,
		OPT_SCATTER_FILE,
		OPT_SCATTER_PARAM,
		OPT_DIAGNOSTIC_FILE,
		OPT_SERIALIZE
	};

	CSimpleOpt::SOption parser_options[] = {
		{ OPT_VERSION,			"--version",				SO_NONE, 	"Library and driver version" },
		{ OPT_VERBOSE,			"-v",						SO_OPT, 	"Verbose level (if no argument, defaults to 1)" },
		{ OPT_VERBOSE,			"--verbose",				SO_OPT, 	"" },
		{ OPT_HELP,				"-?",						SO_NONE,    "This help" },
		{ OPT_HELP,				"-h",						SO_NONE,    "" },
		{ OPT_HELP,				"--help",					SO_NONE,    "" },
		{ OPT_TYPES,			"-t",						SO_REQ_SEP, "List of chemical elements to be assigned to POSCAR atoms" },
		{ OPT_TYPES,			"--elements",				SO_REQ_SEP, "" },
		{ OPT_MAX_STEP,			"-es",						SO_REQ_SEP, "Last step to load (default: all)" },
		{ OPT_MAX_STEP,			"--max-step",				SO_REQ_SEP, ""},
		{ OPT_MAX_STEP,			"--end-step",				SO_REQ_SEP, ""},
		{ OPT_MIN_STEP,			"-ss",						SO_REQ_SEP, "First step to load (default: first)" },
		{ OPT_MIN_STEP,			"--start-step",				SO_REQ_SEP, ""},
		{ OPT_ENERGY_PER_STR,	"-et",						SO_NONE,    "Energy from file is per structure, not per atom" },
		{ OPT_ENERGY_PER_STR,	"--energy-per-structure",	SO_NONE,    "" },
		{ OPT_THRESH,			"-e",						SO_REQ_SEP, "Energy threshold" },
		{ OPT_THRESH,			"--energy-threshold",		SO_REQ_SEP, "" },
		{ OPT_REV_THRESH,		"-r",						SO_REQ_SEP, "Threshold from minimum energy" },
		{ OPT_REV_THRESH,		"--threshold-from-min",		SO_REQ_SEP, "" },
		{ OPT_CUTOFF_DIST,		"-c",						SO_REQ_SEP, "Fingerprint forced cutoff distance" },
		{ OPT_CUTOFF_DIST,		"--cutoff-distance",		SO_REQ_SEP, "" },
		{ OPT_NANO_CLUST,		"-n",						SO_NONE,    "The structures are nanoclusters, not crystals" },
		{ OPT_NANO_CLUST,		"--nano-clusters",			SO_NONE,    "" },
		{ OPT_NANO_CLUST,		"--nanoclusters",			SO_NONE,    "" },
		{ OPT_BIN_SIZE,			"-b",						SO_REQ_SEP, "Bin size for the pseudo-diffraction methods" },
		{ OPT_BIN_SIZE,			"--bin-size",				SO_REQ_SEP, "" },
		{ OPT_PEAK_SIZE,		"-p",						SO_REQ_SEP, "Peak smearing size" },
		{ OPT_PEAK_SIZE,		"--peak-size",				SO_REQ_SEP, "" },
		{ OPT_FORCE_DIM,		"-x",						SO_REQ_SEP,	"Force to this value the total length of the fingerprint (valid only on reloading from checkpoint)" },
		{ OPT_FORCE_DIM,		"--force-dim",				SO_REQ_SEP, "" },
		{ OPT_CHKPT_DIR,		"-chk",						SO_REQ_SEP, "Set fingerprint checkpoint directory" },
		{ OPT_CHKPT_DIR,		"--checkpoint-dir",			SO_REQ_SEP, "" },
		{ OPT_IGNORE_CHKPT,		"-i",						SO_NONE,    "Ignore previous checkpoint" },
		{ OPT_IGNORE_CHKPT,		"--ignore-checkpoint",		SO_NONE,    "" },
		{ OPT_FP_METHOD,		"-f",						SO_REQ_SEP, "Compute fingerprints using the given method ('help' to list all methods)" },
		{ OPT_FP_METHOD,		"--fingerprint-method",		SO_REQ_SEP, "" },
		{ OPT_DIST_METHOD,		"-d",						SO_REQ_SEP, "Compute distances using the given method ('help' to list all methods)"},
		{ OPT_DIST_METHOD,		"--distance-method",		SO_REQ_SEP, ""},
		{ OPT_GROUPING_METHOD,	"-g",						SO_REQ_SEP, "Grouping method ('help' to list all methods)" },
		{ OPT_GROUPING_METHOD,	"--grouping-method",		SO_REQ_SEP, "" },
		{ OPT_GROUPING_THRESH,	"-gt",						SO_REQ_SEP, "Grouping max distance threshold" },
		{ OPT_GROUPING_THRESH,	"--grouping-threshold",		SO_REQ_SEP, "" },
		{ OPT_K_VALUE,			"-k",						SO_REQ_SEP, "K value needed by some grouping algorithm" },
		{ OPT_K_VALUE,			"--k",						SO_REQ_SEP, "" },
		{ OPT_SUMMARY,			"-s", 						SO_REQ_SEP, "Output a summary file to be used by another program" },
		{ OPT_SUMMARY,			"--summary", 				SO_REQ_SEP, "" },
		{ OPT_FLD_FP,			"-fo",						SO_REQ_SEP, "Output fingerprints in FLD AVS format" },
		{ OPT_FLD_FP,			"--fld-fingerprints",		SO_REQ_SEP, "" },
		{ OPT_FLD_DIST,			"-do",						SO_REQ_SEP, "Output distance matrix in FLD AVS format" },
		{ OPT_FLD_DIST,			"--fld-distances",			SO_REQ_SEP, "" },
		{ OPT_DIST_DISTRIB,		"-sd",						SO_REQ_SEP, "Output distances in increasing order" },
		{ OPT_DIST_DISTRIB,		"--sorted-distances",		SO_REQ_SEP, "" },
		{ OPT_REMOVE_DUPL,		"-rd",						SO_OPT,     "Remove duplicates substituting them with a representative structure and writing the map of indices to file" },
		{ OPT_REMOVE_DUPL,		"--remove-dupl",			SO_OPT,     "" },
		{ OPT_ANALYSIS,			"-a",						SO_REQ_SEP,	"Create a chart for the given analysis method (can be: ax,ay or ax)" },
		{ OPT_ANALYSIS,			"--analysis",				SO_REQ_SEP,	"" },
		{ OPT_ANALYSIS_FILE,	"-af",						SO_REQ_SEP, "File output for analysis" },
		{ OPT_ANALYSIS_FILE,	"--analysis-file",			SO_REQ_SEP, "" },
		{ OPT_ANALYSIS_PARAM,	"-ap",						SO_MULTI,	"Pass one parameter for the given analysis method as: code value (code: bins, part, idx)" },
		{ OPT_ANALYSIS_PARAM,	"--analysis-param",			SO_MULTI,	"" },
		{ OPT_SCATTER,			"-sc",						SO_NONE,    "Create a scatterplot" },
		{ OPT_SCATTER,			"--scatterplot",			SO_NONE,    "" },
		{ OPT_SCATTER_FILE,		"-sf",						SO_REQ_SEP, "File output for scatterplot" },
		{ OPT_SCATTER_FILE,		"--scatterplot-file",		SO_REQ_SEP, "" },
		{ OPT_DIAGNOSTIC_FILE,	"-df",						SO_REQ_SEP, "File output for diagnostic chart assocated to the scatterplot" },
		{ OPT_DIAGNOSTIC_FILE,	"--diagnostic-file",		SO_REQ_SEP, "" },
		{ OPT_SCATTER_PARAM,	"-sp",						SO_MULTI,	"Pass one parameter for the scatterplot as: code value (code: retry, energy, iterations, kind, diagnostic, timestep, mass, stiffness, damping, perturb, bins, wobble)" },
		{ OPT_SCATTER_PARAM,	"--scatterplot-param",		SO_MULTI,	"" },
		{ OPT_SERIALIZE,		"-sa",						SO_REQ_SEP, "Save the current CrystalFp class to the given file" },
		{ OPT_SERIALIZE,		"--save-file",				SO_REQ_SEP, "" },
		SO_END_OF_OPTIONS
	};

	// Setup the usage string (to be concatenated with argv[0])
	const char* usage_msg = " [options] POSCARfile [ENERGIESfile]";

    // Declare our options parser, pass in the arguments from main as well as our array of valid options.
    CSimpleOpt args(aCnt, aVal, parser_options, SO_O_NOSLASH|SO_O_EXACT);

    // While there are arguments left to process
    while(args.Next())
	{
		if(args.LastError() == SO_OPT_INVALID)
		{
			std::cerr << "Error: " << getLastErrorText(args) << ": " << args.OptionText() << std::endl;
			throw CmdLineFatal();
        }
        if(args.LastError() != SO_SUCCESS)
		{
            std::cerr << "Error: " << getLastErrorText(args) << " for: " << args.OptionText() << std::endl;
			throw CmdLineFatal();
        }

		int tmp;
		const char* list_methods_tmp;
		char **arg_param;

		switch(args.OptionId())
		{
		case OPT_VERSION:
			// Get the last commit date from Git using: git log -n -pretty=format:"%ci"
			std::cerr << "CrystalFp library and driver: " << "2011-07-03 06:59:30 +0200" << std::endl;
			throw CmdLineSuccess();

		case OPT_VERBOSE:
			if(args.OptionArg())
			{
				tmp = atoi(args.OptionArg());
				if(tmp < 0) std::cerr << "Invalid value for " << args.OptionText() << " using 0" << std::endl;
				mVerboseLevel = static_cast<unsigned int>(tmp);
			}
			else mVerboseLevel = 1;
			break;

		case OPT_HELP:
			std::cerr << "Usage:" << std::endl;
			std::cerr << "    " << aVal[0] << usage_msg << std::endl << std::endl;
			showHelp(parser_options);
			throw CmdLineSuccess();

		case OPT_TYPES:
			convertAtomSymbolsToZ(args.OptionArg(), mAtomZ);
			break;

		case OPT_MIN_STEP:
			tmp = atoi(args.OptionArg());
			if(tmp < 1) {std::cerr << "Invalid value for " << args.OptionText() << " using 1" << std::endl; tmp = 1;}
			mStartStep = static_cast<unsigned int>(tmp);
			break;

		case OPT_MAX_STEP:
			tmp = atoi(args.OptionArg());
			if(tmp < 1) {std::cerr << "Invalid value for " << args.OptionText() << " using 0 (all steps)" << std::endl; tmp = 0;}
			mEndStep = static_cast<unsigned int>(tmp);
			break;

		case OPT_ENERGY_PER_STR:
			mEnergyIsPerAtom = false;
			break;

		case OPT_THRESH:
			mEnergyThreshold = (float)atof(args.OptionArg());
			mHasEnergyThreshold = true;
			break;

		case OPT_REV_THRESH:
			mEnergyThreshold = (float)atof(args.OptionArg());
			mHasReverseEnergyThreshold = true;
			break;

		case OPT_CUTOFF_DIST:
			mCutoffDistance = (float)atof(args.OptionArg());
			break;

		case OPT_NANO_CLUST:
			mIsNanocluster = true;
			break;

		case OPT_BIN_SIZE:
			mDiffrBinSize = (float)atof(args.OptionArg());
			break;

		case OPT_PEAK_SIZE:
			mDiffrPeakSize = (float)atof(args.OptionArg());
			break;

		case OPT_FORCE_DIM:
			tmp = atoi(args.OptionArg());
			if(tmp < 0) {std::cerr << "Invalid value for " << args.OptionText() << " ignoring" << std::endl; tmp = 0;}
			mForcedFpLen = static_cast<unsigned int>(tmp);
			break;

		case OPT_CHKPT_DIR:
			mCheckpointDir = args.OptionArg();
			break;

		case OPT_IGNORE_CHKPT:
			mOverwriteChkptDir = true;
			break;

		case OPT_FP_METHOD:
			list_methods_tmp = args.OptionArg();
			if(!strcasecmp(list_methods_tmp, "help"))
				mListFingerprintingMethods = true;
			else
				mFingerprintingMethod = atoi(list_methods_tmp);
			break;

		case OPT_DIST_METHOD:
			list_methods_tmp = args.OptionArg();
			if(!strcasecmp(list_methods_tmp, "help"))
				mListDistanceMethods = true;
			else
				mDistanceMethod = atoi(list_methods_tmp);
			break;

		case OPT_GROUPING_METHOD:
			list_methods_tmp = args.OptionArg();
			if(!strcasecmp(list_methods_tmp, "help"))
				mListGroupingMethods = true;
			else
				mGroupingMethod = atoi(list_methods_tmp);
			break;

		case OPT_GROUPING_THRESH:
			mMaxDistanceForGrouping = (float)atof(args.OptionArg());
			break;

		case OPT_K_VALUE:
			tmp = atoi(args.OptionArg());
			if(tmp < 0) {std::cerr << "Invalid value for " << args.OptionText() << " ignoring" << std::endl; tmp = 0;}
			mK = static_cast<unsigned int>(tmp);
			break;

		case OPT_SUMMARY:
			mSummaryFile = args.OptionArg();
			break;

		case OPT_FLD_FP:
			mFldOutFp = args.OptionArg();
			break;

		case OPT_FLD_DIST:
			mFldOutDist = args.OptionArg();
			break;

		case OPT_DIST_DISTRIB:
			mSortedDistFile = args.OptionArg();
			break;

		case OPT_REMOVE_DUPL:
			if(args.OptionArg()) mMapFile = args.OptionArg();
			mRemoveDuplicates = true;
			break;

		case OPT_ANALYSIS:
			list_methods_tmp = args.OptionArg();
			if(!strcasecmp(list_methods_tmp, "help"))
			{
				mListAnalysisMethods = true;
			}
			else
			{
				mAnalysisMethod = atoi(list_methods_tmp);
				int i;
				for(i=0; isdigit(list_methods_tmp[i]); ++i) {}
				for(; list_methods_tmp[i] && !isdigit(list_methods_tmp[i]); ++i) {}
				mAnalysisMethod2 = (list_methods_tmp[i]) ? atoi(list_methods_tmp+i) : CmdLine::NO_METHOD_SELECTED;
			}
			break;

		case OPT_ANALYSIS_PARAM:
			arg_param = args.MultiArg(2);
			if(!arg_param || !strncasecmp(arg_param[0], "help", 1))
			{
				std::cerr << "bin       " << "Number of bins for histograms" << std::endl;
				std::cerr << "parts     " << "Fingerprint part to be show" << std::endl;
				std::cerr << "idx       " << "Structure index" << std::endl;

				throw CmdLineSuccess();
			}
			mAnalysisParams.insert(std::pair<std::string,std::string>(arg_param[0], arg_param[1]));
			break;

		case OPT_ANALYSIS_FILE:
			mAnalysisFile = args.OptionArg();
			break;

		case OPT_SCATTER:
			mCreateScatterplot = true;
			break;

		case OPT_SCATTER_FILE:
			mScatterplotFile = args.OptionArg();
			break;

		case OPT_SCATTER_PARAM:
			arg_param = args.MultiArg(2);
			if(!arg_param || !strncasecmp(arg_param[0], "help", 1))
			{
				std::cerr << "retry       " << "Number of retries" << std::endl;
				std::cerr << "mass        " << "Ball mass" << std::endl;
				std::cerr << "stiffness   " << "Of the spring" << std::endl;
				std::cerr << "damping     " << "Damping factor for the movement" << std::endl;
				std::cerr << "perturb     " << "Perturb scale for the retries (it perturb the initial position of the masses)" << std::endl;
				std::cerr << "bins        " << "Number of bins for the binned distances diagnostics" << std::endl;
				std::cerr << "wobble      " << "Wobble scale for the position of the binned points" << std::endl;
				std::cerr << "energy      " << "Max kinetic energy to end iterations" << std::endl;
				std::cerr << "iterations  " << "Max number of iterations" << std::endl;
				std::cerr << "kind        " << "Kind of data to be output" << std::endl;
				std::cerr << "diagnostic  " << "Kind of diagnostic chart to produce" << std::endl;
				std::cerr << "timestep    " << "Timestep for the iterations" << std::endl;

				throw CmdLineSuccess();
			}
			mScatterplotParams.insert(std::pair<std::string,std::string>(arg_param[0], arg_param[1]));
			break;

		case OPT_DIAGNOSTIC_FILE:
			mDiagnosticFile = args.OptionArg();
			break;

		case OPT_SERIALIZE:
			mSerializeFile = args.OptionArg();
			break;
		}
	}
	
	// Parse the file arguments
	switch(args.FileCount())
	{
	case 0:
		// If method list requested, ignore files
		if(mListFingerprintingMethods  || mListDistanceMethods || mListGroupingMethods || mListAnalysisMethods) return;

		std::cerr << "Missing POSCAR file" << std::endl << std::endl;
		std::cerr << "Usage:" << std::endl;
		std::cerr << "    " << aVal[0] << usage_msg << std::endl << std::endl;
		showHelp(parser_options);
		throw CmdLineFatal();

	case 1:
		mPoscarFile = args.File(0);
		break;

	default:
		mPoscarFile = args.File(0);
		mEnergyFile = args.File(1);
		break;
	}

	// Few sanity checks on the parameters
	if(mEndStep > 0 && mEndStep < mStartStep)
	{
		throw CmdLineFatal("End step should be >= start step");
	}
	if(mHasEnergyThreshold && mHasReverseEnergyThreshold)
	{
		throw CmdLineFatal("Specify only one of energy threshold and reverse energy threshold");
	}
	if(mDiffrBinSize <= 0)
	{
		throw CmdLineFatal("Bin size cannot be <= 0");
	}
	if(mDiffrPeakSize < 0)
	{
		throw CmdLineFatal("Peak size cannot be < 0");
	}
	if(mMaxDistanceForGrouping < 0)
	{
		throw CmdLineFatal("Max distance for grouping cannot be < 0");
	}
}

