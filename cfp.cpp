/// @file cfp.cpp
/// Driver program for the CrystalFp library.
///
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <algorithm>
#include <vector>
#include "CmdLine.h"
#include "CrystalFp.h"
#include "ReadPoscar.h"
#include "AtomSymbols.h"
#include "CrystalFpAnalysis.h"
#include "CrystalFpScatterplot.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif

using namespace cfp;

/// Main program for CrystalFp.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-02-09 (initial version)
///     @version 1.0
///
///	@param[in] ac Number of command line parameters
/// @param[in] av Command line parameters
///
int main(int ac, char **av)
{
	try
	{

	// Parse the command line
	CmdLine cmd(ac, av, true);

	// Write out command line parameters
	if(cmd.mVerboseLevel > 0)
	{
		std::cerr << std::endl;
		std::cerr << "Poscar file:   " << cmd.mPoscarFile << std::endl;
		if(cmd.mEnergyFile) std::cerr << "Energy file:   " << cmd.mEnergyFile << std::endl;
		if(cmd.mVerboseLevel > 0) std::cerr << "Verbose level: " << cmd.mVerboseLevel << std::endl;
		if(!cmd.mAtomZ.empty())
		{
			std::cerr << "Elements:      ";
			std::vector<unsigned int>::const_iterator iaz;
			for(iaz=cmd.mAtomZ.begin(); iaz != cmd.mAtomZ.end(); ++iaz) std::cerr << elementZToSymbol(*iaz) << ' ';
			std::cerr << std::endl;
		}
		if(cmd.mStartStep > 1) 	std::cerr << "Start step:    " << cmd.mStartStep << std::endl;
		if(cmd.mEndStep > 0) 	std::cerr << "End step:      " << cmd.mEndStep << std::endl;
		if(cmd.mCutoffDistance > 0) std::cerr << "Cutoff at:     " << cmd.mCutoffDistance << std::endl;
		std::cerr << "Energy:        ";
		if(cmd.mEnergyIsPerAtom) std::cerr << "per atom" << std::endl;
		else                     std::cerr << "per structure" << std::endl;
		std::cerr << "Structures:    ";
		if(cmd.mIsNanocluster) std::cerr << "nanoclusters" << std::endl;
		else                   std::cerr << "crystals" << std::endl;
		std::cerr << "Bin size:      " << cmd.mDiffrBinSize << std::endl;
		std::cerr << "Peak size:     " << cmd.mDiffrPeakSize << std::endl;
		if(cmd.mForcedFpLen > 0) std::cerr << "Forced fp len: " << cmd.mForcedFpLen << std::endl;
		if(cmd.mCheckpointDir)   std::cerr << "Checkpt dir:   " << cmd.mCheckpointDir << std::endl;
		if(cmd.mOverwriteChkptDir && cmd.mCheckpointDir)   std::cerr << "Prev checkpt:  overwrite" << std::endl;
#ifdef _OPENMP
		std::cerr << "Max threads:   " << omp_get_max_threads() << std::endl;
#endif
		std::cerr << std::endl;
	}

	// Construct the classifier object
	CrystalFp cfp(cmd.mVerboseLevel);

	// List methods if requested
	if(cmd.mListFingerprintingMethods || cmd.mListDistanceMethods || cmd.mListGroupingMethods || cmd.mListAnalysisMethods)
	{
		if(cmd.mListFingerprintingMethods)
		{
			std::cerr << std::endl << "Available fingerprinting methods:" << std::endl;
			std::vector<std::string> fnames = cfp.getFingerprintMethodsNames();
			std::vector<std::string>::const_iterator ifn;
			int n=0;
			for(ifn=fnames.begin(); ifn != fnames.end(); ++ifn, ++n) std::cerr << "   " << n << " " << *ifn << std::endl;
		}

		if(cmd.mListDistanceMethods)
		{
			std::cerr << std::endl << "Available distance measures:" << std::endl;
			std::vector<std::string> fnames = cfp.getDistanceMethodsNames();
			std::vector<std::string>::const_iterator ifn;
			int n=0;
			for(ifn=fnames.begin(); ifn != fnames.end(); ++ifn, ++n) std::cerr << "   " << n << " " << *ifn << std::endl;
		}

		if(cmd.mListGroupingMethods)
		{
			std::cerr << std::endl << "Available grouping methods:" << std::endl;
			std::vector<std::string> fnames = cfp.getGroupingMethodsNames();
			std::vector<std::string>::const_iterator ifn;
			int n=0;
			for(ifn=fnames.begin(); ifn != fnames.end(); ++ifn, ++n) std::cerr << "   " << n << " " << *ifn << std::endl;
		}
		if(cmd.mListAnalysisMethods)
		{
			// Start the analysis method
			CrystalFpAnalysis analysis(&cfp);

			std::vector<std::string>::const_iterator ifn;
			int n=0;

			std::cerr << "\nSimple analysis methods (method_x,method_y)" << std::endl;
			const std::vector<std::string> n1 = analysis.getAnalysisMethodsNames(CrystalFpAnalysis::CATEGORY_SIMPLE);
			for(ifn=n1.begin(); ifn != n1.end(); ++ifn, ++n) std::cerr << "   " << n << " " << *ifn << std::endl;

			std::cerr << "\nHistogram analysis methods (method)" << std::endl;
			const std::vector<std::string> n2 = analysis.getAnalysisMethodsNames(CrystalFpAnalysis::CATEGORY_HIST);
			for(ifn=n2.begin(); ifn != n2.end(); ++ifn, ++n) std::cerr << "   " << n << " " << *ifn << std::endl;

			std::cerr << "\nSpecial analysis methods (method)" << std::endl;
			const std::vector<std::string> n3 = analysis.getAnalysisMethodsNames(CrystalFpAnalysis::CATEGORY_SPECIAL);
			for(ifn=n3.begin(); ifn != n3.end(); ++ifn, ++n) std::cerr << "   " << n << " " << *ifn << std::endl;

			std::cerr << "\nAll per structure analysis methods (method)" << std::endl;
			const std::vector<std::string> n4 = analysis.getAnalysisMethodsNames(CrystalFpAnalysis::CATEGORY_ALL);
			for(ifn=n4.begin(); ifn != n4.end(); ++ifn, ++n) std::cerr << "   " << n << " " << *ifn << std::endl;
		}

		std::cerr << std::endl;
		return 0;
	}

	// Check and set methods
	if(cmd.mFingerprintingMethod != CmdLine::NO_METHOD_SELECTED) cfp.setFingerprintMethod(cmd.mFingerprintingMethod);
	if(cmd.mDistanceMethod != CmdLine::NO_METHOD_SELECTED)       cfp.setDistanceMethod(cmd.mDistanceMethod);
	if(cmd.mGroupingMethod != CmdLine::NO_METHOD_SELECTED)       cfp.setGroupingMethod(cmd.mGroupingMethod);

	// Read requested steps
	if(cmd.mVerboseLevel >= 1) std::cerr << "Start file loading" << std::endl;
	readPoscarAndEnergies(cmd.mPoscarFile, cmd.mEnergyFile, cmd.mEnergyIsPerAtom, cmd.mStartStep, cmd.mEndStep, cmd.mAtomZ, cfp);

	// Select structures based on energies
	if(cfp.hasEnergies())
	{
		if(cmd.mVerboseLevel >= 2) std::cerr << "Min energy:       " << cfp.getMinEnergy() << std::endl;

		// Select based on energies
		if(cmd.mHasEnergyThreshold)
		{
			if(cmd.mVerboseLevel >= 2) std::cerr << "Threshold energy: " << cmd.mEnergyThreshold << std::endl;
			cfp.energyThreshold(cmd.mEnergyThreshold);
		}
		else if(cmd.mHasReverseEnergyThreshold)
		{
			cmd.mEnergyThreshold += cfp.getMinEnergy();
			if(cmd.mVerboseLevel >= 2) std::cerr << "Threshold energy: " << cmd.mEnergyThreshold << std::endl;
			cfp.energyThreshold(cmd.mEnergyThreshold);
		}
	}

	// Show the number of loaded structures
	if(cmd.mVerboseLevel >= 1) std::cerr << "Loaded " << cfp.getNumActiveStructures() << " of " << cfp.getNumTotalStructures() << std::endl;

	// Compute the cutoff distance to use and set it
	if(cmd.mCutoffDistance > 0)
	{
		if(cmd.mVerboseLevel >= 2) std::cerr << "Forced cutoff distance: " << std::fixed << std::setprecision(4) << cmd.mCutoffDistance << std::endl;
	}
	else
	{
		cmd.mCutoffDistance = cfp.computeCutoffDistance();
		if(cmd.mVerboseLevel >= 2) std::cerr << "Computed cutoff distance: " << std::fixed << std::setprecision(4) << cmd.mCutoffDistance << std::endl;
	}

	// Force the new dimensionality for the fingerprints (should be called before loadCheckpoint())
	cfp.forceFpLength(cmd.mForcedFpLen);

	// Set checkpointing
	cfp.setCheckpointDir(cmd.mCheckpointDir);

	// Load previous checkpoint if needed
	if(cmd.mCheckpointDir && !cmd.mOverwriteChkptDir)
	{
		if(cmd.mVerboseLevel >= 1) std::cerr << "Reload checkpoint" << std::endl;
		cfp.loadCheckpoint();
	}

	// Compute fingerprints
	if(cmd.mFingerprintingMethod != CmdLine::NO_METHOD_SELECTED)
	{
		if(cmd.mVerboseLevel >= 1)
		{
			std::cerr << std::endl << "Fingerprint computation using method: " << cfp.getFingerprintMethod();
			if(cfp.isDiffractionLike()) std::cerr << " (diffraction like)";
			std::cerr << std::endl;
		}

		// Set the kind of structure under analysis
		if(cmd.mIsNanocluster) cfp.setNanoclusterStructureType();

		// Set fingerprinting parameters
		cfp.setCutoffDistance(cmd.mCutoffDistance);
		cfp.setDiffrBinSize(cmd.mDiffrBinSize);
		cfp.setDiffrPeakSize(cmd.mDiffrPeakSize);

		// Compute fingerprints (not already loaded from the checkpoint dir)
		cfp.computeFingerprints();

		if(cmd.mVerboseLevel >= 1)
		{
			std::cerr << "Done" << std::endl;
			if(!cfp.hasFingerprints()) std::cerr << "Some fingerprints missing!" << std::endl;
			std::cerr << "Dimensionality: " << cfp.getFingerprintNumSections() * cfp.getFingerprintSectionLen() << std::endl;
		}
	}

	// Compute distances matrix
	if(cmd.mDistanceMethod != CmdLine::NO_METHOD_SELECTED)
	{
		if(cmd.mVerboseLevel >= 1) std::cerr << std::endl << "Distances computation using method: " << cfp.getDistanceMethod() << std::endl;

		cfp.computeDistanceMatrix();

		if(cmd.mVerboseLevel >= 1) std::cerr << "Done" << std::endl;
	}

	// Classify the results
	if(cmd.mGroupingMethod != CmdLine::NO_METHOD_SELECTED)
	{
		if(cmd.mVerboseLevel >= 1)
		{
			std::cerr << std::endl << "Start grouping using method: " << cfp.getGroupingMethod();
			if(cfp.groupingNeedsK()) std::cerr << " (needs K)";
			std::cerr << std::endl;
		}

		// Set parameters
    	cfp.setMaxGroupingDistance(cmd.mMaxDistanceForGrouping);
    	cfp.setK(cmd.mK);

		// Group the results
		cfp.groupResults();

		if(cmd.mVerboseLevel >= 1)
		{
			std::cerr << "Done" << std::endl;
			std::cerr << "Num grouped entries: " << cfp.getNgroups() << std::endl;
			std::cerr << "Num single entries:  " << cfp.getNsingle() << std::endl;
		}
	}

	// Remove duplicated structures
	if(cmd.mRemoveDuplicates && cfp.hasFingerprints() && cfp.hasDistanceMatrix() && cfp.getNgroups() > 0)
	{
		if(cmd.mVerboseLevel >= 1) std::cerr << std::endl << "Remove duplicates" << std::endl;

		std::vector<unsigned int> new_idx;
		int new_cnt = cfp.reduceDuplicatesToRepresentative(new_idx);

		// Save the map to a file if a filename has been specified and if mapping should be done
		if(cmd.mMapFile && new_cnt > 0)
		{
			std::ofstream out(cmd.mMapFile, std::ios_base::trunc | std::ios_base::out);
			if(!out.good())
			{
				std::cerr << "Cannot create duplicates map file <" << cmd.mMapFile << ">" << std::endl;
			}
			else
			{
				for(std::vector<unsigned int>::const_iterator ii= new_idx.begin(); ii != new_idx.end(); ++ii) out << *ii << std::endl;
				out.close();
			}
		}

		if(cmd.mVerboseLevel >= 1)
		{
			if(new_cnt > 0) std::cerr << "New count: " << new_cnt << std::endl;
			else            std::cerr << "Nothing done" << std::endl;
		}
	}

	// Output fingerprints in AVS FLD format
	if(cmd.mFldOutFp && cfp.hasFingerprints())
	{
		// Prepare the names of the fld and dat files
		std::string fld_out_fp = cmd.mFldOutFp;
		std::string ffld;
		std::string fdat;
		size_t ext = fld_out_fp.find_last_of('.');
		if(ext == std::string::npos)
		{
			ffld = fld_out_fp + ".fld";
			fdat = fld_out_fp + ".dat";
		}
		else
		{
			ffld = fld_out_fp.substr(0, ext) + ".fld";
			fdat = fld_out_fp.substr(0, ext) + ".dat";
		}
		size_t slash = fdat.find_last_of("/\\");
		std::string fdat_base = (slash != std::string::npos) ? fdat.substr(slash+1) : fdat;

		// Create the description file
		std::ofstream fld(ffld.c_str(), std::ios_base::trunc | std::ios_base::out);
		if(!fld.good())
		{
			std::cerr << "Cannot create fingerprints fld file <" << ffld << ">" << std::endl;
		}
		else
		{
			unsigned int fplen = cfp.getFingerprintNumSections() * cfp.getFingerprintSectionLen();

			fld << "# AVS field file" << std::endl;
			fld << "ndim=2        # number of dimensions in the field" << std::endl;
			fld << "dim1=" << std::left << std::setw(9) << fplen << "# dimension of axis 1 (fp length)" << std::endl;
			fld << "dim2=" << std::left << std::setw(9) << cfp.getNumActiveStructures() << "# dimension of axis 2 (num. structures)" << std::endl;
			fld << "nspace=2      # number of physical coordinates per point" << std::endl;
			fld << "veclen=1      # number of components at each point" << std::endl;
			fld << "data=float    # portable data format" << std::endl;
			fld << "field=uniform # field type (uniform, rectilinear, irregular)" << std::endl;
			fld << "variable 1 file=" << fdat_base << " filetype=binary skip=" << 3*sizeof(int) << std::endl;
			fld.close();
		}

		// Create the binary file (first two integers contain the dim1 and dim2 values)
		std::ofstream dat(fdat.c_str(), std::ios_base::binary | std::ios_base::trunc | std::ios_base::out);
		if(!dat.good())
		{
			std::cerr << "Cannot create fingerprints dat file <" << fdat << ">" << std::endl;
		}
		else
		{
			int fplen = (int)(cfp.getFingerprintNumSections() * cfp.getFingerprintSectionLen());
			dat.write((const char *)&fplen, sizeof(int));
			int numfp = (int)cfp.getNumActiveStructures();
			dat.write((const char *)&numfp, sizeof(int));
			int nparts = (int)cfp.getFingerprintNumSections();
			dat.write((const char *)&nparts, sizeof(int));

			for(unsigned int s=0; s < (unsigned int)numfp; ++s)
			{
				const float* f = cfp.getFingerprint(s);
				dat.write((const char *)f, sizeof(float)*fplen);
			}
			dat.close();
		}
	}

	// Output distances in AVS FLD format
	if(cmd.mFldOutDist && cfp.hasDistanceMatrix())
	{
		// Prepare the names of the fld and dat files
		std::string fld_out_dist = cmd.mFldOutDist;
		std::string ffld;
		std::string fdat;
		size_t ext = fld_out_dist.find_last_of('.');
		if(ext == std::string::npos)
		{
			ffld = fld_out_dist + ".fld";
			fdat = fld_out_dist + ".dat";
		}
		else
		{
			ffld = fld_out_dist.substr(0, ext) + ".fld";
			fdat = fld_out_dist.substr(0, ext) + ".dat";
		}
		size_t slash = fdat.find_last_of("/\\");
		std::string fdat_base = (slash != std::string::npos) ? fdat.substr(slash+1) : fdat;

		// Create the description file
		std::ofstream fld(ffld.c_str(), std::ios_base::trunc | std::ios_base::out);
		if(!fld.good())
		{
			std::cerr << "Cannot create distances fld file <" << ffld << ">" << std::endl;
		}
		else
		{
			unsigned int nstruct = cfp.getNumActiveStructures();

			fld << "# AVS field file" << std::endl;
			fld << "ndim=2        # number of dimensions in the field" << std::endl;
			fld << "dim1=" << std::left << std::setw(9) << nstruct << "# dimension of axis 1 (num. structures)" << std::endl;
			fld << "dim2=" << std::left << std::setw(9) << nstruct << "# dimension of axis 2 (num. structures)" << std::endl;
			fld << "nspace=2      # number of physical coordinates per point" << std::endl;
			fld << "veclen=1      # number of components at each point" << std::endl;
			fld << "data=float    # portable data format" << std::endl;
			fld << "field=uniform # field type (uniform, rectilinear, irregular)" << std::endl;
			fld << "variable 1 file=" << fdat_base << " filetype=binary skip=" << 2*sizeof(int) << std::endl;
			fld.close();
		}

		// Create the binary file (first two integers contain the dim1 and dim2 values)
		std::ofstream dat(fdat.c_str(), std::ios_base::binary | std::ios_base::trunc | std::ios_base::out);
		if(!dat.good())
		{
			std::cerr << "Cannot create distances dat file <" << fdat << ">" << std::endl;
		}
		else
		{
			int nstruct = (int)cfp.getNumActiveStructures();
			dat.write((const char *)&nstruct, sizeof(int));
			dat.write((const char *)&nstruct, sizeof(int));

			unsigned int s, v;
			unsigned int ns = cfp.getNumActiveStructures();
			for(s=0; s < ns; ++s)
			{
				for(v=0; v < ns; ++v)
				{
					float d = cfp.getDistance(s, v);

					dat.write((const char *)&d, sizeof(float));
				}
			}

			dat.close();
		}
	}

	// Output sorted distances
	if(cmd.mSortedDistFile && cfp.hasDistanceMatrix())
	{
		if(cmd.mVerboseLevel >= 1) std::cerr << std::endl << "Start sorting distances" << std::endl;
		std::ofstream outsorted(cmd.mSortedDistFile, std::ios_base::binary | std::ios_base::trunc | std::ios_base::out);
		if(!outsorted.good())
		{
			std::cerr << "Cannot create sorted distances file <" << cmd.mSortedDistFile << ">" << std::endl;
		}
		else
		{
			unsigned int ns = cfp.getNumActiveStructures();
			unsigned int s, v;
			std::vector<float> dist;
			dist.reserve((ns*(ns-1))/2);
			for(s=0; s < ns-1; ++s)
			{
				for(v=s+1; v < ns; ++v)
				{
					float d = cfp.getDistance(s, v);
					dist.push_back(d);
				}
			}
			std::sort(dist.begin(), dist.end());

			int ndist = (int)dist.size();
			outsorted.write((const char *)&ndist, sizeof(int));

			std::vector<float>::const_iterator idist;
			for(idist = dist.begin(); idist != dist.end(); ++idist)
			{
				float d = *idist;
				outsorted.write((const char *)&d, sizeof(float));
			}
			outsorted.close();
		}
		if(cmd.mVerboseLevel >= 1) std::cerr << "End sorting distances" << std::endl;
	}

	// Output a summary of the computed quantities
	if(cmd.mSummaryFile)
	{
		std::ofstream outsummary(cmd.mSummaryFile, std::ios_base::trunc | std::ios_base::out);
		if(!outsummary.good())
		{
			std::cerr << "Cannot create summary <" << cmd.mSummaryFile << ">" << std::endl;
		}
		else
		{
			outsummary << "Loaded:" << cfp.getNumActiveStructures() << std::endl;
			outsummary << "Computed cutoff distance:" << std::fixed << std::setprecision(4) << cfp.computeCutoffDistance() << std::endl;
			outsummary << "Actual cutoff distance:" << std::fixed << std::setprecision(4) << cfp.getCutoffDistance() << std::endl;
			outsummary << "Dimension:" << cfp.getFingerprintNumSections() * cfp.getFingerprintSectionLen() << std::endl;

			outsummary.close();
		}
	}

	// Compute the required analysis
	if(cmd.mAnalysisMethod != CmdLine::NO_METHOD_SELECTED)
	{
		// Start the analysis method
		CrystalFpAnalysis analysis(&cfp);

		// Set the analysis (and verify valididy)
		if(!analysis.setAnalysisMethod(cmd.mAnalysisMethod, cmd.mAnalysisMethod2)) throw CrystalFpFatal("Invalid analysis requested");
#if 0
		// Enable caching (useful only for testing here because makes sense only in interactive usage)
		analysis.enableCaching();
#endif
		// Set parameters
		for(std::map<std::string,std::string>::const_iterator im=cmd.mAnalysisParams.begin(); im != cmd.mAnalysisParams.end(); ++im)
		{
			analysis.setNamedParam(im->first, im->second);
		}
		
		// Compute number of values (first element is the number of X values following)
		std::vector<unsigned int> num_results = analysis.numValues();

		// Output the values on file (well does not make much sense not to have this set)
		if(cmd.mAnalysisFile)
		{
			// Create the analysis results file
			std::ofstream fa(cmd.mAnalysisFile, std::ios_base::trunc | std::ios_base::out);
			if(!fa.good())
			{
				std::cerr << "Cannot create analysis file <" << cmd.mAnalysisFile << ">" << std::endl;
			}
			else
			{
				// Output header
				std::vector<std::string> lbl = analysis.getLabels();
				std::vector<std::string>::const_iterator ilbl;
				bool first = true;
				for(ilbl=lbl.begin(); ilbl != lbl.end(); ++ilbl)
				{
					if(first) first = false; else fa << ',';
					fa << '"' << *ilbl << '"';
				}
				fa << std::endl;

				// Output values
				unsigned int i, j;
				float **out = new float*[num_results.size()-1];
				for(i=0; i < num_results.size()-1; ++i) out[i] = new float[num_results[i+1]];
				analysis.getValues(out);

				for(i=0; i < num_results[1]; ++i)
				{
					first = true;
					for(j=0; j < num_results.size()-1; ++j)
					{
						if(first) first = false; else fa << ',';
						fa << std::scientific << std::setprecision(6) << out[j][i];
					}
					fa << std::endl;
				}

				// Close file and release memory
				fa.close();
				for(i=0; i < num_results.size()-1; ++i) delete [] out[i];
				delete [] out;
			}
		}
	}

	// Creaste the corresponding scatterplot
	if(cmd.mCreateScatterplot)
	{
		if(cmd.mVerboseLevel >= 1) std::cerr << std::endl << "Start creating scatterplot" << std::endl;

		// Create the scatterplot object and set few default values
		CrystalFpScatterplot sp;
		unsigned int num_retries = 1;
		float		 min_energy = 1e-6F;
		unsigned int max_iterations = 600;
		CrystalFpScatterplot::ValueType value_type = CrystalFpScatterplot::VAL_STRESS;
		CrystalFpScatterplot::DiagnosticType diagnostic_type = CrystalFpScatterplot::DIAG_BINNED_DISTANCES;
		float timestep = 0.02F;

		// Set parameters from command line
		for(std::map<std::string,std::string>::const_iterator im=cmd.mScatterplotParams.begin(); im != cmd.mScatterplotParams.end(); ++im)
		{
			if(!strncasecmp(im->first.c_str(), "retry", 1))
			{
				num_retries = atoi(im->second.c_str());

				sp.setNamedParam(im->first, im->second);
			}
			else if(!strncasecmp(im->first.c_str(), "energy", 1))
			{
				min_energy = (float)atof(im->second.c_str());
			}
			else if(!strncasecmp(im->first.c_str(), "iterations", 1))
			{
				max_iterations = atoi(im->second.c_str());
			}
			else if(!strncasecmp(im->first.c_str(), "kind", 1))
			{
				value_type = (CrystalFpScatterplot::ValueType)atoi(im->second.c_str());
			}
			else if(!strncasecmp(im->first.c_str(), "diagnostic", 1))
			{
				diagnostic_type = (CrystalFpScatterplot::DiagnosticType)atoi(im->second.c_str());
			}
			else if(!strncasecmp(im->first.c_str(), "timestep", 1))
			{
				timestep = (float)atof(im->second.c_str());
			}
		}

		// Recap the scatterplot values
		if(cmd.mVerboseLevel >= 1)
		{
			sp.dumpParams();
			std::cerr << "Num. retries:   " << std::setw(12) << num_retries << std::endl;
			std::cerr << "Min energy:     " << std::setw(12) << min_energy << std::endl;
			std::cerr << "Max iterations: " << std::setw(12) << max_iterations << std::endl;
			std::cerr << "Value type:     " << std::setw(12) << value_type << std::endl;
			std::cerr << "Diagnostic:     " << std::setw(12) << diagnostic_type << std::endl;
			std::cerr << "Timestep:       " << std::setw(12) << timestep << std::endl;
		}

		// Initialize the scatterplot and run it the required number of retries
		unsigned int npoints = sp.initScatterplot(&cfp);
		unsigned int i, j;
		for(i=0; i < num_retries; ++i)
		{
			if(cmd.mVerboseLevel >= 1) std::cerr << std::endl << "Start retry " << i << std::endl;

			// Iterate and stop on energy or num. iterations criterias
			for(j=0; j < max_iterations; ++j)
			{
				float energy = sp.stepScatterplot(timestep);
				if(energy < min_energy) break;
			}

			// Perturb the point position to run another retry
			if(num_retries > 1) sp.perturbPositions();
		}

		// Output the values on file (well, does not make much sense not to have this set)
		if(cmd.mScatterplotFile)
		{
			// Create the analysis results file
			std::ofstream fa(cmd.mScatterplotFile, std::ios_base::trunc | std::ios_base::out);
			if(!fa.good())
			{
				std::cerr << "Cannot create scatterplot file <" << cmd.mScatterplotFile << ">" << std::endl;
			}
			else
			{
				float* coords  = new float[2*npoints];
				float* vals    = new float[npoints];

				// Get the values
				sp.getPoints(coords);
				sp.getValues(vals, value_type);

				// Output header
				fa << "x,y,value" << std::endl;

				// Output values
				for(i=0; i < npoints; ++i)
				{
					fa << std::scientific << std::setprecision(6) << coords[2*i+0] << ','
					   << std::scientific << std::setprecision(6) << coords[2*i+1] << ','
					   << std::scientific << std::setprecision(6) << vals[i] << std::endl;
				}

				// Close file and release memory
				fa.close();
				delete [] coords;
				delete [] vals;
			}
		}
		if(cmd.mDiagnosticFile && diagnostic_type != CrystalFpScatterplot::DIAG_DO_NOTHING)
		{
			// Create the analysis results file
			std::ofstream fa(cmd.mDiagnosticFile, std::ios_base::trunc | std::ios_base::out);
			if(!fa.good())
			{
				std::cerr << "Cannot create scatterplot diagnostic file <" << cmd.mDiagnosticFile << ">" << std::endl;
			}
			else
			{
				npoints = sp.initDiagnostic(diagnostic_type);
				float* coords  = new float[2*npoints];
				float* vals    = new float[npoints];

				// Get the values
				sp.getDiagnosticValues(coords, vals);

				// Output header
				fa << "x,y,value" << std::endl;

				// Output values (for binned, only the not empty bins)
				for(i=0; i < npoints; ++i)
				{
					if(diagnostic_type != CrystalFpScatterplot::DIAG_BINNED_DISTANCES || vals[i] > 0)
					{
						fa << std::scientific << std::setprecision(6) << coords[2*i+0] << ','
						   << std::scientific << std::setprecision(6) << coords[2*i+1] << ','
						   << std::scientific << std::setprecision(6) << vals[i] << std::endl;
					}
				}

				// Close file and release memory
				fa.close();
				delete [] coords;
				delete [] vals;
			}
		}
	}

	// Serialize the CrystalFp class to file
	if(cmd.mSerializeFile)
	{
		if(cmd.mVerboseLevel >= 1) std::cerr << std::endl << "Start serializing content" << std::endl;
		std::ofstream serialized(cmd.mSerializeFile, std::ios_base::binary | std::ios_base::trunc | std::ios_base::out);
		if(!serialized.good())
		{
			std::cerr << "Cannot create serialized file <" << cmd.mSerializeFile << ">" << std::endl;
		}
		else
		{
			cfp.serialize(serialized);
			serialized.close();
		}
		if(cmd.mVerboseLevel >= 1) std::cerr << "End serializing content" << std::endl;

		// TEST
		if(cmd.mVerboseLevel >= 1) std::cerr << "Try to deserialize content" << std::endl;

		CrystalFp xcfp(cmd.mVerboseLevel);
		std::ifstream unserialized(cmd.mSerializeFile, std::ios_base::binary | std::ios_base::in);
		if(!unserialized.good())
		{
			std::cerr << "Cannot open serialized file <" << cmd.mSerializeFile << ">" << std::endl;
		}
		else
		{
			xcfp.unserialize(unserialized);
			unserialized.close();
		}
		if(cmd.mVerboseLevel >= 1) std::cerr << "End deserializing content" << std::endl;

		std::cerr << std::endl << "**** Original" << std::endl;
		cfp.dump();

		std::cerr << std::endl << "**** Reloaded" << std::endl;
		xcfp.dump();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// Catch all exceptions
	}
	catch(CmdLineSuccess&)
	{
		return 0;
	}
	catch(CrystalFpFatal& e)
	{
		// Don't print throw messages for which the message has been already printed
		const char* p = e.what();
		if(*p) std::cerr << std::endl << "Error: " << p << std::endl;
		return 1;
	}
	catch(CmdLineFatal& e)
	{
		// Don't print throw messages for which the message has been already printed
		const char* p = e.what();
		if(*p) std::cerr << std::endl << "Error: " << p << std::endl;
		return 1;
	}
	catch(...)
	{
		std::cerr << std::endl << "Error: Default exception" << std::endl;
		return 1;
	}
}

