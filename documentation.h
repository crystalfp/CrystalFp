/// @file documentation.h
/// This file contains only the CrystalFp library documentation to be processed by Doxygen.
///
/// @mainpage CrystalFp: the library for crystal fingerprinting and classification
///
///
/// @section intro_sec CrystalFp: the crystal fingerprinting project
///
/// %CrystalFp started as a way to solve a specific problem with the <a href="http://mysbfiles.stonybrook.edu/~aoganov/USPEX.html">USPEX</a>
/// crystal structure predictor.
/// Every USPEX run produces hundreds or thousands of crystal structures, some of which may be identical.
/// To ease the extraction of unique and potentially interesting structures a method to find and remove duplicated structures has to be found.
///
/// The approach adopted was to apply usual high-dimensional classification concepts to the unusual field of crystallography.
///
/// We adopted a visual design and validation method to develop a classifier library (CrystalFp) and an end-user application
/// to select and validate method choices, to gain users' acceptance and to tap into their domain expertise.
///
/// Using the end-user application with real datasets, we experimented with various crystal structure descriptors,
/// distinct distance measures and tried different clustering methods to identify groups of similar structures.
/// These methods are already applied in combinatorial chemistry to organic molecules for a different goal and
/// in somewhat different forms, but are not widely used for crystal structures classification.
///
/// The use of the classifier has already accelerated the analysis of USPEX output by at least one order of magnitude,
/// promoting some new crystallographic insight and discovery. Furthermore the visual display of key algorithm indicators has led to diverse,
/// unexpected discoveries that will improve the USPEX algorithms.
///
///
/// @section unexp_sec Unexpected discoveries
///
/// Looking at the data we found unexpected correlations and patterns that deserve further investigation. So the research continues.
///
/// 
/// @section res_sec Resources
///
/// Up to date information on the CrystalFp library can be found on the <a href="http://www.cscs.ch/~mvalle/CrystalFp/">project pages</a>.
///
/// Here, the fingerprinting method section collects general information about the project, its philosophy and related publications.
/// The library contains the source code and API documentation of the CrystalFp library.
/// And, finally, the current research section acts as a (restricted) repository of ideas, things to do and so on.
///
/// On the wiki there are pages about the tools we are using to analyze the data with a particular emphasis on HPC methods to leverage CSCS computational
/// resources and speedup the analysis work.
///
///
/// @section contacts_sec Contacts
///
/// Contact us if you want more information on the project, want to collaborate or have new ideas.
///
///    - Ing. <a href="mailto:mvalle@cscs.ch">Mario Valle</a> - Swiss National Supercomputing Centre (CSCS) - Switzerland
///    - Prof. <a href="mailto:aoganov@notes.cc.sunysb.edu">Artem R. Oganov</a> - Dept. of Geosciences and New York Center for Computational Science, State University of New York at Stony Brook. 
///


/// @page cppstd_page C++ Coding Standard
/// Here are collected few rules for coding this project.
///
/// @section cnames_sect Class names
/// Class names are CamelCase with first letter uppercase.
///
/// Ex: CrystalFp
///
/// @section cmemb_sect Class methods
/// Class methods names are CamelCase with the first letter lowercase.
///
/// Ex: addStructureBatch
///
/// There is one exception: when the function mimics some well know one like size(), clean(), serialize().
///
/// @section cmemb_sect Class data members
/// Class member variables names start with 'm' followed by CamelCase name.
///
/// Ex: mDistanceMatrix
///
/// @section carg_sect Function arguments
/// Function arguments names start with 'a' followed by CamelCase name.
///
/// Ex: aEnergyThreshold
///
/// @section const_sect Constants and enumeration
/// Constants and enumerations are all uppercase with words separated by '_'.
/// The first letters specify the kind of constant (like: STS_ status, OPT_ option value).
///
/// Ex: OPT_GROUPING_METHOD
///
/// @section stack_sect Temporary variables
/// All the other variables are all lower case with parts separated by '_'.
///
/// Ex: max_basis_len
///
/// @section cglob_sect Global data
/// Global variables should be avoided, but if present their names should start with 'g' followed by CamelCase name.
///
/// Ex: gFgBranch
///
/// @section misc_sect Miscellaneous rules
/// In case of error main should return 1.
///
/// Counters should be unsigned int. But beware of down counting for loops.
///


/// @page cmd_page Command Line switches
/// The driver program can set various CrystalFp parameters from the command line and produces (if needed) answers as csv files.
/// At a minimum the program requires a POSCAR file containing the structures. To use energy related functionality, a file containing energies
/// (one floating point value per line in the same order as the structures in the POSCAR file) is required.
///
/// To list the available numerical options for \c --fingerprint-method \c --distance-method and \c --grouping-method simply follow the switch by \c help.
///
/// To list \c --scatterplot-param and \c --analysis-param param names simply enter \c h \c h after the switch.
///
/// @code
/// Usage:
///     cfp [options] POSCARfile [ENERGIESfile]
/// 
/// -v  --verbose (optional argument)
///         Verbose level (if no argument, defaults to 1)
/// 
/// -?  -h  --help (no argument)
///         This help
/// 
/// -t  --elements (required argument)
///         List of chemical elements to be assigned to POSCAR atoms
/// 
/// -es  --max-step  --end-step (required argument)
///         Last step to load (default: all)
/// 
/// -ss  --start-step (required argument)
///         First step to load (default: first)
/// 
/// -et  --energy-per-structure (no argument)
///         Energy from file is per structure, not per atom
/// 
/// -e  --energy-threshold (required argument)
///         Energy threshold
/// 
/// -r  --threshold-from-min (required argument)
///         Threshold from minimum energy
/// 
/// -c  --cutoff-distance (required argument)
///         Fingerprint forced cutoff distance
/// 
/// -n  --nano-clusters  --nanoclusters (no argument)
///         The structures are nanoclusters, not crystals
/// 
/// -b  --bin-size (required argument)
///         Bin size for the pseudo-diffraction methods
/// 
/// -p  --peak-size (required argument)
///         Peak smearing size
/// 
/// -x  --force-dim (required argument)
///         Force to this value the total length of the fingerprint (valid only on reloading from checkpoint)
/// 
/// -chk  --checkpoint-dir (required argument)
///         Set fingerprint checkpoint directory
/// 
/// -i  --ignore-checkpoint (no argument)
///         Ignore previous checkpoint
/// 
/// -f  --fingerprint-method (required argument)
///         Compute fingerprints using the given method ('help' to list all methods)
/// 
/// -d  --distance-method (required argument)
///         Compute distances using the given method ('help' to list all methods)
/// 
/// -g  --grouping-method (required argument)
///         Grouping method ('help' to list all methods)
/// 
/// -gt  --grouping-threshold (required argument)
///         Grouping max distance threshold
/// 
/// -k  --k (required argument)
///         K value needed by some grouping algorithm
/// 
/// -s  --summary (required argument)
///         Output a summary file to be used by another program
/// 
/// -fo  --fld-fingerprints (required argument)
///         Output fingerprints in FLD AVS format
/// 
/// -do  --fld-distances (required argument)
///         Output distance matrix in FLD AVS format
/// 
/// -sd  --sorted-distances (required argument)
///         Output distances in increasing order
/// 
/// -rd  --remove-dupl (optional argument)
///         Remove duplicates substituting them with a representative structure and writing the map of indices to file
/// 
/// -a  --analysis (required argument)
///         Create a chart for the given analysis method (can be: ax,ay or ax)
/// 
/// -af  --analysis-file (required argument)
///         File output for analysis
/// 
/// -ap  --analysis-param (multiple arguments)
///         Pass one parameter for the given analysis method as: code value (code: bins, part)
/// 
/// -sc  --scatterplot (no argument)
///         Create a scatterplot
/// 
/// -sf  --scatterplot-file (required argument)
///         File output for scatter plot
/// 
/// -df  --diagnostic-file (required argument)
///         File output for diagnostic chart assocated to the scatterplot
///
/// -sp  --scatterplot-param (multiple arguments)
///         Pass one parameter for the scatterplot as: code value (code: retry, energy, iterations, kind, diagnostic, timestep, mass, stiffness, damping, perturb, bins, wobble)
/// 
/// @endcode


/// @page ex_page Usage example
///
/// @code
/// cfp \
/// 	--verbose=3 \
/// 	--elements="Si O" \
/// 	--cutoff-distance=30 \
/// 	--fingerprint-method=0 \
/// 	--distance-method=0 \
/// 	--grouping-method=1 \
/// 	--grouping-threshold=0.000002 \
/// 	--remove-dupl=map.dat \
/// 	--checkpoint-dir=./chk \
/// 	--summary=summary.dat \
/// 	--fld-fingerprints=fp.fld  \
/// 	--fld-distances=dist.fld \
/// 	--sorted-distances=sorted.dat \
/// 	gatheredPOSCARS \
/// 	enthalpies_nospace.dat
/// @endcode
/// 
/// Verbose: show also the progress of the various computations (on standard error).
/// Compute the fingerprints with a cutoff distance of 30 using the "Per element diffraction" method (number 0).
/// Use distance measure "Cosine distance" (number 0) and grouping method "Hierarchical Single Linkage" (number 1) with a threshold of 2e-6.
/// Based on the grouping, remove duplicates and put the correspondence [new index-old index] in the <tt>map.dat</tt> file.
/// All the fingerprint details are saved in the <tt>./chk</tt> directory that will be reloaded in case of crash or
/// to compute lower dimensionality values (use <tt>--force-dim=newDim</tt>).
/// The \c gatheredPOSCARS file contains the structures and the \c enthalpies_nospace.dat file the corresponding energies.
/// 
/// The above example produces the following files (documentation on the \ref fmt_page page):
/// \li \c summary.dat with computation summary data
/// \li \c fp.fld and \c fp.dat with the fingerprints (the fld file is useful to load the result inside AVS/Express)
/// \li \c dist.fld and \c dist.dat with the distancew matrix (the fld file is useful to load the result inside AVS/Express)
/// \li \c sorted.dat with the sorted distances
///


/// @page fmt_page Output file formats
///
/// @section fmt_fp_sect Binary output for fingerprints
///
/// fplen (int) Length of each fingerprint\n
/// numfp (int) Number of structures (i.e. the number of fingerprints)\n
/// nparts (int) Number of parts of each fingerprint\n
/// fp1 (fplen*float) First fingerprint\n
/// fp2 (fplen*float) Next fingerprint\n
/// ...\n
/// fp(numfp) (fplen*float) Last fingerprint\n
///
///
/// @section fmt_dist_sect Binary output for distance matrix
///
/// nstruct (int) Number of structures\n
/// nstruct (int) Number of structures\n
/// distances (nstruct*nstruct*float) distance matrix by row\n
///
///
/// @section fmt_sort_sect Binary output for sorted distances
///
/// ndist (int) Number of distances\n
/// distances (ndist*float) ascending sorted distances\n
///
///
/// @section fmt_chkpt_sect Checkpoint file
///
/// Filename encodes the StepId:\n
///		FPP%010d.dat (positive StepId)
///		FPN%010d.dat (negative StepId)
/// The StepId absolute value goes in the ten digits zero padded field.
///
///	step (int) StepId value\n
/// fingerprint_total_length (unsigned int) Total length of the fingerprint\n
/// fingerprint_parts (unsigned int) Number of parts\n
/// fingerprint (float*fingerprint_total_length) The fingerprint proper\n
///
///


/// @page bld_page Build
///
/// The library uses CMAKE as its build system.
/// There is a switch that can be set during configuration to enable OpenMP usage to speedup fingerprint and distances computation.
/// As usual the kind of build (Release or Debug) should be selected on Linux. On Windows the build types are available in the generated
/// Visual C++ project.
///
/// CMAKE generates other two targets:
/// "doc" generates the library interface documentation. Instead "devdoc" generates documentation for the developers for all sources.
///


/// @page future_page Future and to do
///
/// Create the serialize/unserialize mechanism for CrystalFp class.
///
/// Enable reading of VASP poscar + elements
/// 
/// Put more formats readers
///