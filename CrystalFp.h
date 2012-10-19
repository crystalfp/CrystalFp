/// @file CrystalFp.h
/// Main interface to CrystalFp library.
///
///
#ifndef CRYSTALFP_H
#define CRYSTALFP_H

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif
#include <vector>
#include <set>
#include <string>
#include <climits>
#include "CrystalFpExceptions.h"

/// @namespace cfp
/// Public namespace for the library
///
namespace cfp
{
/// CrystalFp public interface.
///
///	@author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///	@date 2008-06-06 (initial version)
///	@date 2009-09-01 (version 1.0)
///	@version 1.0
///
class CrystalFp
{
public:
	//###################################################################################
	/// @name Constructor, destructor and structure reset
	//@{

	/// Constructor
	/// @param aVerboseLevel Set the verbosity level of the library
	/// - 0: No messages
	/// - 1: Moderate
	/// - 2: Something more
	/// - 3: Trace all the structure operations (which step is processed)
	///
	CrystalFp(unsigned int aVerboseLevel);

	/// Destructor
	///
	~CrystalFp();

	/// Reset everything to the status just after the construction
	///
	void resetAll(void);

	//@}


	//###################################################################################
	/// @name Available methods names
	/// Names of the fingerprinting, distance and classification methods
	//@{

	/// Return the names of the implemented fingerprint computation methods
	/// @return The vector of fingerprinting methods names
	///
	const std::vector<std::string> getFingerprintMethodsNames(void) const;

	/// Return the names of the implemented distance computation methods
	/// @return The list of distance measure methods names
	///
	const std::vector<std::string> getDistanceMethodsNames(void) const;

	/// Return the names of the implemented grouping methods
	/// @return The list of grouping method names
	///
	const std::vector<std::string> getGroupingMethodsNames(void) const;

	//@}


	//###################################################################################
	/// @name Load phase methods
	/// Load structures to be processed
	//@{

	/// Load one new crystal structure optionally with an energy
	///
	/// @param aStep A numeric identifier of the structure
	/// @param aNumAtoms Number of atoms of the structure to be loaded
	/// @param aCoords The coordinates of the atoms. The array contains contiguous x1, y1, z1, x2, y2, z2, ...
	/// @param aZ The Z value of the corresponding atom
	/// @param aUnitCell The unit cell given as a vector[16]: [0..2] a; [4..6] b;  [8..10] c; [12..14] origin; [15] set to 1.0 and the rest set to 0.0
	/// @param aHasEnergy If true, the structure has associated energy
	/// @param aEnergy The internal energy or enthalpy of the structure
	/// @param aEnergyIsPerAtom If true the energy is per atom and not per structure
	///
	void  			addStructure(int aStep, unsigned int aNumAtoms, const float *aCoords, const unsigned int *aZ, const float *aUnitCell, bool aHasEnergy=false, float aEnergy=0.0F, bool aEnergyIsPerAtom=true);

	/// Load one new crystal structure optionally with an energy, but it does not update the selected structures list till addStructureBatchFinish() is called.
	///
	/// After a set of AddStructureBatch() calls have been made, the AddStructureBatchFinish() call should be issued.
	/// @param aStep A numeric identifier of the structure
	/// @param aNumAtoms Number of atoms of the structure to be loaded
	/// @param aCoords The coordinates of the atoms. The array contains contiguous x1, y1, z1, x2, y2, z2, ...
	/// @param aZ The Z value of the corresponding atom
	/// @param aUnitCell The unit cell given as a vector[16]: [0..2] a; [4..6] b;  [8..10] c; [12..14] origin; [15]  set to 1.0 and the rest set to 0.0
	/// @param aHasEnergy If true, the structure has associated energy
	/// @param aEnergy The internal energy or enthalpy of the structure
	/// @param aEnergyIsPerAtom If true the energy is per atom and not per structure
	///
	void  			addStructureBatch(int aStep, unsigned int aNumAtoms, const float *aCoords, const unsigned int *aZ, const float *aUnitCell, bool aHasEnergy=false, float aEnergy=0.0F, bool aEnergyIsPerAtom=true);

	/// Finish a batch loading
	///
	void			addStructureBatchFinish(void);

	/// Return the number of structures loaded that will be used for the computational phases.
	///
	/// @return The number of active structures
	///
	size_t			getNumActiveStructures(void) const;

	/// Return the total number of structures loaded.
	///
	/// @return The number of loaded structures (active and ignored)
	///
	size_t			getNumTotalStructures(void) const;

	//@}


	//###################################################################################
	/// @name Select active structures methods
	/// Select structures based on energies
	//@{

	/// Mark as active only the structures with energy less than the given value.
	///
	/// @param aEnergyThreshold The threshold energy value
	///
	void			energyThreshold(float aEnergyThreshold);

	/// Mark as active all the loaded structures regardless of their energy.
	///
	void			noEnergyThreshold(void);

	//@}

	
	//###################################################################################
	/// @name Methods to access global quantities
	/// Return quantities related to the full set of selected structures
	//@{

	/// Test if all the loaded structures have an associated energy value.
	///
	/// @return True if all the loaded structures have an associated energy value
	///
	bool			hasEnergies(void) const;

	/// Get the minimum energy value between all the loaded structures.
	///
	/// @return The minimum energy value
	///
	float			getMinEnergy(void) const;
	
	/// Check if the structures have unit cell.
	///
	/// @return True if all structures have unit cell
	///
	bool			hasUnitCell(void) const;

	//@}


	//###################################################################################
	/// @name Methods to access structure quantities at given index
	/// Return quantities related to the given structure
	//@{

	/// Get the step number for the selected structure.
	///
	/// @param aIdx Index of the structure
	/// @return The step number
	///
	int				idxToStep(size_t aIdx) const;

	/// Get the total energy value for the selected structure.
	///
	/// @param aIdx Index of the structure
	/// @return The energy value
	///
	float			getTotalEnergy(size_t aIdx) const;

	/// Get the energy per atom value for the selected structure.
	///
	/// @param aIdx Index of the structure
	/// @return The energy value
	///
	float			getPerAtomEnergy(size_t aIdx) const;
	
	/// Get the number of atoms for the selected structure.
	///
	/// @param aIdx Index of the structure
	/// @return The number of atoms
	///
	unsigned int	getNatoms(size_t aIdx) const;

	/// Get the unit cell for the selected structure.
	///
	/// @param aIdx Index of the structure
	/// @return The unit cell
	///
	const float*	getUnitCell(size_t aIdx) const;
	
	/// Get the kind of atoms for the selected structure
	/// @param aIdx Index of the structure
	/// @return The array of atom's types
	///
	const unsigned int* getAtomZ(size_t aIdx) const;

	/// Get the coordinates for the selected structure.
	///
	/// @param aIdx Index of the structure
	/// @return The coordinates array
	///
	const float*	getCoords(size_t aIdx) const;

	//@}


	//###################################################################################
	/// @name Compute fingerprints
	/// Compute fingerprints for all the active structures
	//@{

	/// Return the suggested cutoff distance.
	///
	/// @param aMargin The percentage to increase the cutoff distance to avoid border cases (default 2%)
	/// @return The suggested cutoff distance
	///
	float			computeCutoffDistance(float aMargin=0.02F) const;

	/// Set the fingerprinting method to be used.
	///
	/// @param aFingerprintType The kind of fingerprint to compute (index in the list of names)
	///
	/// @exception CrystalFpFatal If aFingerprintType is invalid
	///
	void   			setFingerprintMethod(unsigned int aFingerprintType);

	/// Check if the method is "per element diffraction" and similar.
	///
	/// @return True if the selected method is diffraction like
	///
	bool			isDiffractionLike(void) const;

	/// Get the description of the currently used fingerprinting method.
	///
	/// @return The description
	///
	const std::string getFingerprintMethod(void) const;

	/// Set the cutoff distance to be used.
	///
	/// @param aCutoff The cutoff distance to use
	///
	void			setCutoffDistance(float aCutoff);

	/// Modify the method for nanoclusters.
	/// If this function is called, the unit cell is ignored
	///
	void			setNanoclusterStructureType(void);

	/// Set diffraction bin width
	/// @param aBinSize The diffraction bin size
	///
	void			setDiffrBinSize(float aBinSize);

	/// Set the gaussian peak smoothing width.
	///
	/// @param aPeakSize The gaussian peak width to be used (if zero the peak is not smoothed)
	///
	void			setDiffrPeakSize(float aPeakSize);

	/// Force the fingerprint dimensionality to the given value.
	///
	/// @param aDim The forced dimensionality value
	///
	void			forceFpLength(unsigned int aDim);

	/// Set the checkpoint directory.
	/// In this directory one file is created for each fingerprint.
	///
	/// @param aDir The directory path
	///
	void			setCheckpointDir(const char* aDir);

	/// Load the data contained in the checkpoint directory.
	///
	void			loadCheckpoint(void);

	/// Compute the fingerprints.
	///
	void   			computeFingerprints(void);

	/// Get the number of sections composing the fingerprint
	/// @return The number of sections
	///
	unsigned int   	getFingerprintNumSections(void) const;

	/// The fingerprint section length.
	///
	/// @return The length of one fingerprint section
	///
	unsigned int	getFingerprintSectionLen(void) const;

	/// Check if the fingerprints for the loaded structures have been computed.
	///
	/// @return True if all fingerprints have been computed (or loaded from checkpoint dir)
	///
	bool  			hasFingerprints(void) const;

	/// Get the computed or set cutoff distance.
	///
	/// @return The cutoff distance or zero if not set
	///
	float			getCutoffDistance(void) const;

	/// Access the fingerprint for the selected structure.
	///
	/// @param aStructureIdx  Index of the structure
	/// @return The fingerprint (length is num_sections * section_length)
	///
	const float*	getFingerprint(size_t aStructureIdx) const;

	/// Get the gaussian peak smoothing width.
	///
	/// @return The gaussian peak width used (could be zero)
	///
	float			getDiffrPeakSize(void) const;

	/// Get the bin width for diffraction like fingerprints smoothing.
	///
	/// @return The diffraction bin size
	///
	float			getDiffrBinSize(void) const;

	/// Return the interatomic distances from the given atom in the given structure.
	///
	/// @param aIdx Index of the structure
	/// @return Vector (one for each atom) of vector of distances to all other atoms in the extended unit cell
	///
	const std::vector< std::vector<float> >& getInteratomicDistances(unsigned int aIdx) const;
	
	/// Get the weights for the selected structure.
	///
	/// @param aIdx Index of the structure
	/// @return The weights (one for each fingerprint part)
	///
	const float*	getWeights(size_t aIdx) const;

	//@}

	
	//###################################################################################
	/// @name Compute distances
	/// Compute distances between fingerprints
	//@{

	/// Set the distance measure to be used.
	///
	/// @param aMeasureType The kind of distance measure (index in the list of names)
	///
	/// @exception CrystalFpFatal If aMeasureType is invalid
	///
	void   			setDistanceMethod(unsigned int aMeasureType);

	/// Get the description of the currently used distance measure.
	///
	/// @return The description
	///
	const std::string getDistanceMethod(void) const;

	/// Compute distances between all fingerprints.
	///
	/// @exception CrystalFpFatal If fingerprints has been not computed before.
	///
	void  			computeDistanceMatrix(void);

	/// Check if the distances for the loaded structures have been computed.
	///
	/// @return True if the distance matrix has been computed
	///
	bool  			hasDistanceMatrix(void) const;

	/// Returns the distance between two structures.
	///
	/// @param aIdx1 Index of the first structure
	/// @param aIdx2 Index of the second structure
	/// @return The distance between the two structures
	///
	/// @exception CrystalFpFatal On invalid indices.
	///
	float 			getDistance(size_t aIdx1, size_t aIdx2) const;

	/// Get the maximum distance between two fingerprints.
	///
	/// @return The maximum distance
	///
	float			getMaxDistance(void) const;

	//@}


	//###################################################################################
	/// @name Classify structures methods
	/// Classify structures into groups
	//@{

	/// Set the grouping method to be used.
	///
	/// @param aGroupingMethod The kind of grouping method (index in the list of names)
	///
	/// @exception CrystalFpFatal If aGroupingMethod is invalid
	///
	void   			setGroupingMethod(unsigned int aGroupingMethod);

	/// Get the description of the currently used grouping method.
	///
	/// @return The description
	///
	const std::string getGroupingMethod(void) const;

	/// Do the grouping.
	///
	/// @exception CrystalFpFatal If distances has not been computed before.
	///
	void			groupResults(void);

	/// Set the distance threshold for grouping.
	///
	/// @param aMaxDistance The threshold distance
	///
	void			setMaxGroupingDistance(float aMaxDistance);

	/// Get the distance threshold for grouping.
	///
	/// @return The threshold distance
	///
	float			getMaxGroupingDistance(void) const;

	/// Set the common neighbors count.
	/// It is used by some grouping methods
	///
	/// @param aK The common neighbors count
	///
	void			setK(unsigned int aK);

	/// Return the number of groups found.
	///
	/// @return The number of groups found
	///
	unsigned int   	getNgroups(void) const;

	/// Return the list of groups found.
	///
	/// @return Array of groups each represented as array of the index of the structures in this group
	///
	const std::vector< std::set<unsigned int> >& getGroups(void) const;

	/// Return the number of ungrouped entries found.
	///
	/// @return The number of ungrouped entries (called single) found
	///
	unsigned int   	getNsingle(void) const;

	/// Test if the grouping method selected needs K parameter.
	///
	/// @return True if K needed
	///
	bool			groupingNeedsK(void) const;

	/// Remove grouped structures leaving only one representative structure.
	///
	/// @param[out] aNewIndexList The routine fills the given array with the original indices of the new, reduced, set of structures
	/// @return The new number of active structures or zero if nothing has been done
	///
	unsigned int	reduceDuplicatesToRepresentative(std::vector<unsigned int>& aNewIndexList);

	//@}


	//###################################################################################
	/// @name Serialize/unserialize class
	///  Serialize/unserialize class
	//@{

	/// Serialize the class to the given binary stream
	///
	/// @param[in] aStream The stream to which the class will be serialized
	///
	/// @exception CrystalFpFatal On write error.
	///
	void serialize(std::ofstream& aStream) const;

	/// Unserialize the class from the given binary stream
	///
	/// @param[in] aStream The stream from which the class should be deserialized
	/// @param[in] aAppend Append the read structures to the list of already read structures
	/// @param[in] aStepOffset Apply this offset to the id of the appended structures
	///
	/// @exception CrystalFpFatal On read error and on format validation error.
	///
	void unserialize(std::ifstream& aStream, bool aAppend=false, int aStepOffset=10000);

	//@}

	//###################################################################################
	/// @name Debugging support
	///  Class debugging support
	//@{

	void dump(void) const;
	//@}

private:
	CrystalFp(const CrystalFp&);
	CrystalFp& operator=(const CrystalFp&);

private:
	static void unitCellInverse(const float *aUnitCell, float aUnitCellInverse[3][3]);
	void computeExpansion(const float* aUnitCell, unsigned int* aExpansion) const;

private:
	struct CrystalFpImpl;
	struct CrystalFpImpl* mPimpl;
};

}

#endif

