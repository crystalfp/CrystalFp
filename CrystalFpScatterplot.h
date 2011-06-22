/// @file CrystalFpScatterplot.h
/// Computes the 2D projection of the multidimensional CrystalFp object.
///
#ifndef CRYSTALFPSCATTERPLOT_H
#define CRYSTALFPSCATTERPLOT_H

#include "CrystalFp.h"
namespace cfp
{

/// Computes the 2D projection of the multidimensional CrystalFp object.
/// The points are projected in 2D trying to preserve the original distances.
///	Normal usage is:
/// @code
///		CrystalFpScatterplot sp;
///		npoints = sp.initScatterplot(&cfp);
///		allocate npoints*2 float for cooordinates (coords) and  for values (vals)
///		sp.getPoints(coords);
///		sp.getValues(vals, CrystalFpScatterplot::VAL_TOTAL_ENERGY);
///
///		for(num. retries) {
///			for(num. iterations) {
///				energy = sp.stepScatterplot(timestep);
///				sp.getPoints(coords);
///				sp.getValues(vals, CrystalFpScatterplot::VAL_TOTAL_ENERGY);
///				if(energy < min_energy) break;
///			}
///			sp.perturbPositions();
///			sp.getPoints(coords);
///			sp.getValues(vals, CrystalFpScatterplot::VAL_TOTAL_ENERGY);
///		}
///		sp.getPoints(coords);
///		sp.getValues(vals, CrystalFpScatterplot::VAL_TOTAL_ENERGY);
///
/// @endcode
///		
///	@author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///	@date 2011-06-01 (initial version)
///
///	@version 1.0
///
class CrystalFpScatterplot
{
public:
	/// Constructor.
	///
	CrystalFpScatterplot();

	/// Destructor.
	///
	~CrystalFpScatterplot();
	

	/// Set the parameter identified by the name.
	/// Mainly used to set parameters from the command line.
	/// The accepted codes are:
	///	\li \c retry		Number of retries (used only here to compute the decay coefficient)
	///	\li \c mass			Ball mass
	///	\li \c stiffness	Of the spring
	///	\li \c damping		Damping factor for the movement
	///	\li \c perturb		Perturb scale for the retries (it perturb the initial position of the masses)
	///	\li \c bins			Number of bins for the binned distances diagnostics
	///	\li \c wobble		Wobble scale for the position of the binned points
	///
	/// @param aName Name of the parameter to be set (see above)
	/// @param aValue String representation of the value to be set (internally it is converted to the needed data type)
	///
	/// @exception CrystalFpFatal If the param name is invalid
	///
	void setNamedParam(const std::string& aName, const std::string& aValue);

	/// Initialize the scatterplot computation for the given CrystalFp object.
	/// The method used is the ball and spring simple method.
	///
	/// @param[in] aCfp Pointer to the already computed CrystalFp structure.
	///
	/// @exception CrystalFpFatal If CrystalFp has no distances available.
	///
	unsigned int initScatterplot(const CrystalFp *aCfp);

	/// Type of value associated to the points.
	///
	enum ValueType
	{
		VAL_TOTAL_ENERGY,		///< The total energy associated to the structure
		VAL_PER_ATOM_ENERGY,	///< Per atom energy associated to the structure
		VAL_STRESS,				///< The stress computed by the multidimensional scaling algorithm
		VAL_GROUP,				///< Group to which the point pertains
		VAL_STEP				///< The step value
	};

    /// Get the points coordinates.
    /// If no points, then do nothing.
    ///
    /// @param[out] aCoords Array of x,y coordinates filled on output.
    ///
    void getPoints(float* aCoords) const;

    /// Return the values associated to the points.
	/// It is an independent routine so you can change the visualized value without redisplaying the points.
    ///
    /// @param[out] aValues An array where the values will be stored.
    /// @param[in] aValueType The value to be returned.
	///
	/// @exception CrystalFpFatal If aValueType is invalid.
    ///
    void getValues(float* aValues, ValueType aValueType) const;
	
	/// Move the points one step ahead.
	///
	/// @param[in] aTimestep The time from the previous step.
	///
	///	@return The cinetic energy of the points (to decide when they are sufficiently "calm")
	///
	float stepScatterplot(float aTimestep);

	/// Check if the current solution is the best.
	/// If the current solution is not an improvement, then reload the previous optimum solution.
	/// Perturb the point positions to escape a local minima.
	/// The perturbation is proportional to the point stress and is slowly reduced every restart.
	///
	void perturbPositions(void);

	/// Type of diagnostic requested.
	///
	enum DiagnosticType
	{
		DIAG_DISTANCES,			///< Return all inter-point distances
		DIAG_BINNED_DISTANCES,	///< Bin the interpoint distances
		DIAG_DO_NOTHING			///< To simplify visualization of the scatterplot only
	};

	/// Compute scatterplot diagnostic chart.
	/// The chart shows distances on the scatterplot vs. distances in real space (both distances are normalized between 0 and 1).
	///
	/// @param[in] aDiagnostic The diagnostic to be computed
	///
	/// @return The number of points
	///
	/// @exception CrystalFpFatal If aDiagnostic is invalid.
	///
	unsigned int initDiagnostic(DiagnosticType aDiagnostic);

	/// Return the diagnostic chart points
	///
	/// @param[out] aCoords A preallocated array that will be filled by the point coordinates (x, y)
	/// @param[out] aValues A preallocated array that will be filled with the point's values
	///
	void getDiagnosticValues(float* aCoords, float* aValues) const;

	/// Prints the params values (the one set with setNamedParam()).
	/// The print is on cerr
	///
	void dumpParams(void) const;

private:
	void colorByGroup(unsigned int aNumPoints, float *aResult) const;
	float computeCost(void) const;

private:
	const CrystalFp* mCrystalFp;

private:
	struct CrystalFpScatterplotImpl;
	struct CrystalFpScatterplotImpl* mPimpl;
};

}
#endif

