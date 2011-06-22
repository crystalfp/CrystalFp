/// @file CrystalFpExceptions.h
/// The exceptions raised by CrystalFp library.
///
#ifndef CRYSTALFPEXCEPTIONS_H
#define CRYSTALFPEXCEPTIONS_H

#include <stdexcept>

namespace cfp
{

/// Fatal exception in CrystalFp.
///
/// The message explain the reason. If absent then print nothing in the exception handler
/// (useful if a more complex message has already been printed)
///
class CrystalFpFatal : public std::runtime_error
{
public:
	/// Create the exception.
	///
	/// @param aMsg An optional message string.
	///
	CrystalFpFatal(const char *aMsg="") : runtime_error(aMsg)
	{}
};

}
#endif

