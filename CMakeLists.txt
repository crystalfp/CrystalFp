cmake_minimum_required(VERSION 2.6)
project(CrystalFp)

# Source code
set(DRIVER_SRC
	cfp.cpp
	CmdLine.cpp
	AtomSymbols.cpp
	ReadPoscar.cpp
)

set(LIB_SRC
	CrystalFp.cpp
	Structure.cpp
	StructureList.cpp
	FingerprintingMethod.cpp
	SmoothPeak.cpp
	DistanceMethod.cpp
	GroupingMethod.cpp
	CrystalFpAnalysis.cpp
	AnalysisMethod.cpp
	ResultsCache.cpp
	CrystalFpScatterplot.cpp
)

find_package(OpenMP)

if(OPENMP_FOUND)
	SET(USE_OPENMP YES CACHE BOOL "Use OpenMP")
endif(OPENMP_FOUND)

# Set OpenMP and other platform specific flags
if(WIN32)
	add_definitions(/D_SECURE_SCL=0)
	add_definitions(/D_CRT_SECURE_NO_WARNINGS)

	if(USE_OPENMP AND OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif(USE_OPENMP AND OPENMP_FOUND)

else(WIN32)

	if(USE_OPENMP AND OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif(USE_OPENMP AND OPENMP_FOUND)

endif(WIN32)

# Create library
add_library(CrystalFp ${LIB_SRC})

# Executable
add_executable(cfp ${DRIVER_SRC})

# Link the executable to the CrystalFp library
target_link_libraries(cfp CrystalFp) 

# Build the documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")
   set(DOXY_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")
  endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")

  add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG})
  
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_dev")
   set(DOXY_DEV_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_dev")
  endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_dev")

  add_custom_target(devdoc ${DOXYGEN_EXECUTABLE} ${DOXY_DEV_CONFIG})
  
endif(DOXYGEN_FOUND)

