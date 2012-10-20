cmake_minimum_required (VERSION 2.6)
project (CircalPP)

set(CMAKE_CXX_FLAGS "-Wall -Wshadow -Weffc++")

IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

FIND_PACKAGE( OpenMP REQUIRED )

if(OPENMP_FOUND)
 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


# --- the IMPROVED_FIND_LIBRARY macro (taken from the Bio++ CMakeLists.txt) ---
macro(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
  FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})

  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} found here:")
    MESSAGE("    ${${lib_name}_INCLUDE_DIR}/")
    MESSAGE("    ${${lib_name}_LIBRARY}")
  ELSE(${lib_name}_LIBRARY)
    MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
  ENDIF(${lib_name}_LIBRARY)
  
  #add the dependency:
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)

#--- Libraries & includes ---
INCLUDE_DIRECTORIES(src)
IMPROVED_FIND_LIBRARY(LIBS_BPP bpp-seq "Bpp/Seq/Sequence.h")
IMPROVED_FIND_LIBRARY(LIBS_BPP bpp-core "Bpp/Clonable.h")

add_executable(CircalPP src/PseudoCircularAlignmentFactory.cpp
						src/MultipleAlignmentFactory.cpp
						src/MatrixHelper.cpp
						src/AlignmentSymbol.cpp
						src/GenomeAlphabet.cpp
						src/MultipleCircularAlignmentFactory.cpp
						src/AlignmentFactory.cpp
						src/RotatedSequence.cpp
						src/CircularAlignmentFactory.cpp
						src/MultiplePseudoCircularAlignmentFactory.cpp
						src/CircalPP.cpp
						src/Alignment.cpp
						src/ScoringModel.cpp
						src/PseudoRotatedSequence.cpp
						src/ErrorClasses.cpp
						src/Output.cpp
						src/VertebrateMitochondrialGenomeAlphabet.cpp
						src/RandomSequence.cpp
						src/WhitespaceFasta.cpp
						src/CorrectedFasta.cpp
						src/SequenceProxy.cpp
)
target_link_libraries(${PROJECT_NAME} ${LIBS_BPP})
