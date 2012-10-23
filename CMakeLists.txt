cmake_minimum_required (VERSION 2.6)
project (CircalPP)

add_definitions(-Wall -Weffc++ -Wshadow)
#-Werror -pedantic 

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
INCLUDE_DIRECTORIES(include)
IMPROVED_FIND_LIBRARY(LIBS_BPP bpp-seq "Bpp/Seq/Sequence.h")
IMPROVED_FIND_LIBRARY(LIBS_BPP bpp-core "Bpp/Clonable.h")

add_subdirectory (src)

find_package(CxxTest)
if(CXXTEST_FOUND)
	include_directories(${CXXTEST_INCLUDE_DIR})
	enable_testing()
	add_definitions(-O0 -ggdb -fprofile-arcs -ftest-coverage)
	 	
	CXXTEST_ADD_TEST(unittest_MatrixHelper MatrixHelperTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/MatrixHelperTest.h
	${CMAKE_CURRENT_SOURCE_DIR}/src/PseudoCircularAlignmentFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MultipleAlignmentFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MatrixHelper.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/AlignmentSymbol.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/GenomeAlphabet.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MultipleCircularAlignmentFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/AlignmentFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RotatedSequence.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/CircularAlignmentFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MultiplePseudoCircularAlignmentFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Alignment.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/ScoringModel.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/PseudoRotatedSequence.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/ErrorClasses.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Output.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/VertebrateMitochondrialGenomeAlphabet.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/RandomSequence.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/WhitespaceFasta.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/CorrectedFasta.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/SequenceProxy.cpp)
		
	target_link_libraries(unittest_MatrixHelper ${LIBS_BPP} gcov)
endif()

find_program(LCOV_TOOL "lcov")
find_program(LCOV_HTML "genhtml")

if(LCOV_TOOL AND LCOV_HTML AND CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")

    if(NOT DEFINED "LCOV_OUTPUT_DIR")
        set(LCOV_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doc)
    endif()

    if(NOT DEFINED "LCOV_SCAN_DIR")
        set(LCOV_SCAN_DIR ${CMAKE_CURRENT_BINARY_DIR})
    endif()

	add_custom_target(lcov)
    add_dependencies(lcov check)

    add_custom_command(TARGET lcov POST_BUILD
        COMMAND ${LCOV_TOOL} --capture --directory ${LCOV_SCAN_DIR} --output-file ${CMAKE_CURRENT_BINARY_DIR}/coverage.lcov
		COMMAND ${LCOV_TOOL} --remove ${CMAKE_CURRENT_BINARY_DIR}/coverage.lcov '/usr/include/*' --output-file ${CMAKE_CURRENT_BINARY_DIR}/coverage.lcov
		COMMAND ${LCOV_TOOL} --remove ${CMAKE_CURRENT_BINARY_DIR}/coverage.lcov '/usr/local/include/*' --output-file ${CMAKE_CURRENT_BINARY_DIR}/coverage.lcov
		COMMAND ${LCOV_TOOL} --remove ${CMAKE_CURRENT_BINARY_DIR}/coverage.lcov '*Test.*' --output-file ${CMAKE_CURRENT_BINARY_DIR}/coverage.lcov
        COMMAND ${LCOV_HTML} ${CMAKE_CURRENT_BINARY_DIR}/coverage.lcov --no-branch-coverage --output-directory ${LCOV_OUTPUT_DIR}/lcov
    )
endif()
