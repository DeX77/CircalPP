cmake_minimum_required (VERSION 2.6)
project (CircalPP)

add_definitions(-Wall -Wshadow -Werror -pedantic)

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
	
	# Find all the input files
	FILE(GLOB srcFiles "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
		
	list( REMOVE_ITEM srcFiles "${CMAKE_CURRENT_SOURCE_DIR}/src/CircalPP.cpp")
	
	CXXTEST_ADD_TEST(unittest_MatrixHelper MatrixHelperTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/MatrixHelperTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_RandomSequenceTest RandomSequenceTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/RandomSequenceTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_RotatedSequenceTest RotatedSequenceTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/RotatedSequenceTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_PseudoRotatedSequenceTest PseudoRotatedSequenceTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/PseudoRotatedSequenceTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_AlignmentTest AlignmentTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/AlignmentTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_CorrectedFastaTest CorrectedFastaTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/CorrectedFastaTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_WhitespaceFastaTest WhitespaceFastaTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/WhitespaceFastaTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_ScoringModelTest ScoringModelTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/ScoringModelTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_VertebrateMitochondrialGenomeAlphabetTest VertebrateMitochondrialGenomeAlphabetTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/VertebrateMitochondrialGenomeAlphabetTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_MultiplePseudoCircularAlignmentFactoryTest MultiplePseudoCircularAlignmentFactoryTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/MultiplePseudoCircularAlignmentFactoryTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_MultipleCircularAlignmentFactoryTest MultipleCircularAlignmentFactoryTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/MultipleCircularAlignmentFactoryTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_CircularAlignmentFactoryTest CircularAlignmentFactoryTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/CircularAlignmentFactoryTest.h ${srcFiles})
	
	CXXTEST_ADD_TEST(unittest_OutputTest OutputTest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/OutputTest.h ${srcFiles})
	
	
 	target_link_libraries(unittest_MatrixHelper ${LIBS_BPP} gcov)
	target_link_libraries(unittest_RandomSequenceTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_RotatedSequenceTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_PseudoRotatedSequenceTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_AlignmentTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_CorrectedFastaTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_WhitespaceFastaTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_ScoringModelTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_VertebrateMitochondrialGenomeAlphabetTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_MultiplePseudoCircularAlignmentFactoryTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_MultipleCircularAlignmentFactoryTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_CircularAlignmentFactoryTest ${LIBS_BPP} gcov)
	target_link_libraries(unittest_OutputTest ${LIBS_BPP} gcov)
	
	
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
