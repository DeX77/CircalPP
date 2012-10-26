add_executable(CircalPP PseudoCircularAlignmentFactory.cpp
						MultipleAlignmentFactory.cpp
						MatrixHelper.cpp
						AlignmentSymbol.cpp
						AlignmentFactory.cpp
						RotatedSequence.cpp
						CircularAlignmentFactory.cpp
						MultiplePseudoCircularAlignmentFactory.cpp
						CircalPP.cpp
						Alignment.cpp
						ScoringModel.cpp
						PseudoRotatedSequence.cpp
						ErrorClasses.cpp
						Output.cpp
						VertebrateMitochondrialGenomeAlphabet.cpp
						RandomSequence.cpp
						WhitespaceFasta.cpp
						CorrectedFasta.cpp
						SequenceProxy.cpp
)
target_link_libraries(${PROJECT_NAME} ${LIBS_BPP})