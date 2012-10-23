/*
 Circal++  - Calculates Multiple conservative Alignments of Circular Sequences
 Copyright (C) 2007  Daniel Exner
 <dex@dragonslave.de>

 This program is free software; you can redistribute it
 and/or modify it under the terms of the GNU General Public License as published
 by the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
 St, Fifth Floor, Boston, MA 02110, USA
 */

#ifndef MATRIXHELPER_TEST_H_
#define MATRIXHELPER_TEST_H_

#include "MatrixHelper.h"
#include "VertebrateMitochondrialGenomeAlphabet.h"
#include "ScoringModel.h"
#include "SequenceProxy.h"
#include "PseudoRotatedSequence.h"

#include <Bpp/Seq/Alphabet/DNA.h>
#include <cxxtest/TestSuite.h>

class MatrixHelperTest: public CxxTest::TestSuite
  {
public:

  Circal::MatrixHelper mhelper;
  Circal::ScoringModel* scoreM;
  bpp::Alphabet* alpha;
  Circal::SequenceProxy A, B;
  Circal::PseudoRotatedSequence C;

  MatrixHelperTest() :
      scoreM(new Circal::ScoringModel()), alpha(new bpp::DNA()), A("A", "AA",
          alpha), B("B", "BB", alpha), C(B)
    {
    }
  ;
  ~MatrixHelperTest()
    {
      delete scoreM;
      delete alpha;

    }
  ;

  void testInitializeScoreMatrixDistances(void)
    {
      TS_TRACE("Starting InitializeScoreMatrixDistances test");
      Circal::ScoreMatrix empty;
      Circal::ScoreMatrix output = mhelper.InitializeScoreMatrixDistances(A, B,
          scoreM);
      TS_ASSERT_DIFFERS( output, empty);
      TS_ASSERT_EQUALS(output.at(0).at(0), 0);
      TS_TRACE("Finishing InitializeScoreMatrixDistances test");
    }
  ;

  void testInitScoreMatrixWith(void)
    {
      TS_TRACE("Starting InitScoreMatrixWith test");
      Circal::ScoreMatrix empty;
      Circal::ScoreMatrix output = mhelper.InitScoreMatrixWith(A, B, 0.5);
      TS_ASSERT_DIFFERS( output, empty);

      TS_ASSERT_EQUALS(output.at(0).at(0), 0.5);
      TS_TRACE("Finishing InitScoreMatrixWith test");
    }
  ;

  void testInitScoreMatrix3DWith(void)
    {
      TS_TRACE("Starting InitScoreMatrix3DWith test");
      Circal::ScoreMatrix3D empty;
      Circal::ScoreMatrix3D output = mhelper.InitScoreMatrix3DWith(A, C, 2,
          0.5);
      TS_ASSERT_EQUALS(output.at(0).at(0).at(0), 0.5);

      TS_ASSERT_DIFFERS( output, empty);

      TS_TRACE("Finishing InitScoreMatrix3DWith test");
    }
  ;
//
//      void testCreateAdjacenceGraph(void)
//        {
//          BoolMatrix CreateAdjacenceGraph(Alignment* pairWiseAlignments,
//              int biggestSequenceSize);
//        }
//      ;
//
  void testCutRowFromTo(void)
    {
      TS_TRACE("Starting CutRowFromTo test");

      Circal::ScoreMatrix sm = mhelper.InitScoreMatrixWith(A, B, 0.5);
      int sizebefore = sm.size();

      mhelper.CutRowFromTo(&sm, 0, 1);
      TS_ASSERT_EQUALS(sm.size(), sizebefore-1);

      TS_TRACE("Finishing CutRowFromTo test");
    }
  ;

  void testCutColumnFromTo(void)
    {
      TS_TRACE("Starting CutColumnFromTo test");

      Circal::ScoreMatrix sm = mhelper.InitScoreMatrixWith(A, B, 0.5);
      int sizebefore = sm.at(0).size();

      mhelper.CutColumnFromTo(&sm, 0, 1);
      TS_ASSERT_EQUALS(sm.at(0).size(), sizebefore-1);

      TS_TRACE("Finishing CutColumnFromTo test");
    }
  ;
//
//      void testSearchBestPositionFrom(void)
//        {
//          double SearchBestPositionFrom(const ScoreMatrix* M, uint &i, uint &j,
//              const ScoringModel* scoreM);
//        }
//      ;
//
//      void testSearchBestInRow(void)
//        {
//          int SearchBestInRow(const ScoreMatrix* M, const ScoringModel* scoreM,
//              const uint &start, const uint &row);
//        }
//      ;
//      void testSearchBestInRow3D(void)
//        {
//          int SearchBestInRow3D(const ScoreMatrix3D* M,
//              const ScoringModel* scoreM, const uint &start, const uint &row,
//              uint &k);
//        }
//      ;
//
//      void testSearchBestInColumn(void)
//        {
//          int SearchBestInColumn(const ScoreMatrix* M,
//              const ScoringModel* scoreM, const uint &start,
//              const uint &column);
//        }
//      ;
//
//      void testSearchBestInColumn3D(void)
//        {
//          int SearchBestInColumn3D(const ScoreMatrix3D* M,
//              const ScoringModel* scoreM, const uint &start, const uint &column,
//              uint &k);
//        }
//      ;

  };

#endif /* MATRIXHELPER_TEST_H_ */
