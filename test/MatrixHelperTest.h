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
#include "Alignment.h"

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
  Circal::Alignment* pairwiseAlign;

  MatrixHelperTest() :
      scoreM(new Circal::ScoringModel()), alpha(new bpp::DNA()), A("A", "AAAA",
          alpha), B("B", "BBBB", alpha), C(B), pairwiseAlign(
          new Circal::Alignment(alpha))
    {
      pairwiseAlign->addSequence(A);
      pairwiseAlign->addSequence(B);
    }
  ;
  ~MatrixHelperTest()
    {
      delete pairwiseAlign;
      delete alpha;
      delete scoreM;

    }
  ;

  void testMatrixHelperDestructor(void)
    {
      TS_TRACE("Starting MatrixHelperDestructor test");

      Circal::MatrixHelper* temp = new Circal::MatrixHelper();
      delete temp;

      TS_TRACE("Finishing MatrixHelperDestructor test");
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

      TS_ASSERT_EQUALS(output.size(), A.size()+1);
      TS_ASSERT_EQUALS(output.at(0).size(), B.size()+1);

      TS_ASSERT_EQUALS(output.at(0).at(0), 0.5);
      TS_ASSERT_EQUALS(output.at(A.size()).at(B.size()), 0.5);
      TS_TRACE("Finishing InitScoreMatrixWith test");
    }
  ;

  void testInitScoreMatrix3DWith(void)
    {
      TS_TRACE("Starting InitScoreMatrix3DWith test");
      Circal::ScoreMatrix3D empty;
      int delta = 2;
      Circal::ScoreMatrix3D output = mhelper.InitScoreMatrix3DWith(A, C, delta,
          0.5);

      TS_ASSERT_EQUALS(output.size(), A.size()+1);
      TS_ASSERT_EQUALS(output.at(0).size(), C.size()+1);
      TS_ASSERT_EQUALS(output.at(0).at(0).size(), (C.size()/2)/delta);

      TS_ASSERT_EQUALS(output.at(0).at(0).at(0), 0.5);
      TS_ASSERT_EQUALS(output.at(A.size()).at(C.size()).at((C.size()/2/delta)-1), 0.5);

      TS_ASSERT_DIFFERS( output, empty);

      TS_TRACE("Finishing InitScoreMatrix3DWith test");
    }
  ;

  void testCreateAdjacenceGraph(void)
    {
      TS_TRACE("Starting CreateAdjacenceGraph test");

      Circal::BoolMatrix output = mhelper.CreateAdjacenceGraph(pairwiseAlign,
          2);

      TS_ASSERT(output.size() > 0);
      TS_TRACE("Finishing CreateAdjacenceGraph test");
    }
  ;

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

  void testInvalidCutRowFromTo(void)
    {
      TS_TRACE("Starting CutRowFromTo with invalid data test");

      Circal::ScoreMatrix sm = mhelper.InitScoreMatrixWith(A, B, 0.5);
      int sizebefore = sm.size();

      mhelper.CutRowFromTo(&sm, 0, 100);
      TS_ASSERT_EQUALS(sm.size(), sizebefore);

      TS_TRACE("Finishing CutRowFromTo  with invalid data  test");
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

  void testInvalidCutColumnFromTo(void)
    {
      TS_TRACE("Starting CutColumnFromTo with invalid range test");

      Circal::ScoreMatrix sm = mhelper.InitScoreMatrixWith(A, B, 0.5);
      int sizebefore = sm.at(0).size();

      mhelper.CutColumnFromTo(&sm, 0, 100);
      TS_ASSERT_EQUALS(sm.at(0).size(), sizebefore);

      TS_TRACE("Finishing CutColumnFromTo with invalid range  test");
    }
  ;

  void testSearchBestPositionFrom(void)
    {
      TS_TRACE("Starting BestPositionFrom test");

      double empty = -10;

      Circal::ScoreMatrix sm = mhelper.InitScoreMatrixWith(A, B, 0.5);
      uint temp = 1;

      double output = mhelper.SearchBestPositionFrom(&sm, temp, temp, scoreM);

      TS_ASSERT_DIFFERS( output, empty);

      TS_TRACE("Finishing BestPositionFrom test");

    }
  ;

  void testSearchBestInRow(void)
    {
      TS_TRACE("Starting SearchBestInRow test");

      Circal::ScoreMatrix sm = mhelper.InitScoreMatrixWith(A, B, 0.5);
      uint temp1 = 1;
      uint temp0 = 0;

      int output = mhelper.SearchBestInRow(&sm, scoreM, temp0, temp1);
      TS_ASSERT_DIFFERS( output, 1);

      TS_TRACE("Finishing SearchBestInRow test");
    }
  ;
  void testSearchBestInRow3D(void)
    {
      TS_TRACE("Starting SearchBestInRow3D test");

      Circal::ScoreMatrix3D sm3d = mhelper.InitScoreMatrix3DWith(A, B, 0, 0.8);
      uint temp1 = 1;
      uint temp0 = 0;
      uint temp2 = 2;
      int output = mhelper.SearchBestInRow3D(&sm3d, scoreM, temp0, temp1,
          temp2);
      TS_ASSERT_DIFFERS( output, 9999);

      TS_TRACE("Finishing SearchBestInRow3D test");
    }
  ;

  void testSearchBestInColumn(void)
    {
      TS_TRACE("Starting SearchBestInColumn test");

      Circal::ScoreMatrix sm = mhelper.InitScoreMatrixWith(A, B, 0.5);
      uint temp1 = 1;
      uint temp0 = 0;

      int output = mhelper.SearchBestInColumn(&sm, scoreM, temp0, temp1);
      TS_ASSERT_DIFFERS( output, 999);

      TS_TRACE("Finishing SearchBestInColumn test");
    }
  ;

  void testSearchBestInColumn3D(void)
    {
      TS_TRACE("Starting SearchBestInColumn3D test");

      Circal::ScoreMatrix3D sm3d = mhelper.InitScoreMatrix3DWith(A, B, 0, 0.8);
      uint temp1 = 1;
      uint temp0 = 0;
      uint temp2 = 2;
      int output = mhelper.SearchBestInColumn3D(&sm3d, scoreM, temp0, temp1,
          temp2);
      TS_ASSERT_DIFFERS( output, 9999);

      TS_TRACE("Finishing SearchBestInColumn3D test");
    }
  ;

  };

#endif /* MATRIXHELPER_TEST_H_ */
