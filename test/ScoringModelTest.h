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

#ifndef SCORINGMODELTEST_H_
#define SCORINGMODELTEST_H_

#include "ScoringModel.h"

#include <cxxtest/TestSuite.h>

class ScoringModelTest: public CxxTest::TestSuite
  {
public:

  std::string sm;

  ScoringModelTest() :
      sm("p 10  1 : A B C D E F G H I J K L M N O P Q R S T U V W X Y Z\n")
    {
    }
  ;

  void testConstructor(void)
    {
      TS_TRACE("Starting constructor test");
      Circal::ScoringModel test;

      TS_TRACE("Finishing constructor test");
    }
  ;

  void testConstructorIstream(void)
    {
      TS_TRACE("Starting constructor with istream test");
      std::istringstream input(sm);

      Circal::ScoringModel test(input);

      TS_TRACE("Finishing constructor with istream  test");
    }
  ;

  void testDestructor(void)
    {
      TS_TRACE("Starting constructor test");
      Circal::ScoringModel* test = new Circal::ScoringModel();

      delete test;

      TS_TRACE("Finishing constructor test");
    }
  ;

  void testGetModelSize(void)
    {
      TS_TRACE("Starting GetModelSize test");
      std::istringstream input(sm);
      Circal::ScoringModel test(input);

      TS_ASSERT_EQUALS(test.GetModelSize(), 52);

      TS_TRACE("Finishing GetModelSize test");

    }
  ;

  void testScoreOf(void)
    {
      TS_TRACE("Starting ScoreOf test");
      std::istringstream input(sm);

      Circal::ScoringModel test(input);

      TS_ASSERT_EQUALS(test.ScoreOf("A", "A"), test.ScoreOf("B", "B"));

      TS_TRACE("Finishing ScoreOf test");
    }
  ;

  void testScoreOfGapOpen(void)
    {
      TS_TRACE("Starting ScoreOfGapOpen test");
      std::istringstream input(sm);

      Circal::ScoringModel test(input);

      TS_ASSERT_EQUALS(test.ScoreOfGapOpen("A"), test.ScoreOfGapOpen("B"));

      TS_TRACE("Finishing ScoreOfGapOpen test");
    }
  ;

  void testScoreOfGapExtend(void)
    {
      TS_TRACE("Starting ScoreOfGapExtend test");
      std::istringstream input(sm);

      Circal::ScoringModel test(input);

      TS_ASSERT_EQUALS(test.ScoreOfGapExtend("A"), test.ScoreOfGapExtend("B"));

      TS_TRACE("Finishing ScoreOfGapExtend test");
    }
  ;

  void testBestOfTwo(void)
    {
      TS_TRACE("Starting BestOfTwo test");
      std::istringstream input(sm);
      Circal::ScoringModel test(input);

      double a = 0;
      double b = 1;

      TS_ASSERT_EQUALS(test.BestOfTwo(a,b), b);

      TS_TRACE("Finishing BestOfTwo test");
    }
  ;

  void testBestOfThree(void)
    {
      TS_TRACE("Starting BestOfThree test");
      std::istringstream input(sm);
      Circal::ScoringModel test(input);

      double a = 0;
      double b = 1;
      double c = 99;

      TS_ASSERT_EQUALS(test.BestOfThree(a,b,c), c);

      TS_TRACE("Finishing BestOfThree test");
    }
  ;

  void testNormalizeScores(void)
    {
      TS_TRACE("Starting NormalizeScores test");
      std::istringstream input(sm);

      Circal::ScoringModel test(input);

      double before = test.ScoreOf("A", "B");

      test.NormalizeScores(before);

      TS_ASSERT_EQUALS(test.ScoreOf("A", "B"), 0);

      TS_TRACE("Finishing NormalizeScores test");
    }
  ;

  };

#endif /* SCORINGMODELTEST_H_ */
