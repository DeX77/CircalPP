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

#ifndef OUTPUTTEST_H_
#define OUTPUTTEST_H_

#include "Output.h"
#include "Alignment.h"

#include <Bpp/Seq/Alphabet/DNA.h>

#include <cxxtest/TestSuite.h>

class OutputTest: public CxxTest::TestSuite
  {
public:

  void testSequencePrettyPrint(void)
    {
      TS_TRACE("Starting SequencePrettyPrint test");
      Circal::Output tester;
      Circal::SequenceProxy A("A", "AAAA", new bpp::DNA());

      std::string result = tester.SequencePrettyPrint(A);
      TS_TRACE("Finishing SequencePrettyPrint test");
    }
  ;

  void testAlignmentPrettyPrint(void)
    {
      TS_TRACE("Starting AlignmentPrettyPrint test");
      Circal::Output tester;
      Circal::SequenceProxy A("A", "AAAA", new bpp::DNA());
      Circal::SequenceProxy B("B", "BBBB", new bpp::DNA());
      Circal::Alignment* aln = new Circal::Alignment(new bpp::DNA());

      aln->addSequence(A);
      aln->addSequence(B);

      std::string result = tester.AlignmentPrettyPrint(aln);
      TS_TRACE("Finishing AlignmentPrettyPrint test");
    }
  ;

  void testScoreMatrixPrettyPrint(void)
    {
      TS_TRACE("Starting ScoreMatrixPrettyPrint test");
      Circal::Output tester;
      Circal::SequenceProxy A("A", "AAAA", new bpp::DNA());
      Circal::SequenceProxy B("B", "BBBB", new bpp::DNA());
      Circal::ScoreMatrix sm;

      std::string result = tester.ScoreMatrixPrettyPrint(A, B, sm);

      TS_TRACE("Finishing ScoreMatrixPrettyPrint test");
    }
  ;

  void testAdjacenceMatrixPrettyPrint(void)
    {
      TS_TRACE("Starting AdjacenceMatrixPrettyPrint test");
      Circal::Output tester;
      Circal::BoolMatrix bm;

      std::string result = tester.AdjacenceMatrixPrettyPrint(bm);

      TS_TRACE("Finishing AdjacenceMatrixPrettyPrint test");
    }
  ;

  void testTCoffeeLibHeader(void)
    {
      TS_TRACE("Starting TCoffeeLibHeader test");
      Circal::Output tester;
      Circal::SequenceProxy A("A", "AAAA", new bpp::DNA());
      Circal::SequenceProxy B("B", "BBBB", new bpp::DNA());
      Circal::Alignment* aln = new Circal::Alignment(new bpp::DNA());

      aln->addSequence(A);
      aln->addSequence(B);

      std::string result = tester.TCoffeeLibHeader(aln);
      TS_TRACE("Finishing TCoffeeLibHeader test");
    }
  ;

  void testTCoffeeAlignFormat(void)
    {
      TS_TRACE("Starting TCoffeeAlignFormat test");
      Circal::Output tester;
      Circal::SequenceProxy A("A", "AAAA", new bpp::DNA());
      Circal::SequenceProxy B("B", "BBBB", new bpp::DNA());
      Circal::Alignment* aln = new Circal::Alignment(new bpp::DNA());

      aln->addSequence(A);
      aln->addSequence(B);

      std::string result = tester.TCoffeeAlignFormat(aln, aln);
      TS_TRACE("Finishing TCoffeeAlignFormat test");
    }
  ;

  void testTCoffeeLibFooter(void)
    {
      TS_TRACE("Starting TCoffeeLibFooter test");
      Circal::Output tester;

      TS_ASSERT_EQUALS(tester.TCoffeeLibFooter(), "CPU 0\n! SEQ_1_TO_N\n");
      TS_TRACE("Finishing TCoffeeLibFooter test");
    }
  ;
  void testTCoffeeLibFormat(void)
    {
      TS_TRACE("Starting TCoffeeLibFormat test");
      Circal::Output tester;
      Circal::SequenceProxy A("A", "AAAA", new bpp::DNA());
      Circal::SequenceProxy B("B", "BBBB", new bpp::DNA());
      Circal::Alignment* aln = new Circal::Alignment(new bpp::DNA());

      aln->addSequence(A);
      aln->addSequence(B);

      std::string result = tester.TCoffeeLibFormat(aln, aln);
      TS_TRACE("Finishing TCoffeeLibFormat test");
    }
  ;

  };

#endif /* OUTPUTTEST_H_ */
