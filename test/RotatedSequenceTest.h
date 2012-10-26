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

#ifndef ROTATEDSEQUENCETEST_H_
#define ROTATEDSEQUENCETEST_H_

#include "RotatedSequence.h"

#include <Bpp/Seq/Alphabet/DNA.h>

#include <cxxtest/TestSuite.h>

class RotatedSequenceTest: public CxxTest::TestSuite
  {
public:
  bpp::Alphabet* alpha;
  Circal::SequenceProxy seq;

  RotatedSequenceTest() :
      alpha(new bpp::DNA()), seq("A", "AAAA", alpha)
    {
    }
  ;

  ~RotatedSequenceTest()
    {
      delete alpha;
    }
  ;

  void testConstructor(void)
    {
      TS_TRACE("Starting constructor test");
      Circal::RotatedSequence output(seq, 3);

      TS_TRACE("Finishing constructor test");
    }
  ;

  void testDestructor(void)
    {
      TS_TRACE("Starting destructor test");
      Circal::RotatedSequence* output = new Circal::RotatedSequence(seq, 3);

      delete output;
      TS_TRACE("Finishing destructor test");
    }
  ;
  void testtoString(void)
    {
      TS_TRACE("Starting toString test");
      Circal::RotatedSequence output(seq, 3);

      TS_ASSERT_EQUALS(output.toString(), "AAAA");

      TS_TRACE("Finishing toString test");
    }
  ;
  void testgetValue(void)
    {
      TS_TRACE("Starting getValue test");
      Circal::RotatedSequence output(seq, 3);

      TS_ASSERT_EQUALS(output.getValue(0), 0);

      TS_TRACE("Finishing getValue test");
    }
  ;
  void testgetChar(void)
    {
      TS_TRACE("Starting getChar test");
      Circal::RotatedSequence output(seq, 3);

      TS_ASSERT_EQUALS(output.getChar(0), "A");

      TS_TRACE("Finishing getChar test");
    }
  ;

  void testSize(void)
    {
      TS_TRACE("Starting size test");
      Circal::RotatedSequence output(seq, 3);

      TS_ASSERT_EQUALS(output.size(), 4);

      TS_TRACE("Finishing size test");
    }
  ;

  };

#endif /* RotatedSequenceTEST_H_ */
