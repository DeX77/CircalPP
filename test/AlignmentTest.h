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

#ifndef ALIGNMENTTEST_H_
#define ALIGNMENTTEST_H_

#include "Alignment.h"

#include <Bpp/Seq/Alphabet/DNA.h>

#include <cxxtest/TestSuite.h>

class AlignmentTest: public CxxTest::TestSuite
  {
public:

  bpp::Alphabet* alpha;

  AlignmentTest() :
      alpha(new bpp::DNA())
    {

    }
  ;

  ~AlignmentTest()
    {
      delete alpha;
    }
  ;

  void testConstructor(void)
    {
      TS_TRACE("Starting constructor test");
      Circal::Alignment output(alpha);

      TS_TRACE("Finishing constructor test");

    }
  ;

  void testDestructor(void)
    {
      TS_TRACE("Starting destructor test");
      Circal::Alignment* output = new Circal::Alignment(alpha);

      delete output;
      TS_TRACE("Finishing destructor test");
    }

  void testoffsetA(void)
    {
      TS_TRACE("Starting offsetA test");
      Circal::Alignment output(alpha);

      uint offset = 3;
      output.set_offsetA(offset);

      TS_ASSERT_EQUALS(output.get_offsetA(), offset);

      TS_TRACE("Finishing offsetB test");
    }

  void testoffsetB(void)
    {
      TS_TRACE("Starting offsetB test");
      Circal::Alignment output(alpha);

      uint offset = 3;
      output.set_offsetB(offset);

      TS_ASSERT_EQUALS(output.get_offsetB(), offset);

      TS_TRACE("Finishing offsetB test");
    }

  void testScore(void)
    {
      TS_TRACE("Starting Score test");
      Circal::Alignment output(alpha);

      double score = 3.14;
      output.set_Score(score);

      TS_ASSERT_EQUALS(output.get_Score(), score);

      TS_TRACE("Finishing Score test");
    }
  };

#endif /* ALIGNMENTTEST_H_ */
