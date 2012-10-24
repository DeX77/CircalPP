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

#ifndef RANDOMSEQUENCETEST_H_
#define RANDOMSEQUENCETEST_H_

#include "RandomSequence.h"

#include <Bpp/Seq/Alphabet/DNA.h>

#include <cxxtest/TestSuite.h>

class RandomSequenceTest: public CxxTest::TestSuite
  {
public:

  bpp::Alphabet* alpha;

  RandomSequenceTest() :
      alpha(new bpp::DNA())
    {

    }
  ;

  void testRandomSequenceConstructor(void)
    {
      TS_TRACE("Starting RandomSequenceConstructor test");

      Circal::RandomSequence test = Circal::RandomSequence(10, alpha);

      TS_TRACE("Finishing RandomSequenceConstructor test");
    }
  ;

  void testRandomSequenceDestructor(void)
    {
      TS_TRACE("Starting RandomSequenceDestructor test");
      Circal::RandomSequence* test = new Circal::RandomSequence(10, alpha);
      delete test;

      TS_TRACE("Finishing RandomSequenceDestructor test");
    }
  ;

  };
#endif /* RANDOMSEQUENCETEST_H_ */
