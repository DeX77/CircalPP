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

#ifndef PSEUDOROTATEDSEQUENCETEST_H_
#define PSEUDOROTATEDSEQUENCETEST_H_

#include "PseudoRotatedSequence.h"

#include <Bpp/Seq/Alphabet/DNA.h>

#include <cxxtest/TestSuite.h>

class PseudoRotatedSequenceTest: public CxxTest::TestSuite
  {
public:
  bpp::Alphabet* alpha;
  Circal::SequenceProxy seq;

  PseudoRotatedSequenceTest() :
      alpha(new bpp::DNA()), seq("A", "AAA", alpha)
    {
    }
  ;

  ~PseudoRotatedSequenceTest()
    {
      delete alpha;
    }
  ;

  void testConstructor(void)
    {
      TS_TRACE("Starting constructor test");
      Circal::PseudoRotatedSequence output(seq);

      TS_TRACE("Finishing constructor test");
    }
  ;

  void testDestructor(void)
    {
      TS_TRACE("Starting destructor test");
      Circal::PseudoRotatedSequence* output = new Circal::PseudoRotatedSequence(
          seq);

      delete output;
      TS_TRACE("Finishing destructor test");
    }
  ;

  void testSize(void)
    {
      TS_TRACE("Starting size test");
      Circal::PseudoRotatedSequence output(seq);
      TS_ASSERT_EQUALS(output.size(), seq.size()*2);

      TS_TRACE("Finishing size test");
    }
  ;

//      unsigned int PseudoRotatedSequence::size() const
//        {
//          //Since this Sequence is metarotated return twice the size minus one
//          return ((2 * content_.size()));
//        }
//      int PseudoRotatedSequence::getValue(unsigned int pos) const
//          throw (bpp::IndexOutOfBoundsException)
//        {
//          pos = pos % content_.size();
//          if (pos > this->size())
//            throw bpp::IndexOutOfBoundsException(
//                "SymbolList::getChar. Invalid position.", pos, 0, size() - 1);
//          return content_[pos];
//        }
//
//      std::string PseudoRotatedSequence::getChar(unsigned int pos) const
//          throw (bpp::IndexOutOfBoundsException)
//        {
//
//          pos = pos % content_.size();
//
//          return bpp::BasicSequence::getChar(pos);
//        }
//      const int & PseudoRotatedSequence::operator[](unsigned int i) const
//        {
//          i = i % content_.size();
//          return content_[i];
//        }
//      int & PseudoRotatedSequence::operator[](unsigned int i)
//        {
//          i = i % content_.size();
//          return content_[i];
//        }

  };

#endif /* PSEUDOROTATEDSEQUENCETEST_H_ */
