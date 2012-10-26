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

#ifndef MULTIPLEPSEUDOCIRCULARALIGNMENTFACTORYTEST_H_
#define MULTIPLEPSEUDOCIRCULARALIGNMENTFACTORYTEST_H_

#include "MultiplePseudoCircularAlignmentFactory.h"
#include "ScoringModel.h"

#include <Bpp/Seq/Alphabet/DNA.h>

#include <cxxtest/TestSuite.h>

class MultiplePseudoCircularAlignmentFactoryTest: public CxxTest::TestSuite
  {
public:

  void testGotohalignMultiple(void)
    {
      TS_TRACE("Starting GotohalignMultiple test");

      Circal::MultiplePseudoCircularAlignmentFactory tester;
      bpp::DNA* alpha = new bpp::DNA();

      bpp::VectorSequenceContainer* input = new bpp::VectorSequenceContainer(
          alpha);

      Circal::SequenceProxy seqA("A", "AAA", alpha);
      Circal::SequenceProxy seqB("B", "BBB", alpha);

      input->addSequence(seqA);
      input->addSequence(seqB);

      Circal::Alignment output = tester.GotohalignMultiple(input,
          new Circal::ScoringModel(), 0, false);

      TS_TRACE("Finishing GotohalignMultiple test");
    }
  ;

  void testNMWalignMultiple(void)
    {
      TS_TRACE("Starting NMWalignMultiple test");

      Circal::MultiplePseudoCircularAlignmentFactory tester;
      bpp::DNA* alpha = new bpp::DNA();

      bpp::VectorSequenceContainer* input = new bpp::VectorSequenceContainer(
          alpha);

      Circal::SequenceProxy seqA("A", "AAA", alpha);
      Circal::SequenceProxy seqB("B", "BBB", alpha);

      input->addSequence(seqA);
      input->addSequence(seqB);

      Circal::Alignment output = tester.NMWalignMultiple(input,
          new Circal::ScoringModel(), 3, false);

      TS_TRACE("Finishing NMWalignMultiple test");
    }
  ;

  };

#endif /* MULTIPLEPSEUDOCIRCULARALIGNMENTFACTORYTEST_H_ */
