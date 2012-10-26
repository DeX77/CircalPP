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

#ifndef CORRECTEDFASTATEST_H_
#define CORRECTEDFASTATEST_H_

#include "CorrectedFasta.h"
#include "ScoringModel.h"
#include "VertebrateMitochondrialGenomeAlphabet.h"

#include <cxxtest/TestSuite.h>

class CorrectedFastaTest: public CxxTest::TestSuite
  {

public:
  std::string sm;
  std::stringstream smstream;
  Circal::ScoringModel* scorem;
  bpp::Alphabet* alpha;
  std::string fasta_input;

  CorrectedFastaTest() :
      sm("p 10  1 : A B C D E F G H I J K L M N O P Q R S T U V W X Y Z\n"), smstream(
          sm), scorem(new Circal::ScoringModel(smstream)), alpha(
          new Circal::VertebrateMitochondrialGenomeAlphabet(scorem))
    {
      fasta_input =
          ">Lorem  Ipsum   1\n"
              "Sed     ut      perspiciatis    unde    omnis   iste    natus   error   sit     voluptatem      accusantium     doloremque      laudantium      totam   rem     aperiam eaque   ipsa    quae    ab      illo    inventore       veritatis  et       quasi   architecto      beatae  vitae   dicta   sunt    explicabo\n"
              ">Lorem  Ipsum   2\n"
              "ut      perspiciatis    unde    omnis   iste    natus   error   sit     voluptatem      accusantium     doloremque      laudantium      totam   rem     aperiam eaque   ipsa    quae    ab      illo    inventore       veritatis       et quasi    architecto      beatae  vitae   dicta   sunt    explicabo       Sed\n"
              ">Lorem  Ipsum   3\n"
              "perspiciatis    unde    omnis   iste    natus   error   sit     voluptatem      accusantium     doloremque      laudantium      totam   rem     aperiam eaque   ipsa    quae    ab      illo    inventore       veritatis       et quasi    architecto      beatae  vitae   dicta   sunt    explicabo       Sed        ut\n"
              ">Lorem  Ipsum   4\n"
              "unde    omnis   iste    natus   error   sit     voluptatem      accusantium     doloremque      laudantium      totam   rem     aperiam eaque   ipsa    quae    ab      illo    inventore       veritatis       et quasi    architecto      beatae  vitae   dicta   sunt    explicabo       Sed        ut      perspiciatis\n"
              ">Lorem  Ipsum   5\n"
              "omnis   iste    natus   error   sit     voluptatem      accusantium     doloremque      laudantium      totam   rem     aperiam eaque   ipsa    quae    ab      illo    inventore       veritatis       et quasi    architecto      beatae  vitae   dicta   sunt    explicabo       Sed        ut      perspiciatis    unde\n"
              ">Lorem  Ipsum   6\n"
              "iste    natus   error   sit     voluptatem      accusantium     doloremque      laudantium      totam   rem     aperiam eaque   ipsa    quae    ab      illo    inventore       veritatis       et quasi    architecto      beatae  vitae   dicta   sunt    explicabo       Sed        ut      perspiciatis    unde    omnis\n"
              ">Lorem  Ipsum   7\n"
              "beatae  vitae   dicta   sunt    explicabo       Sed     ut      perspiciatis    unde    omnis   iste    natus   error   sit     voluptatem      accusantium     doloremque      laudantium      totam   rem     aperiam eaque   ipsa    quae    ab      illo    inventore       veritatis  et       quasi   architecto\n";
    }
  ;
  ~CorrectedFastaTest()
    {
      delete alpha;
    }
  ;

  void testappendFromStream(void)
    {
      TS_TRACE("Starting appendFromStream test");

      Circal::CorrectedFasta tester;
      std::istringstream input(fasta_input);
      bpp::VectorSequenceContainer vsc(alpha);

      tester.appendFromStream(input, vsc);

      TS_ASSERT_EQUALS(vsc.getNumberOfSequences(), 7);

      TS_TRACE("Finishing appendFromStream test");

    }
  ;

  };

#endif /* CORRECTEDFASTATEST_H_ */
