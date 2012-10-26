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

#ifndef VERTEBRATEMITOCHONDRIALGENOMEALPHABETTEST_H_
#define VERTEBRATEMITOCHONDRIALGENOMEALPHABETTEST_H_

#include "VertebrateMitochondrialGenomeAlphabet.h"
#include "ScoringModel.h"

#include <cxxtest/TestSuite.h>

class VertebrateMitochondrialGenomeAlphabetTest: public CxxTest::TestSuite
  {
public:

  void testConstrucor(void)
    {
      TS_TRACE("Starting constructor test");

      Circal::VertebrateMitochondrialGenomeAlphabet test;

      TS_TRACE("Finishing constructor test");
    }
  ;

  void testConstrucorScoringModel(void)
    {
      TS_TRACE("Starting ScoringModel test");

      const std::string sm =
          "p 10  1 : A B C D E F G H I J K L M N O P Q R S T U V W X Y Z\n";
      std::stringstream stringsm(sm);

      Circal::ScoringModel* scm = new Circal::ScoringModel(stringsm);
      Circal::VertebrateMitochondrialGenomeAlphabet test(scm);
      delete scm;

      TS_TRACE("Finishing ScoringModel test");
    }
  ;

  void testDestrucor(void)
    {
      TS_TRACE("Starting Destrucor test");

      Circal::VertebrateMitochondrialGenomeAlphabet* test =
          new Circal::VertebrateMitochondrialGenomeAlphabet();

      TS_TRACE("Finishing Destrucor test");
    }
  ;

//  explicit VertebrateMitochondrialGenomeAlphabet(const ScoringModel* scm);

  void testgetSize(void)
    {
      TS_TRACE("Starting getSize test");

      Circal::VertebrateMitochondrialGenomeAlphabet test;
      TS_ASSERT_EQUALS(test.getSize(), 26);

      TS_TRACE("Finishing getSize test");
    }
  ;

  void testgetNumberOfTypes(void)
    {
      TS_TRACE("Starting getNumberOfTypes test");

      Circal::VertebrateMitochondrialGenomeAlphabet test;
      TS_ASSERT_EQUALS(test.getNumberOfTypes(), 27);

      TS_TRACE("Finishing getNumberOfTypes test");
    }
  ;
  void testgetUnknownCharacterCode(void)
    {
      TS_TRACE("Starting getUnknownCharacterCode test");

      Circal::VertebrateMitochondrialGenomeAlphabet test;
      TS_ASSERT_EQUALS(test.getUnknownCharacterCode(), 0);

      TS_TRACE("Finishing getUnknownCharacterCode test");
    }
  ;

  void testgetAlphabetType(void)
    {
      TS_TRACE("Starting getAlphabetType test");

      Circal::VertebrateMitochondrialGenomeAlphabet test;
      TS_ASSERT_EQUALS(test.getAlphabetType(),
          "Vertebrate Mitochondrial Genome");

      TS_TRACE("Finishing getAlphabetType test");
    }
  ;

  void testisUnresolvedInt(void)
    {
      TS_TRACE("Starting isUnresolved int test");

      Circal::VertebrateMitochondrialGenomeAlphabet test;
      TS_ASSERT_EQUALS(test.isUnresolved(0), false);

      TS_TRACE("Finishing isUnresolved int test");
    }
  ;

  void testisUnresolvedString(void)
    {
      TS_TRACE("Starting isUnresolved String test");

      Circal::VertebrateMitochondrialGenomeAlphabet test;
      TS_ASSERT_EQUALS(test.isUnresolved("dasdasdsa"), false);

      TS_TRACE("Finishing isUnresolved String test");
    }
  ;
  void testconvertModeltoAlphabet(void)
    {
      TS_TRACE("Starting convertModeltoAlphabet test");

      const std::string sm =
          "p 10  1 : A B C D E F G H I J K L M N O P Q R S T U V W X Y Z\n";
      std::stringstream stringsm(sm);

      Circal::ScoringModel* scm = new Circal::ScoringModel(stringsm);
      Circal::VertebrateMitochondrialGenomeAlphabet test(scm);

      TS_ASSERT_EQUALS(test.isUnresolved("dasdasdsa"), false);

      delete scm;
      TS_TRACE("Finishing convertModeltoAlphabet test");
    }
  ;

//  void convertModeltoAlphabet(const ScoringModel* scm);
  };

#endif /* VERTEBRATEMITOCHONDRIALGENOMEALPHABETTEST_H_ */
