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

#include "VertebrateMitochondrialGenomeAlphabet.h"
#include "AlignmentSymbol.h"
#include "ScoringModel.h"

#include <Bpp/Seq/Alphabet/AlphabetState.h>

namespace Circal
  {
    VertebrateMitochondrialGenomeAlphabet::VertebrateMitochondrialGenomeAlphabet()
      {
      }

    VertebrateMitochondrialGenomeAlphabet::VertebrateMitochondrialGenomeAlphabet(
        const ScoringModel* scm)
      {
        convertModeltoAlphabet(scm);
      }

    VertebrateMitochondrialGenomeAlphabet::~VertebrateMitochondrialGenomeAlphabet()
      {
      }
    unsigned int VertebrateMitochondrialGenomeAlphabet::getSize() const
      {
        return GenomeAlphabet::getSize();
      }

    unsigned int VertebrateMitochondrialGenomeAlphabet::getNumberOfTypes() const
      {
        return GenomeAlphabet::getNumberOfTypes();
      }

    int VertebrateMitochondrialGenomeAlphabet::getUnknownCharacterCode() const
      {
        return 0;
      }

    std::string VertebrateMitochondrialGenomeAlphabet::getAlphabetType() const
      {
        return "Vertebrate Mitochondrial Genome";
      }

    bool VertebrateMitochondrialGenomeAlphabet::isUnresolved(int state) const
      {
        return ((unsigned int) state > getNumberOfTypes());
      }

    bool VertebrateMitochondrialGenomeAlphabet::isUnresolved(
        const std::string &state) const
      {
        return isUnresolved(charToInt(state));
      }

    void VertebrateMitochondrialGenomeAlphabet::convertModeltoAlphabet(
        const ScoringModel* scm)
      {
//        alphabet_->resize(scm->GetModelSize() + 1);
//
//        int i = 1;
//        unsigned int size;
//        std::string temp = "";
//
//        for (ModelValues::const_iterator itr = scm->constitStart();
//            itr != scm->constitEnd(); itr++)
//          {
//            temp = itr->first;
//            alphabet_[i].num = i - 1;
//            alphabet_[i].letter = temp;
//            alphabet_[i].abbr = temp;
//            alphabet_[i].name = temp;
//            size = itr->second.symbol.length();
//            i++;
////            std::cout << "Chevron <" << temp << "> added" << std::endl;
//          }
        //Add gap code
        bpp::AlphabetState gap(-1, "_", "Gap");
        registerState(gap);

      }
  }

