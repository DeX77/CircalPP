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
        return alphabet.size();
      }

    unsigned int VertebrateMitochondrialGenomeAlphabet::getNumberOfTypes() const
      {
        return alphabet.size()-1;
      }

    int VertebrateMitochondrialGenomeAlphabet::getUnknownCharacterCode() const
      {
        return 0;
      }

    string VertebrateMitochondrialGenomeAlphabet::getAlphabetType() const
      {
        return "Vertebrate Mitochondrial Genome";
      }

    bool VertebrateMitochondrialGenomeAlphabet::isUnresolved(int state) const
      {
        return ((uint)state > getNumberOfTypes());
      }

    bool VertebrateMitochondrialGenomeAlphabet::isUnresolved(
        const string &state) const
      {
        return isUnresolved(charToInt(state));
      }

    void VertebrateMitochondrialGenomeAlphabet::convertModeltoAlphabet(
        const ScoringModel* scm)
      {
        alphabet.resize(scm->GetModelSize()+1);

        int i=1;
        uint size;
        string temp = "";

        for (ModelValues::const_iterator itr = scm->constitStart(); itr
            != scm->constitEnd(); itr++)
          {
            temp = itr->first;
            alphabet[i].num = i-1;
            alphabet[i].letter = temp;
            alphabet[i].abbr = temp;
            alphabet[i].name = temp;
            size = itr->second.symbol.length();
            i++;
          }
        //Add gap code
        alphabet[0].num = -1;
        alphabet[0].letter = string(size, '_');
        alphabet[0].abbr = string(size, '_');
        alphabet[0].name = "Gap";

      }
  }

