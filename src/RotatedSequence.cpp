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

#include "RotatedSequence.h"
#include <iostream>

namespace Circal
  {

    RotatedSequence::RotatedSequence(const bpp::Sequence* a, const uint &i) :
      bpp::Sequence(*a)
      {
        this->offset = i;
      }

    RotatedSequence::~RotatedSequence()
      {
      }

    int RotatedSequence::getValue(unsigned int pos) const
        throw (bpp::IndexOutOfBoundsException)
      {

        pos = ((pos+offset) % this->size());
        if (pos > this->size())
          throw bpp::IndexOutOfBoundsException("SymbolList::getChar. Invalid position.", pos, 0, size() - 1);
        return _content[pos];
      }

    string RotatedSequence::getChar(unsigned int pos) const
        throw (bpp::IndexOutOfBoundsException)
      {
        pos = ((pos+offset) % this->size());
        if (pos > this->size())
          throw bpp::IndexOutOfBoundsException("SymbolList::getChar. Invalid position.", pos, 0, size() - 1);
        string c = "";
        try
          {
            c = _alphabet -> intToChar(_content[pos]);
          }
        catch(bpp::BadIntException bie)
          {

            //This should never happen!
          }
        return c;
      }

    const int & RotatedSequence::operator[](unsigned int i) const
      {
        i = ((i+offset) % this->size());
        return _content[i];
      }
    int & RotatedSequence::operator[](unsigned int i)
      {
        i = ((i+offset) % this->size());
        return _content[i];
      }
    string RotatedSequence::toString() const
      {
        string result = "";
        for(unsigned int i = 0; i < this->size(); i++)
          {
            result += _alphabet -> intToChar(getValue(i));
          }
        return result;
      };

  }
