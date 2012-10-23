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

#ifndef RANDOMSEQUENCE_H_
#define RANDOMSEQUENCE_H_

#include "SequenceProxy.h"

namespace Circal
  {

    class RandomSequence: public SequenceProxy
      {
    public:
      RandomSequence(const uint &sequence_size, const bpp::Alphabet* alpha);
      virtual ~RandomSequence();

      //Virtual overrides
      void setToSizeL(unsigned int newSize)
        {
          BasicSequence::setToSizeL(newSize);
        }
      ;
      std::string toString() const
        {
          return BasicSymbolList::toString();
        }
      ;
      std::string getChar(unsigned int pos) const
          throw (bpp::IndexOutOfBoundsException)
        {
          return BasicSymbolList::getChar(pos);
        }
      ;
      void shuffle()
        {
          BasicSymbolList::shuffle();
        }
      ;
      void setElement(unsigned int pos, const std::string& c)
          throw (bpp::BadCharException, bpp::IndexOutOfBoundsException)
        {
          BasicSymbolList::setElement(pos, c);
        }
      ;
      void setElement(unsigned int pos, int v) throw (bpp::BadIntException,
          bpp::IndexOutOfBoundsException)
        {
          BasicSymbolList::setElement(pos, v);
        }
      ;
      void addElement(const std::string& c) throw (bpp::BadCharException)
        {
          BasicSymbolList::addElement(c);
        }
      ;
      void addElement(unsigned int pos, const std::string& c)
          throw (bpp::BadCharException, bpp::IndexOutOfBoundsException)
        {
          BasicSymbolList::addElement(pos, c);
        }
      ;
      const int& operator[](unsigned int i) const
        {
          return content_[i];
        }
      ;
      int& operator[](unsigned int i)
        {
          return content_[i];
        }
      ;
      const bpp::Alphabet* getAlphabet() const
        {
          return BasicSymbolList::getAlphabet();
        }
      ;
      unsigned int size() const
        {
          return BasicSymbolList::size();
        }
      ;
      const std::vector<int>& getContent() const
        {
          return BasicSymbolList::getContent();
        }
      ;
      void deleteElement(unsigned int pos)
          throw (bpp::IndexOutOfBoundsException)
        {
          BasicSequence::deleteElement(pos);
        }
      ;
      void deleteElements(unsigned int pos, unsigned int len)
          throw (bpp::IndexOutOfBoundsException)
        {
          BasicSequence::deleteElements(pos, len);
        }
      ;
      void addElement(int v) throw (bpp::BadIntException)
        {
          return BasicSequence::addElement(v);
        }
      ;
      void addElement(unsigned int pos, int v) throw (bpp::BadIntException,
          bpp::IndexOutOfBoundsException)
        {
          return BasicSequence::addElement(pos, v);
        }
      ;
      int getValue(unsigned int pos) const
          throw (bpp::IndexOutOfBoundsException)
        {
          return BasicSequence::getValue(pos);
        }
      ;
      };

  }

#endif /*RANDOMSEQUENCE_H_*/
