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

#ifndef SEQUENCEPROXY_H_
#define SEQUENCEPROXY_H_

#include <Bpp/Seq/Sequence.h>

/*
 *
 */
class SequenceProxy: public bpp::BasicSequence
  {
public:
  SequenceProxy(const std::string& name, const std::vector<int>& sequence,
      const bpp::Alphabet* alpha) throw (bpp::BadIntException);
  SequenceProxy(const bpp::Sequence& sequence);
  SequenceProxy(const std::string& name, const std::string& sequence, const bpp::Alphabet* alpha);

  virtual ~SequenceProxy();

  void shuffle()
    {
      bpp::BasicSequence::shuffle();
    }
  ;
  void deleteElement(unsigned int pos) throw (bpp::IndexOutOfBoundsException)
    {
      bpp::BasicSequence::deleteElement(pos);
    }
  ;
  void deleteElements(unsigned int pos, unsigned int len)
      throw (bpp::IndexOutOfBoundsException)
    {
      bpp::BasicSequence::deleteElements(pos, len);
    }
  ;
  const bpp::Alphabet* getAlphabet() const
    {
      return bpp::BasicSequence::getAlphabet();
    }
  ;
  void setContent(const std::vector<int>& list) throw (bpp::BadIntException)
    {
      bpp::BasicSequence::setContent(list);
    }
  ;
  const std::vector<int>& getContent() const
    {
      return bpp::BasicSequence::getContent();
    }
  ;
  void setElement(unsigned int pos, const std::string& c)
      throw (bpp::BadCharException, bpp::IndexOutOfBoundsException)
    {
      bpp::BasicSequence::setElement(pos, c);
    }
  ;
  void setElement(unsigned int pos, int v) throw (bpp::BadIntException,
      bpp::IndexOutOfBoundsException)
    {
      bpp::BasicSequence::setElement(pos, v);
    }
  ;
  void addElement(const std::string& c) throw (bpp::BadCharException)
    {
      bpp::BasicSequence::addElement(c);
    }
  ;

  void addElement(unsigned int pos, const std::string& c)
      throw (bpp::BadCharException, bpp::IndexOutOfBoundsException)
    {
      bpp::BasicSequence::addElement(pos, c);
    }
  ;
  void addElement(int v) throw (bpp::BadIntException)
    {
      bpp::BasicSequence::addElement(v);
    }
  ;

  void addElement(unsigned int pos, int v) throw (bpp::BadIntException,
      bpp::IndexOutOfBoundsException)
    {
      bpp::BasicSequence::addElement(pos, v);
    }
  ;
  unsigned int size() const
    {
      return bpp::BasicSequence::size();
    }
  ;
  std::string toString() const
    {
      return bpp::BasicSequence::toString();
    }
  ;
  std::string getChar(unsigned int pos) const
      throw (bpp::IndexOutOfBoundsException)
    {
      return bpp::BasicSequence::getChar(pos);
    }
  ;
  const int& operator[](unsigned int i) const
    {
      return bpp::BasicSequence::operator [](i);
    }
  ;
  int& operator[](unsigned int i)
    {
      return bpp::BasicSequence::operator [](i);
    }
  ;
  int getValue(unsigned int pos) const throw (bpp::IndexOutOfBoundsException)
    {
      return bpp::BasicSequence::getValue(pos);
    }
  ;

  };

#endif /* SEQUENCEPROXY_H_ */
