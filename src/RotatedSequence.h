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

#ifndef ROTATEDSEQUENCE_H_
#define ROTATEDSEQUENCE_H_

#include <Bpp/Seq/Sequence.h>

namespace bpp
  {
    class IndexOutOfBoundsException;
  }
namespace Circal
  {

    class RotatedSequence : public bpp::Sequence
      {
      uint offset;
  public:

      RotatedSequence(const bpp::Sequence* a, const uint &i);
      virtual ~RotatedSequence();

      /**
       * @brief Get the element at position 'pos' as an int.
       *
       * @param pos The position of the character to retrieve.
       */
      virtual int getValue(unsigned int pos) const
          throw (bpp::IndexOutOfBoundsException);

      /**
       * @brief Get the element at position 'pos' as a character.
       *
       * @param pos The position of the character to retrieve.
       */
      virtual string getChar(unsigned int pos) const
          throw (bpp::IndexOutOfBoundsException);

      /**
       * @name Provide direct access to the list content.
       *
       * @warning These operators allow you to modifiy the list content.
       * No alphabet checking is performed for your modifications, so use with care, or
       * consider using the setContent() method.
       *
       * @{
       */

      /**
       * @brief Operator [] overloaded for quick access to a character in list.
       *
       * @param i The position to retrieve.
       * @return The integer value of character at position i.
       */
      virtual const int & operator[](unsigned int i) const;

      /**
       * @brief Operator [] overloaded for quick access to a character in list.
       *
       * @param i The position to retrieve.
       * @return The integer value of character at position i.
       */
      virtual int & operator[](unsigned int i);

      /** @} */

      /**
       * @brief Convert the list as a string.
       *
       * This method is useful for dumping a list to a file or to the screen for display.
       *
       * @return The whole list as a string.
       */
      std::string toString() const;

      };

  }

#endif /*ROTATEDSEQUENCE_H_*/
