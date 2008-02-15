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

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

//Includes from Bio++ 
#include <Seq/AbstractSequenceContainer.h>
#include <Seq/ioseq>
#include <Utils/FileTools.h>

namespace bpp
  {
    class Sequence;
  }

namespace Circal
  {

    class Alignment : public bpp::AlignedSequenceContainer
      {
      uint origSize;
      double Score;

  public:
      Alignment(const bpp::Alphabet* alpha);
      Alignment();
      virtual ~Alignment();
      uint get_origSize() const;
      double get_Score() const;
      void set_Score(const double &s);
      void set_origSize(const uint &orig);
      void addSequence(const bpp::Sequence & sequence, bool checkNames)
          throw(bpp::Exception);
      void addSequence(const bpp::Sequence & sequence) throw(bpp::Exception);
      };
  }
#endif /*ALIGNMENT_H_*/
