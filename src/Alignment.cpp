/*
 Circal++  - Calculates Multiple conservative Alignments of Circular Sequences
 Copyright (C) 2007  Daniel Exner
 <dex@dragonslave.de>>

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

#include "Alignment.h"
#include "ScoringModel.h"
#include "MatrixHelper.h"
#include "Output.h"

namespace Circal
  {
    Alignment::Alignment(const bpp::Alphabet* alpha) :
      VectorSequenceContainer(alpha), AlignedSequenceContainer(alpha)

      {
        this->set_Score(-numeric_limits<int>::infinity());

      }

    Alignment::~Alignment()
      {

      }
    uint Alignment::get_origSize()
      {
        return this->origSize;
      }
    void Alignment::set_origSize(const uint &orig)
      {
        this->origSize = orig;
      }

    double Alignment::get_Score(void)
      {
        return this->Score;
      }

    void Alignment::set_Score(const double &s)
      {
        this->Score = s;
      }

    void Alignment::addSequence(const bpp::Sequence* sequence)
      {
        //Override stupid Memory leaks
        //bpp::VectorSequenceContainer::addSequence(sequence, checkNames);

        _sequences.push_back(sequence->clone());

      }

  }
