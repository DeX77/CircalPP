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

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <vector>
#include <valarray>
#include <string>

namespace bpp
  {
    class Sequence;
  }

namespace Circal
  {
    class Alignment;
    typedef std::vector< std::vector<double> > ScoreMatrix;
    typedef std::valarray< std::valarray <bool> > BoolMatrix;

    class Output
      {
  public:
      Output();
      virtual ~Output();
      std::string SequencePrettyPrint(bpp::Sequence* A);
      std::string ScoreMatrixPrettyPrint(bpp::Sequence* A, bpp::Sequence* B,
          const ScoreMatrix &D);
      std::string AlignmentPrettyPrint(Alignment* aln);
      std::string AdjacenceMatrixPrettyPrint(const BoolMatrix &G);
      std::string TCoffeeLibFormat(Alignment* aln);

      };
  }
#endif /*OUTPUT_H_*/
