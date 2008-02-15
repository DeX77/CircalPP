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

#ifndef ALIGNMENTFACTORY_H_
#define ALIGNMENTFACTORY_H_

#include "Alignment.h"
#include <vector>

namespace Circal
  {
    class Alignment;
    class ScoringModel;
    class MatrixHelper;
    typedef std::vector< std::vector<double> > ScoreMatrix;

    class AlignmentFactory
      {
  protected:

      MatrixHelper* matrix;

      //Needleman Wunsch Alignment


      virtual void ForwardRecursionNMW(const bpp::Sequence* A,
          const bpp::Sequence* B, const ScoringModel* scoreM, ScoreMatrix* D);
      virtual void BacktrackingNMW(Alignment* out, const bpp::Sequence* outA,
          const bpp::Sequence* outB, const ScoringModel* scoreM,
          const ScoreMatrix* D);

      //Gotoh Alignment

      virtual void ForwardRecursionGotoh(const bpp::Sequence* A,
          const bpp::Sequence* B, const ScoringModel* scoreM, ScoreMatrix* D,
          ScoreMatrix* P, ScoreMatrix* Q);
      virtual void BacktrackingGotoh(Alignment* out, const bpp::Sequence* outA,
          const bpp::Sequence* outB, const ScoringModel* scoreM,
          const ScoreMatrix* D, const ScoreMatrix* P, const ScoreMatrix* Q);
  public:
      AlignmentFactory();
      virtual ~AlignmentFactory();

      virtual void NeedlemanWunschAlignment(Alignment* out, const bpp::Sequence* inA,
          const bpp::Sequence* inB, const ScoringModel* scoreM);
      virtual void GotohAlignment(Alignment* out, const bpp::Sequence* inA,
          const bpp::Sequence* inB, const ScoringModel* scoreM);

      };

  }

#endif /*ALIGNMENTFACTORY_H_*/
