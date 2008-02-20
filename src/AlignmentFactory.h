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
    class ScoringModel;
    class MatrixHelper;
    class Output;
    typedef std::vector< std::vector<double> > ScoreMatrix;

    class AlignmentFactory
      {
  protected:

      MatrixHelper* matrix;
      Output* prettyPrint;

      //Needleman Wunsch Alignment
      virtual void ForwardRecursionNMW(const bpp::Sequence* A,
          const bpp::Sequence* B, const ScoringModel* scoreM, ScoreMatrix* D);
      virtual Alignment* BacktrackingNMW(const bpp::Sequence* outA,
          const bpp::Sequence* outB, const ScoringModel* scoreM,
          const ScoreMatrix* D, int &i, int &j);

      //Gotoh Alignment

      virtual void ForwardRecursionGotoh(const bpp::Sequence* A,
          const bpp::Sequence* B, const ScoringModel* scoreM, ScoreMatrix* D,
          ScoreMatrix* P, ScoreMatrix* Q);
      virtual Alignment* BacktrackingGotohLocal(const bpp::Sequence* outA,
          const bpp::Sequence* outB, const ScoringModel* scoreM,
          const ScoreMatrix* D, const ScoreMatrix* P, const ScoreMatrix* Q,
          int &i, int &j);
      virtual Alignment* BacktrackingGotohGlobal(const bpp::Sequence* outA,
          const bpp::Sequence* outB, const ScoringModel* scoreM,
          const ScoreMatrix* D, const ScoreMatrix* P, const ScoreMatrix* Q,
          int &i, int &j);
  public:
      AlignmentFactory();
      virtual ~AlignmentFactory();

      virtual Alignment* NeedlemanWunschAlignment(const bpp::Sequence* inA,
          const bpp::Sequence* inB, const ScoringModel* scoreM);
      virtual Alignment* GotohAlignment(const bpp::Sequence* inA,
          const bpp::Sequence* inB, const ScoringModel* scoreM);

      };

  }

#endif /*ALIGNMENTFACTORY_H_*/
