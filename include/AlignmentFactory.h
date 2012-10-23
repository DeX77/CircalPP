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
#include "MatrixHelper.h"
#include "Output.h"
#include <vector>

namespace Circal
  {
    class ScoringModel;    
    class Output;

    class AlignmentFactory
      {
  protected:

      MatrixHelper matrix;
      Output prettyPrint;

      //Needleman Wunsch Alignment
      virtual void ForwardRecursionNMW(const SequenceProxy A,
          const SequenceProxy B, ScoringModel* scoreM, ScoreMatrix* D);
      virtual Alignment BacktrackingNMW(const SequenceProxy outA,
          const SequenceProxy outB, ScoringModel* scoreM,
          const ScoreMatrix* D, uint &i, uint &j);

      //Gotoh Alignment
      virtual void ForwardRecursionGotoh(const SequenceProxy A,
          const SequenceProxy B, ScoringModel* scoreM, ScoreMatrix* D,
          ScoreMatrix* P, ScoreMatrix* Q);
      virtual Alignment BacktrackingGotohGlocal(const SequenceProxy outA,
          const SequenceProxy outB, ScoringModel* scoreM,
          const ScoreMatrix* D, const ScoreMatrix* P, const ScoreMatrix* Q,
          uint &i, uint &j);
      virtual Alignment BacktrackingGotohGlobal(const SequenceProxy outA,
          const SequenceProxy outB, ScoringModel* scoreM,
          const ScoreMatrix* D, const ScoreMatrix* P, const ScoreMatrix* Q,
          uint &i, uint &j);

      //Smith-Waterman Alignment
      virtual double ForwardRecursionSmithWaterman(const SequenceProxy A,
          const SequenceProxy B, ScoringModel* scoreM, ScoreMatrix* D,
          uint &i, uint &j);
      virtual Alignment BacktrackingSmithWaterman(const SequenceProxy A,
          const SequenceProxy B, ScoringModel* scoreM, const ScoreMatrix* D,
          uint &i, uint &j);

      //Smith-Waterman Alignment-Affin
      virtual double ForwardRecursionSmithWatermanAffin(const SequenceProxy A,
          const SequenceProxy B, ScoringModel* scoreM, ScoreMatrix* D,
          ScoreMatrix* P, ScoreMatrix* Q, uint &bi, uint &bj);
      virtual Alignment BacktrackingSmithWatermanAffin(const SequenceProxy A,
          const SequenceProxy B, ScoringModel* scoreM, const ScoreMatrix* D,
          const ScoreMatrix* P, const ScoreMatrix* Q, uint &i, uint &j);

  public:
      AlignmentFactory();
      virtual ~AlignmentFactory();

      virtual Alignment NeedlemanWunschAlignment(const SequenceProxy inA,
          const SequenceProxy inB, ScoringModel* scoreM, bool verbose=false);
      virtual Alignment GotohAlignment(const SequenceProxy inA,
          const SequenceProxy inB, ScoringModel* scoreM, bool verbose=false);

      virtual Alignment SmithWaterman(const SequenceProxy inA,
          const SequenceProxy inB, ScoringModel* scoreM, bool verbose=false);
      virtual Alignment SmithWatermanAffin(const SequenceProxy inA,
          const SequenceProxy inB, ScoringModel* scoreM, bool verbose=false);

      };

  }

#endif /*ALIGNMENTFACTORY_H_*/
