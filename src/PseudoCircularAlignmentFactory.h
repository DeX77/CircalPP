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

#ifndef PSEUDOCIRCULARALIGNMENT_H_
#define PSEUDOCIRCULARALIGNMENT_H_

#include "AlignmentFactory.h"

namespace Circal
  {
    class PseudoRotatedSequence;

    class PseudoCircularAlignmentFactory : public virtual AlignmentFactory
      {
      //Needleman Wunsch Alignment
      virtual double ForwardRecursionSmithWaterman(const SequenceProxy A,
          const PseudoRotatedSequence B, ScoringModel* scoreM,
          const int &delta, ScoreMatrix3D* D, uint &bi, uint &bj);
      virtual Alignment BacktrackingSmithWaterman(const SequenceProxy A,
          const PseudoRotatedSequence B, ScoringModel* scoreM,
          const int &delta, const ScoreMatrix3D* D, uint &i, uint &j);

      //Gotoh Alignment
      virtual double ForwardRecursionSmithWatermanAffin(const SequenceProxy A,
          const PseudoRotatedSequence B, ScoringModel* scoreM,
          const int &delta, ScoreMatrix3D* D, ScoreMatrix3D* P,
          ScoreMatrix3D* Q, uint &bi, uint &bj);

      virtual Alignment BacktrackingSmithWatermanAffin(const SequenceProxy A,
          const PseudoRotatedSequence B, ScoringModel* scoreM,
          const int &delta, const ScoreMatrix3D* D, const ScoreMatrix3D* P,
          const ScoreMatrix3D* Q, uint &i, uint &j, bool verbose=false);

  public:
      PseudoCircularAlignmentFactory();
      virtual ~PseudoCircularAlignmentFactory();

      Alignment NeedlemanWunschAlignment(const SequenceProxy A,
          const SequenceProxy B, ScoringModel* scoreM, const int &delta,
          bool verbose=false);

      Alignment GotohAlignment(const SequenceProxy inA,
          const SequenceProxy inB, ScoringModel* scoreM, const int &delta,
          bool verbose);
      };
  }

#endif /*PSEUDOCIRCULARALIGNMENT_H_*/
