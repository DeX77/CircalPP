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
    typedef std::vector< std::vector< std::vector<double> > > ScoreMatrix3D;

    class PseudoCircularAlignmentFactory : public virtual AlignmentFactory
      {
      //Needleman Wunsch Alignment
      virtual double ForwardRecursionSmithWaterman(const bpp::Sequence* A,
          const PseudoRotatedSequence* B, const ScoringModel* scoreM,
          const int &delta, ScoreMatrix3D* D, uint &bi, uint &bj);
      virtual Alignment BacktrackingSmithWaterman(const bpp::Sequence* A,
          const PseudoRotatedSequence* B, const ScoringModel* scoreM,
          const int &delta, const ScoreMatrix3D* D, uint &i, uint &j);

      //Gotoh Alignment
      virtual double ForwardRecursionSmithWatermanAffin(const bpp::Sequence* A,
          const PseudoRotatedSequence* B, const ScoringModel* scoreM,
          const int &delta, ScoreMatrix3D* D, ScoreMatrix3D* P,
          ScoreMatrix3D* Q, uint &bi, uint &bj);

      virtual Alignment BacktrackingSmithWatermanAffin(
          const bpp::Sequence* A, const PseudoRotatedSequence* B,
          const ScoringModel* scoreM, const int &delta, const ScoreMatrix3D* D,
          const ScoreMatrix3D* P, const ScoreMatrix3D* Q, uint &i, uint &j);

  public:
      PseudoCircularAlignmentFactory();
      virtual ~PseudoCircularAlignmentFactory();

      Alignment NeedlemanWunschAlignment(const bpp::Sequence* A,
          const bpp::Sequence* B, const ScoringModel* scoreM, const int &delta);

      Alignment GotohAlignment(const bpp::Sequence* A, const bpp::Sequence* B,
          const ScoringModel* scoreM, const int &delta);

      };
  }

#endif /*PSEUDOCIRCULARALIGNMENT_H_*/
