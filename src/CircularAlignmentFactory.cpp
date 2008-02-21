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

#include "CircularAlignmentFactory.h"
#include "ScoringModel.h"
#include "MatrixHelper.h"
#include "RotatedSequence.h"

namespace Circal
  {
    CircularAlignmentFactory::CircularAlignmentFactory()
      {
      }

    CircularAlignmentFactory::~CircularAlignmentFactory()
      {
      }

    Alignment CircularAlignmentFactory::GotohAlignment(const bpp::Sequence* A,
        const bpp::Sequence* B, const ScoringModel* scoreM)
      {

        double bestScore = 0;
        double tempScore = 0;
        uint offset = 0;

        ScoreMatrix D;
        ScoreMatrix P;
        ScoreMatrix Q;

        //#ifdef _OPENMP            
        //#pragma omp parallel for shared(offset)
        //#endif      
        for (uint i=1; i<=A->size(); i++)
          {
            D = matrix->InitScoreMatrixWith(A, B, 0);
            P = matrix->InitScoreMatrixWith(A, B, 0);
            Q = matrix->InitScoreMatrixWith(A, B, 0);

            RotatedSequence rotB(B, i);

            AlignmentFactory::ForwardRecursionGotoh(A, &rotB, scoreM, &D, &P, &Q);
            tempScore = D.at(D.size()-1).at(D.at(0).size()-1);

            if (scoreM->BestOfTwo(tempScore, bestScore) != bestScore)
              {
                bestScore = tempScore;
                offset = i;
              }

          }
        //Stupid but seems the only way
        D = matrix->InitScoreMatrixWith(A, B, 0);
        P = matrix->InitScoreMatrixWith(A, B, 0);
        Q = matrix->InitScoreMatrixWith(A, B, 0);

        RotatedSequence rotB(B, offset);
        AlignmentFactory::ForwardRecursionGotoh(A, &rotB, scoreM, &D, &P, &Q);

        uint i = D.size()-1;
        uint j = D.at(0).size()-1;

        return BacktrackingGotohGlobal(A, &rotB, scoreM, &D, &P, &Q, i, j);

      }

    Alignment CircularAlignmentFactory::NeedlemanWunschAlignment(
        const bpp::Sequence* A, const bpp::Sequence* B,
        const ScoringModel* scoreM)
      {

        double bestScore = 0;
        double tempScore = 0;
        uint offset = 0;

        ScoreMatrix D;

        //#ifdef _OPENMP            
        //#pragma omp parallel for shared(offset)
        //#endif      
        for (uint i=1; i<=A->size(); i++)
          {
            D = matrix->InitScoreMatrixWith(A, B, 0);

            RotatedSequence rotB(B, i);

            AlignmentFactory::ForwardRecursionNMW(A, &rotB, scoreM, &D);
            tempScore = D.at(D.size()-1).at(D.at(0).size()-1);

            if (scoreM->BestOfTwo(tempScore, bestScore) != bestScore)
              {
                bestScore = tempScore;
                offset = i;
              }

          }
        //Stupid but seems the only way
        D = matrix->InitScoreMatrixWith(A, B, 0);

        RotatedSequence rotB(B, offset);
        AlignmentFactory::ForwardRecursionNMW(A, &rotB, scoreM, &D);

        uint i = D.size()-1;
        uint j = D.at(0).size()-1;

        return BacktrackingNMW(A, &rotB, scoreM, &D, i, j);

      }
  }
