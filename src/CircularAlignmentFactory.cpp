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

    Alignment* CircularAlignmentFactory::GotohAlignment(const bpp::Sequence* A,
        const bpp::Sequence* B, const ScoringModel* scoreM)
      {

        double bestScore = 0;
        uint offset = 0;
        Alignment* temp = new Alignment(A->getAlphabet());
        Alignment* out = new Alignment(A->getAlphabet());

//#ifdef _OPENMP            
//#pragma omp parallel for shared(offset)
//#endif      
        for (uint i=1; i<=A->size(); i++)
          {
            temp = AlignmentFactory::GotohAlignment(A,
                dynamic_cast<bpp::Sequence*>(new RotatedSequence(B, i)), scoreM);
            if (scoreM->BestOfTwo(temp->get_Score(), bestScore) != bestScore)
              {
                offset = i;
                bestScore = temp->get_Score();
                std::cout << "Neue Highscore " << bestScore << " bei Offset: "
                    << i << std::endl;
                out = temp;

              }

          }
        return out;

      }

    Alignment* CircularAlignmentFactory::NeedlemanWunschAlignment(
        const bpp::Sequence* A, const bpp::Sequence* B,
        const ScoringModel* scoreM)
      {
        double bestScore = 0;
        uint offset = 0;
        Alignment* temp = new Alignment(A->getAlphabet());
        Alignment* out = new Alignment(A->getAlphabet());

#ifdef _OPENMP            
#pragma omp parallel for
#endif      
        for (uint i=1; i<=A->size(); i++)
          {
            temp = AlignmentFactory::NeedlemanWunschAlignment(A,
                dynamic_cast<bpp::Sequence*>(new RotatedSequence(B, i)), scoreM);
            if (scoreM->BestOfTwo(temp->get_Score(), bestScore) != bestScore)
              {
                offset = i;
                bestScore = temp->get_Score();
                std::cout << "Neue Highscore " << bestScore << " bei Offset: "
                    << i << std::endl;
                delete out;
                out = temp;

              }

          }
        return out;

      }
  }
