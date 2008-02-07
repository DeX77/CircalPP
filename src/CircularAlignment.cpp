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

#include "CircularAlignment.h"
#include "ScoringModel.h"
#include "MatrixHelper.h"
#include "RotatedSequence.h"
#include "Output.h"

namespace Circal
  {
    CircularAlignment::CircularAlignment(const bpp::Alphabet* alpha) :
      bpp::VectorSequenceContainer(alpha), Alignment(alpha)
      {
      }

    CircularAlignment::~CircularAlignment()
      {
      }

    Alignment* CircularAlignment::GotohAlignment(const bpp::Sequence* A,
        const bpp::Sequence* B, const ScoringModel* scoreM)
      {
        Alignment* temp = Alignment::GotohAlignment(A, B, scoreM);
        Alignment* out = temp;
        int score = temp->get_Score();

        for (uint i=1; i<=A->size(); i++)
          {
            temp = Alignment::GotohAlignment(A,
                dynamic_cast<bpp::Sequence*>(new RotatedSequence(B, i)), scoreM);
            if (scoreM->BestOfTwo(temp->get_Score(), score) != score)
              {
                out = temp;
                score = temp->get_Score();
              }
          }
        return out;
      }
    Alignment* CircularAlignment::NeedlemanWunschAlignment(
        const bpp::Sequence* A, const bpp::Sequence* B,
        const ScoringModel* scoreM)
      {
        Alignment* temp = Alignment::NeedlemanWunschAlignment(A, B, scoreM);
        Alignment* out = temp;
        int score = temp->get_Score();

        for (uint i=1; i<A->size(); i++)
          {
            temp = Alignment::NeedlemanWunschAlignment(A,
                dynamic_cast<bpp::Sequence*>(new RotatedSequence(B, i)), scoreM);
            if (temp->get_Score() < score)
              {
                out = temp;
                score = temp->get_Score();
              }
          }
        return out;
      }
  }
