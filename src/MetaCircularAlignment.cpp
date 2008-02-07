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

#include "MetaCircularAlignment.h"
#include "ScoringModel.h"
#include "MatrixHelper.h"
#include "MetaRotatedSequence.h"
#include "Output.h"

namespace Circal
  {
    MetaCircularAlignment::MetaCircularAlignment(const bpp::Alphabet* alpha) :
      VectorSequenceContainer(alpha), Alignment(alpha),
          CircularAlignment(alpha)
      {
      }

    MetaCircularAlignment::~MetaCircularAlignment()
      {
      }

    Alignment* MetaCircularAlignment::GotohAlignment(const bpp::Sequence* inA,
        const bpp::Sequence* inB, const ScoringModel* scoreM)
      {
        
        MetaRotatedSequence* A = new MetaRotatedSequence(inA);

        ScoreMatrix D = matrix->InitScoreMatrixWith(A, inB, 0);
        ScoreMatrix P = matrix->InitScoreMatrixWith(A, inB, 0);
        ScoreMatrix Q = matrix->InitScoreMatrixWith(A, inB, 0);

        //Forward Iteration
        ForwardRecursionGotoh(A, inB, scoreM, &D, &P, &Q);

        return BacktrackingGotoh(A, inB, scoreM, &D, &P, &Q);

      }
  }
