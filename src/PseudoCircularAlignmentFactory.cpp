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

#include "PseudoCircularAlignmentFactory.h"
#include "ScoringModel.h"
#include "MatrixHelper.h"
#include "PseudoRotatedSequence.h"
#include "Output.h"

namespace Circal
  {
    PseudoCircularAlignmentFactory::PseudoCircularAlignmentFactory()
      {
      }

    PseudoCircularAlignmentFactory::~PseudoCircularAlignmentFactory()
      {
      }

    void PseudoCircularAlignmentFactory::GotohAlignment(Alignment* out,
        const bpp::Sequence* inA, const bpp::Sequence* inB,
        const ScoringModel* scoreM)
      {
        PseudoRotatedSequence* A = new PseudoRotatedSequence(inA);

        ScoreMatrix D = matrix->InitScoreMatrixWith(A, inB, 0);
        ScoreMatrix P = matrix->InitScoreMatrixWith(A, inB, 0);
        ScoreMatrix Q = matrix->InitScoreMatrixWith(A, inB, 0);

        //Forward Iteration
        ForwardRecursionGotoh(A, inB, scoreM, &D, &P, &Q);

        BacktrackingGotoh(out, A, inB, scoreM, &D, &P, &Q);

      }
    void PseudoCircularAlignmentFactory::NeedlemanWunschAlignment(
        Alignment* out, const bpp::Sequence* inA, const bpp::Sequence* inB,
        const ScoringModel* scoreM)
      {
        PseudoRotatedSequence* A = new PseudoRotatedSequence(inA);
        ScoreMatrix D = matrix->InitializeScoreMatrixDistances(A, inB, scoreM);

        //Forward Iteration
        ForwardRecursionNMW(A, inB, scoreM, &D);

        BacktrackingNMW(out, A, inB, scoreM, &D);
      }
  }
