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

    Alignment* PseudoCircularAlignmentFactory::NeedlemanWunschAlignment(
        const bpp::Sequence* inA, const bpp::Sequence* inB,
        const ScoringModel* scoreM)
      {
        PseudoRotatedSequence* A = new PseudoRotatedSequence(inA);

        return AlignmentFactory::NeedlemanWunschAlignment(A, inB, scoreM);

      }

    Alignment* PseudoCircularAlignmentFactory::GotohAlignment(
        const bpp::Sequence* inA, const bpp::Sequence* inB,
        const ScoringModel* scoreM)
      {
        PseudoRotatedSequence* A = new PseudoRotatedSequence(inA);

        ScoreMatrix D = matrix->InitScoreMatrixWith(A, inB, 0);
        ScoreMatrix P = matrix->InitScoreMatrixWith(A, inB, 0);
        ScoreMatrix Q = matrix->InitScoreMatrixWith(A, inB, 0);
        ScoreMatrix L = matrix->InitScoreMatrixWith(A, inB, 0);

        //Forward Iteration
        ForwardRecursionGotoh(A, inB, scoreM, &D, &P, &Q, &L);

        //Search Starting Node
        int i = D.size()-1;
        int j = D.at(0).size()-1;
        i = matrix->SearchBestInColumn(&D, scoreM, i, j);

        int horizontalStart = i;
        Alignment* temp = BacktrackingGotohLocal(A, inB, scoreM, &D, &P, &Q, i,
            j);
        int horizontalEnd = i;

        if ((horizontalStart-horizontalEnd) > inA->size())
          {
            std::cout << "Alignment Laenge=" << (horizontalStart-horizontalEnd)
            << " richtig waere aber: "
            << inA->size()
            << std::endl;
          }
        return temp;

      }
  }
