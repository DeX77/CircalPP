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

#ifndef MATRIXHELPER_H_
#define MATRIXHELPER_H_

#include <vector>
#include <valarray>

namespace bpp
  {
    class Sequence;
  }
namespace Circal
  {
    class Alignment;
    class ScoringModel;
    typedef std::valarray< std::valarray <bool> > BoolMatrix;
    typedef std::vector< std::vector<double> > ScoreMatrix;

    class MatrixHelper
      {
  public:
      MatrixHelper();
      virtual ~MatrixHelper();

      ScoreMatrix InitializeScoreMatrixDistances(const bpp::Sequence* A,
          const bpp::Sequence* B, const ScoringModel* scoreM);

      ScoreMatrix InitScoreMatrixWith(const bpp::Sequence* A,
          const bpp::Sequence* B, const double &init);

      BoolMatrix CreateAdjacenceGraph(Alignment* pairWiseAlignments,
          int biggestSequenceSize);
      void CutRowFromTo(ScoreMatrix* D, const int &start, const int &end);
      void CutColumnFromTo(ScoreMatrix* D, const int &start, const int &end);

      double SearchBestPositionFrom(const ScoreMatrix* M, int &i, int &j,
          const ScoringModel* scoreM);
      int SearchBestInRow(const ScoreMatrix* M, const ScoringModel* scoreM,
          const int &start, const int &row);
      int SearchBestInColumn(const ScoreMatrix* M, const ScoringModel* scoreM,
          const int &start, const int &column);
      };
  }
#endif /*MATRIXHELPER_H_*/
