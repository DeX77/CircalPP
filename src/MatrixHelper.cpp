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

#include "MatrixHelper.h"
#include "ScoringModel.h"
#include "Alignment.h"
#include <Seq/sequences>
#include "ErrorClasses.h"

namespace Circal
  {
    MatrixHelper::MatrixHelper()
      {
      }

    MatrixHelper::~MatrixHelper()
      {
      }
    void MatrixHelper::CutRowFromTo(ScoreMatrix* D, const int &start,
        const int &end)
      {

        try
          {
            if ( (start> end )|| (end> D->size()))
            throw new MatrixOutofBoundsError(start,end,(int)D->size());
            D->erase(D->begin()+start, D->begin()+end);
          }
        catch (MatrixOutofBoundsError)
          {

          }
      }

    void MatrixHelper::CutColumnFromTo(ScoreMatrix* D, const int &start,
        const int &end)
      {

        for (int i=0; i<=D->size()-1; i++)
          {
            try
              {
                if ( (start> end )|| (end> D->at(i).size()))
                throw new MatrixOutofBoundsError(start,end,D->at(i).size());
                D->at(i).erase(D->at(i).begin()+start, D->at(i).begin()+end);
              }
            catch (MatrixOutofBoundsError)
              {

              }
          }
      }

    ScoreMatrix MatrixHelper::InitializeScoreMatrixDistances(
        const bpp::Sequence* A, const bpp::Sequence* B,
        const ScoringModel* scoreM)
      {
        //Initialize Score Matrix
        ScoreMatrix D = InitScoreMatrixWith(A, B,
            numeric_limits<double>::infinity());

        D.at(0).at(0) = 0;

        D.at(0).at(1) = scoreM->ScoreOfGapOpen(B->getChar(0));
            + scoreM->ScoreOfGapExtend(B->getChar(0));
        //#ifdef _OPENMP            
        //#pragma omp parallel for
        //#endif 
        for (int i=2; i<D.at(0).size(); i++)
          D.at(0).at(i) = D.at(0).at(i-1) + scoreM->ScoreOfGapExtend(B->getChar(i-1));

        D.at(1).at(0) = scoreM->ScoreOfGapOpen(A->getChar(0))
            + scoreM->ScoreOfGapExtend(A->getChar(0));
        //#ifdef _OPENMP            
        //#pragma omp parallel for
        //#endif         
        for (int i=2; i<D.size(); i++)
          D.at(i).at(0) = D.at(i-1).at(0) + scoreM->ScoreOfGapExtend(A->getChar(i-1));

        return D;
      }

    ScoreMatrix MatrixHelper::InitScoreMatrixWith(const bpp::Sequence* A,
        const bpp::Sequence* B, const double &init)
      {
        //Initialize Score Matrix with inf
        ScoreMatrix P(A->size()+1, vector<double>(B->size() +1, init));

        return P;
      }

    BoolMatrix MatrixHelper::CreateAdjacenceGraph(
        Alignment* pairWiseAlignments, int biggestSequenceSize)
      {
        //Construct Graph as Adjacence Matrix
        valarray <bool> tmp(false, biggestSequenceSize);
        BoolMatrix G(tmp, pairWiseAlignments->getNumberOfSequences()/2);

        cout << "constructing Graph"<< endl;
        for (int i=1; i<pairWiseAlignments->getNumberOfSequences(); i +=2)
          {
#ifdef _OPENMP            
#pragma omp parallel for shared(G)
#endif 
            for (int j=0; j<pairWiseAlignments->getSequence(i)->size(); j++)
              {
                //Check for match or not mismatch := !gap
                if (pairWiseAlignments->getSequence(i)->getValue(j) != -1)
                  G[i/2][j] = true;
                else
                  G[i/2][j] = false;
              }
          }

        return G;
      }

    int MatrixHelper::SearchBestInRow(const ScoreMatrix* M,
        const ScoringModel* scoreM, const int &start, const int &row)
      {

        int i = start;
        double minScore = M->at(row).at(start);

        for (int k=1; k<start; k++)
          {
            if (scoreM->BestOfTwo(M->at(row).at(k), minScore) != minScore)
              {
                i = k;
                minScore = M->at(row).at(k);
              }
          }
        return i;

      }
    int MatrixHelper::SearchBestInColumn(const ScoreMatrix* M,
        const ScoringModel* scoreM, const int &start, const int &column)
      {

        int j = start;
        double minScore = M->at(start).at(column);

        for (int k=1; k<start; k++)
          {
            if (scoreM->BestOfTwo(M->at(k).at(column), minScore) != minScore)
              {
                j = k;
                minScore = M->at(k).at(column);
              }
          }
        return j;
      }

    double MatrixHelper::SearchBestPositionFrom(const ScoreMatrix* M, int &i,
        int &j, const ScoringModel* scoreM)
      {

        double minScore = M->at(i).at(j);
        int tempI = SearchBestInRow(M, scoreM, j, i);
        int tempJ = SearchBestInColumn(M, scoreM, i, j);

        if (scoreM->BestOfTwo(minScore, M->at(tempI).at(j)) != minScore)
          if (scoreM->BestOfTwo(M->at(tempI).at(j) , M->at(i).at(tempJ)) == M->at(i).at(tempJ))
            {
              j = tempJ;
              minScore = M->at(i).at(tempJ);
            }
          else
            {
              i = tempI;
              minScore = M->at(tempI).at(j);
            }

        return minScore;
      }
  }
