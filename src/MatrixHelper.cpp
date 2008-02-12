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
    void MatrixHelper::CutRowFromTo(ScoreMatrix* D, const uint &start,
        const uint &end)
      {
        std::cout << "Start Row: " << start << " End Row: " << end << " Size: "
            << D->size() << std::endl;
        try
          {
            if ( (start> end )|| (end> D->size()))
            throw new MatrixOutofBoundsError(start,end,(uint)D->size());
            D->erase(D->begin()+start, D->begin()+end);
          }
        catch (MatrixOutofBoundsError)
          {

          }
      }

    void MatrixHelper::CutColumnFromTo(ScoreMatrix* D, const uint &start,
        const uint &end)
      {
        std::cout << "Start Column: " << start << " End Column: " << end
            << " Size: " << D->size() << std::endl;
        for (uint i=0; i<=D->size()-1; i++)
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

        D[0][0] = 0;

        D[0][1] = scoreM->ScoreOfGapOpen(B->getChar(0))
            + scoreM->ScoreOfGapExtend(B->getChar(0));
#ifdef _OPENMP            
#pragma omp parallel for
#endif 
        for (uint i=2; i<D[0].size(); i++)
          D[0][i] = D[0][i-1] + scoreM->ScoreOfGapExtend(B->getChar(i-1));

        D[1][0] = scoreM->ScoreOfGapOpen(A->getChar(0))
            + scoreM->ScoreOfGapExtend(A->getChar(0));
#ifdef _OPENMP            
#pragma omp parallel for
#endif         
        for (uint i=2; i<D.size(); i++)
          D[i][0] = D[i-1][0] + scoreM->ScoreOfGapExtend(A->getChar(i-1));

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
        Alignment* pairWiseAlignments, uint biggestSequenceSize)
      {
        //Construct Graph as Adjacence Matrix
        valarray <bool> tmp(false, biggestSequenceSize);
        BoolMatrix G(tmp, pairWiseAlignments->getNumberOfSequences()/2);

        cout << "constructing Graph"<< endl;
        for (uint i=1; i<pairWiseAlignments->getNumberOfSequences(); i +=2)
          {
#ifdef _OPENMP            
#pragma omp parallel for shared(G)
#endif 
            for (uint j=0; j<pairWiseAlignments->getSequence(i)->size(); j++)
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

    int MatrixHelper::SearchMinimuminlastRow(const ScoreMatrix* M,
        const ScoringModel* scoreM)
      {

        int i = M->size()-1;
        int j = M->at(0).size()-1;
        double minScore = M->at(i).at(j);

        for (uint k=1; k<=M->size()-1; k++)
          {
            if (scoreM->BestOfTwo(M->at(k).at(M->at(0).size()-1), minScore) != minScore)
              {
                i = k;
                minScore = M->at(k).at(M->at(0).size()-1);
              }
          }
        return i;

      }
    int MatrixHelper::SearchMinimuminlastColumn(const ScoreMatrix* M,
        const ScoringModel* scoreM)
      {
        int i = M->size()-1;
        int j = M->at(0).size()-1;
        double minScore = M->at(i).at(j);

        for (uint k=1; k<=M->at(0).size()-1; k++)
          {
            if (scoreM->BestOfTwo(M->at(M->size()-1).at(k), minScore) != minScore)
              {
                j = k;
                minScore = M->at(M->size()-1).at(k);
              }
          }
        return j;

      }
    double &MatrixHelper::SearchMinimumPosition(const ScoreMatrix* M, int &i,
        int &j, double &minScore, const ScoringModel* scoreM)
      {
        minScore = M->at(i).at(j);
        int tempI = SearchMinimuminlastRow(M, scoreM);
        int tempJ = SearchMinimuminlastColumn(M, scoreM);

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
