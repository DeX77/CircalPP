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

#include "AlignmentFactory.h"
#include "ScoringModel.h"
#include "MatrixHelper.h"

namespace Circal
  {

    AlignmentFactory::AlignmentFactory()
      {
        this->matrix = new Circal::MatrixHelper();
      }

    AlignmentFactory::~AlignmentFactory()
      {
      }
    void AlignmentFactory::ForwardRecursionGotoh(const bpp::Sequence* A,
        const bpp::Sequence* B, const ScoringModel* scoreM, ScoreMatrix* D,
        ScoreMatrix* P, ScoreMatrix* Q)
      {
        double gapOpenP;
        double gapExtendP;
        double gapOpenQ;
        double gapExtendQ;
        double diagScore;

        //Forward
        for (uint i=1; i<A->size()+1; i++)
          for (uint j=1; j<B->size()+1; j++)
            {
              //Score of Match eg. Mismatch
              diagScore = D->at(i-1).at(j-1) + scoreM->ScoreOf(A->getChar(i-1), B->getChar(j
                  -1));

              //Score of Open Gap in A
              gapOpenP = D->at(i-1).at(j) + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                  + scoreM->ScoreOfGapExtend(B->getChar(j-1));
              //Score of continue Gap in A
              gapExtendP = P->at(i-1).at(j)+ scoreM->ScoreOfGapExtend(B->getChar(j-1));

              //Score of Open Gap in B
              gapOpenQ = D->at(i).at(j-1) + scoreM->ScoreOfGapOpen(A->getChar(i-1))
                  + scoreM->ScoreOfGapExtend(A->getChar(i-1));
              //Score of continue Gap in B
              gapExtendQ = Q->at(i).at(j-1) + scoreM->ScoreOfGapExtend(A->getChar(i-1));

              //Set Helper Matrices Values
              P->at(i).at(j) = scoreM->BestOfTwo(gapOpenP, gapExtendP);
              Q->at(i).at(j) = scoreM->BestOfTwo(gapOpenQ, gapExtendQ);

              //Set Distance Matrix Value
              D->at(i).at(j) = scoreM->BestOfThree(diagScore, P->at(i).at(j), Q->at(i).at(j) );

            }
      }

    void AlignmentFactory::BacktrackingGotoh(Alignment* out,
        const bpp::Sequence* A, const bpp::Sequence* B,
        const ScoringModel* scoreM, const ScoreMatrix* D, const ScoreMatrix* P,
        const ScoreMatrix* Q)
      {

        int i = D->size()-1;
        int j = D->at(0).size()-1;
        double minScore = D->at(i).at(j);

        vector<int> outA;
        vector<int>::iterator itA = outA.begin();
        vector<int> outB;
        vector<int>::iterator itB = outB.begin();

        minScore = matrix->SearchMinimumPosition(D, i, j, minScore, scoreM);

        minScore = 0;
        //        minScore = std::min(D->at(matrix->SearchMinimuminlastColumn(D)).at(j), std::min(P->at(matrix->SearchMinimuminlastColumn(P)).at(j), Q->at(matrix->SearchMinimuminlastColumn(Q)).at(j) ));

        //        std::cout << "Gefunden in: " << i << ":" << j << " Wert=" << minScore
        //            << std::endl;


        while ( (i>0) && (j>0))
          {
            //                        std::cout << "Aktuelle Score: " << minScore << std::endl;

            //Change to Q
            if (scoreM->BestOfTwo(D->at(i).at(j), Q->at(i).at(j) ) == Q->at(i).at(j))
              {

                //Lonely Gap in A
                if (Q->at(i).at(j) == D->at(i).at(j-1) + scoreM->ScoreOfGapOpen(A->getChar(i -1))
                    + scoreM->ScoreOfGapExtend(A->getChar(i-1)))
                  {
                    minScore +=scoreM->ScoreOfGapOpen(A->getChar(i-1));
                    //                                        std::cout << "Left"<< endl;
                    //                    std::cout << "Score + " << scoreM->ScoreOfGapOpen(A->getChar(i-1))
                    //                    << std::endl;
                  }
                //Continued Gap in A
                else if (Q->at(i).at(j) == Q->at(i).at(j-1) + scoreM->ScoreOfGapExtend(A->getChar(i-1)))
                  {
                    //                                        std::cout << "Left cont "<< endl;
                  }
                else
                  {
                    //                    std::cout << "Changed to Q but not possible?" << std::endl;
                  }
                itA = outA.insert(itA, -1);
                itB = outB.insert(itB, B->getValue(j-1));
                j--;
                minScore +=scoreM->ScoreOfGapExtend(A->getChar(i-1));
                //                std::cout << "Score + " << scoreM->ScoreOfGapExtend(A->getChar(i-1))
                //                << std::endl;
                continue;
              }

            //Change to P
            if (scoreM->BestOfTwo(D->at(i).at(j) , P->at(i).at(j) ) == P->at(i).at(j))
              {
                //Lonely Gap in B
                if (P->at(i).at(j) == D->at(i-1).at(j) + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                    + scoreM->ScoreOfGapExtend(B->getChar(j-1)))
                  {
                    //                                        std::cout << "Up"<< endl;
                    minScore +=scoreM->ScoreOfGapOpen(B->getChar(j-1));
                    //                    std::cout << "Score + " << scoreM->ScoreOfGapOpen(B->getChar(j-1))
                    //                    << std::endl;
                  }
                //Continued Gap in B
                else if (P->at(i).at(j) == P->at(i-1).at(j) + scoreM->ScoreOfGapExtend(B->getChar(j -1)))
                  {
                    //                                        std::cout << "Up cont"<< endl;
                  }
                else
                  {
                    //                    std::cout << "Changed to P but not possible?" << std::endl;
                  }
                minScore +=scoreM->ScoreOfGapExtend(B->getChar(j-1));
                //                std::cout << "Score + " << scoreM->ScoreOfGapExtend(B->getChar(j-1))
                //                << std::endl;
                itA = outA.insert(itA, A->getValue(i-1));
                itB = outB.insert(itB, -1);
                i--;
                continue;

              }

            //Diagonal
            else if (D->at(i).at(j) == D->at(i-1).at(j-1)
                + scoreM->ScoreOf(A->getChar(i-1), B->getChar(j-1)))
              {
                minScore += scoreM->ScoreOf(A->getChar(i-1), B->getChar(j-1));
                //                                                std::cout << "Diag"<< endl;
                //                std::cout << "Score + " << scoreM->ScoreOf(A->getChar(i-1), B->getChar(j-1))
                //                << std::endl;
                itA = outA.insert(itA, A->getValue(i-1));
                itB = outB.insert(itB, B->getValue(j-1));
                i--;
                j--;
                continue;
              }

            else
              {

                std::cout << "42!!!"<< std::endl;
                break;
              }

          }

        //                while (i > 0)
        //          {
        //            //            std::cout << "Overlapp A"<< endl;
        //
        //            if (i == 1)
        //              minScore += scoreM->ScoreOfGapOpen(A->getChar(i-1));
        //            minScore += scoreM->ScoreOfGapExtend(A->getChar(i-1));
        //            itA = outA.insert(itA, A->getValue(i-1));
        //            itB = outB.insert(itB, -1);
        //            i--;
        //          }
        //
        //        while (j > 0)
        //          {
        //            //            std::cout << "Overlapp B"<< endl;
        //            if (j == 1)
        //              minScore += scoreM->ScoreOfGapOpen(B->getChar(j-1));
        //            minScore += scoreM->ScoreOfGapExtend(B->getChar(j-1));
        //            itB = outB.insert(itB, B->getValue(j-1));
        //            itA = outA.insert(itA, -1);
        //            j--;
        //          }

        bpp::Sequence seqA(A->getName(), outA, A->getAlphabet());
        bpp::Sequence seqB(B->getName(), outB, B->getAlphabet());

        //Save Score to Container
        out->set_Score(minScore);

        //Add Alinged Sequences to Container
        out->addSequence(seqA);
        out->addSequence(seqB);

      }

    void AlignmentFactory::ForwardRecursionNMW(const bpp::Sequence* A,
        const bpp::Sequence* B, const ScoringModel* scoreM, ScoreMatrix* D)
      {
        double gapOpenP;
        double gapOpenQ;
        double diagScore;

        //Forward
        for (uint i=1; i<A->size()+1; i++)
          for (uint j=1; j<B->size()+1; j++)
            {
              //Score of Match eg. Mismatch
              diagScore = D->at(i-1).at(j-1) + scoreM->ScoreOf(A->getChar(i-1), B->getChar(j
                  -1));

              //Score of Open Gap in A
              gapOpenP = D->at(i-1).at(j) + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                  + scoreM->ScoreOfGapExtend(B->getChar(j-1));

              //Score of Open Gap in B
              gapOpenQ = D->at(i).at(j-1) + scoreM->ScoreOfGapOpen(A->getChar(i-1))
                  + scoreM->ScoreOfGapExtend(A->getChar(i-1));

              //Set Distance Matrix Value
              D->at(i).at(j) = scoreM->BestOfThree(diagScore, gapOpenP, gapOpenQ);

            }
      }
    void AlignmentFactory::BacktrackingNMW(Alignment* out,
        const bpp::Sequence* A, const bpp::Sequence* B,
        const ScoringModel* scoreM, const ScoreMatrix* D)
      {

        int i = D->size()-1;
        int j = D->at(0).size()-1;

        vector<int> outA;
        vector<int>::iterator itA = outA.begin();
        vector<int> outB;
        vector<int>::iterator itB = outB.begin();

        double minScore = double(0);

        /*
         //Search minimum Score in last row
         for (uint k=1; k<=A->size(); k++)
         {

         if (T[k][B->size()] < minScore)
         {
         i = k;
         j = B->size();
         minScore = T[i][j];
         }
         }

         //Search minimum Score in last column    
         for (uint k=1; k<=B->size(); k++)
         {
         if (T[A->size()][k] < minScore)
         {
         i = A->size();
         j = k;
         minScore = T[i][j];
         }
         }

         std::cout << "Gefunden in: "<< i << ":"<< j << " Wert: "<< T[i][j]<< endl;
         */

        double ScoreD;
        double ScoreDiag;
        double ScoreUpD;
        double ScoreLeftD;

        while ( (i>0) && (j>0))
          {
            ScoreD = D->at(i).at(j);
            ScoreDiag = D->at(i-1).at(j-1);
            ScoreLeftD = D->at(i).at(j-1);
            ScoreUpD = D->at(i-1).at(j);

            //Lonely Gap in A
            if (ScoreD == ScoreLeftD + scoreM->ScoreOfGapOpen(A->getChar(i -1))
                + scoreM->ScoreOfGapExtend(A->getChar(i-1)))
              {
                //std::cout << "Left"<< endl;
                minScore += scoreM->ScoreOfGapOpen(A->getChar(i-1));
                minScore += scoreM->ScoreOfGapExtend(A->getChar(i-1));
                itA = outA.insert(itA, -1);
                itB = outB.insert(itB, B->getValue(j-1));
                j--;
              }

            //Lonely Gap in B
            else if (ScoreD == ScoreUpD
                + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                + scoreM->ScoreOfGapExtend(B->getChar(j-1)))
              {
                //std::cout << "Up"<< endl;
                minScore += scoreM->ScoreOfGapOpen(B->getChar(j-1));
                minScore += scoreM->ScoreOfGapExtend(B->getChar(j-1));
                itA = outA.insert(itA, A->getValue(i-1));
                itB = outB.insert(itB, -1);
                i--;
              }
            //Diagonal
            else if (ScoreD == ScoreDiag + scoreM->ScoreOf(A->getChar(i-1),
                B->getChar(j-1)))
              {
                minScore += scoreM->ScoreOf(A->getChar(i-1), B->getChar(j-1));
                //std::cout << "Diag"<< endl;
                itA = outA.insert(itA, A->getValue(i-1));
                itB = outB.insert(itB, B->getValue(j-1));
                i--;
                j--;
              }

            else
              {
                std::cout << "42!!! Lost in D"<< endl;
              }

          }

        while (i > 0)
          {
            //std::cout << "Overlapp A"<< endl;
            minScore += scoreM->ScoreOfGapOpen(A->getChar(i-1));
            itA = outA.insert(itA, A->getValue(i-1));
            itB = outB.insert(itB, -1);
            i--;
          }

        while (j > 0)
          {
            //std::cout << "Overlapp B"<< endl;
            minScore += scoreM->ScoreOfGapOpen(B->getChar(j-1));
            itB = outB.insert(itB, B->getValue(j-1));
            itA = outA.insert(itA, -1);
            j--;
          }

        bpp::Sequence seqA(A->getName(), outA, A->getAlphabet());
        bpp::Sequence seqB(B->getName(), outB, B->getAlphabet());

        //Save Score to Container
        out->set_Score(minScore);

        //Add Alinged Sequences to Container
        out->addSequence(seqA);
        out->addSequence(seqB);

      }

    void AlignmentFactory::NeedlemanWunschAlignment(Alignment* out,
        const bpp::Sequence* inA, const bpp::Sequence* inB,
        const ScoringModel* scoreM)
      {
        ScoreMatrix D =
            matrix->InitializeScoreMatrixDistances(inA, inB, scoreM);

        //Forward Iteration
        ForwardRecursionNMW(inA, inB, scoreM, &D);

        BacktrackingNMW(out, inA, inB, scoreM, &D);

      }

    void AlignmentFactory::GotohAlignment(Alignment* out,
        const bpp::Sequence* inA, const bpp::Sequence* inB,
        const ScoringModel* scoreM)
      {

        ScoreMatrix D = matrix->InitScoreMatrixWith(inA, inB, 0);
        ScoreMatrix P = matrix->InitScoreMatrixWith(inA, inB, 0);
        ScoreMatrix Q = matrix->InitScoreMatrixWith(inA, inB, 0);

        //Forward Iteration
        ForwardRecursionGotoh(inA, inB, scoreM, &D, &P, &Q);

        BacktrackingGotoh(out, inA, inB, scoreM, &D, &P, &Q);

      }
  }
