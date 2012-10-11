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
#include "Output.h"
#include <algorithm>

namespace Circal
  {

    AlignmentFactory::AlignmentFactory()
      {
        this->matrix = new Circal::MatrixHelper();
      }

    AlignmentFactory::~AlignmentFactory()
      {
        delete matrix;
      }
    void AlignmentFactory::ForwardRecursionGotoh(const bpp::Sequence* A,
        const bpp::Sequence* B, ScoringModel* scoreM, ScoreMatrix* D,
        ScoreMatrix* P, ScoreMatrix* Q)
      {

        double gapOpenP = 0;
        double gapExtendP = 0;
        double gapOpenQ = 0;
        double gapExtendQ = 0;
        double diagScore = 0;

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

    Alignment AlignmentFactory::BacktrackingGotohGlocal(const bpp::Sequence* A,
        const bpp::Sequence* B, ScoringModel* scoreM, const ScoreMatrix* D,
        const ScoreMatrix* P, const ScoreMatrix* Q, uint &i, uint &j)
      {

        Alignment out(A->getAlphabet());

        std::vector<int> outA;
        std::vector<int>::iterator itA = outA.begin();
        std::vector<int> outB;
        std::vector<int>::iterator itB = outB.begin();

        double minScore = 0;

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
                    //                                        std::cout << "Left"<< std::endl;
                    //                    std::cout << "Score + " << scoreM->ScoreOfGapOpen(A->getChar(i-1))
                    //                    << std::endl;
                  }
                //Continued Gap in A
                else if (Q->at(i).at(j) == Q->at(i).at(j-1) + scoreM->ScoreOfGapExtend(A->getChar(i-1)))
                  {
                    //                                        std::cout << "Left cont "<< std::endl;
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
                    //                                        std::cout << "Up"<< std::endl;
                    minScore +=scoreM->ScoreOfGapOpen(B->getChar(j-1));
                    //                    std::cout << "Score + " << scoreM->ScoreOfGapOpen(B->getChar(j-1))
                    //                    << std::endl;
                  }
                //Continued Gap in B
                else if (P->at(i).at(j) == P->at(i-1).at(j) + scoreM->ScoreOfGapExtend(B->getChar(j -1)))
                  {
                    //                                        std::cout << "Up cont"<< std::endl;
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
                //                                                std::cout << "Diag"<< std::endl;
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

        //        while (j> 0)
        //          {
        //            //            std::cout << "Overlapp B"<< std::endl;
        //            if (j == 1)
        //              minScore += scoreM->ScoreOfGapOpen(B->getChar(j-1));
        //            minScore += scoreM->ScoreOfGapExtend(B->getChar(j-1));
        //            itB = outB.insert(itB, B->getValue(j-1));
        //            itA = outA.insert(itA, -1);
        //            j--;
        //          }

        bpp::Sequence* seqA = new bpp::Sequence(A->getName(), outA, A->getAlphabet());
        bpp::Sequence* seqB = new bpp::Sequence(B->getName(), outB, B->getAlphabet());

        //Save Score to Container
        out.set_Score(minScore);

        //Add Alinged Sequences to Container
        out.addSequence(seqA);
        out.addSequence(seqB);

        delete seqA;
        delete seqB;

        return out;

      }
    Alignment AlignmentFactory::BacktrackingGotohGlobal(const bpp::Sequence* A,
        const bpp::Sequence* B, ScoringModel* scoreM, const ScoreMatrix* D,
        const ScoreMatrix* P, const ScoreMatrix* Q, uint &i, uint &j)
      {

        Alignment out(A->getAlphabet());

        std::vector<int> outA;
        std::vector<int>::iterator itA = outA.begin();
        std::vector<int> outB;
        std::vector<int>::iterator itB = outB.begin();

        double minScore = 0;

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
                    //                                        std::cout << "Left"<< std::endl;
                    //                    std::cout << "Score + " << scoreM->ScoreOfGapOpen(A->getChar(i-1))
                    //                    << std::endl;
                  }
                //Continued Gap in A
                else if (Q->at(i).at(j) == Q->at(i).at(j-1) + scoreM->ScoreOfGapExtend(A->getChar(i-1)))
                  {
                    //                                        std::cout << "Left cont "<< std::endl;
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
                    //                                        std::cout << "Up"<< std::endl;
                    minScore +=scoreM->ScoreOfGapOpen(B->getChar(j-1));
                    //                    std::cout << "Score + " << scoreM->ScoreOfGapOpen(B->getChar(j-1))
                    //                    << std::endl;
                  }
                //Continued Gap in B
                else if (P->at(i).at(j) == P->at(i-1).at(j) + scoreM->ScoreOfGapExtend(B->getChar(j -1)))
                  {
                    //                                        std::cout << "Up cont"<< std::endl;
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
                //                                                std::cout << "Diag"<< std::endl;
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

        while (i > 0)
          {

            //            std::cout << "Overlapp A"<< std::endl;

            if (i == 1)
              minScore += scoreM->ScoreOfGapOpen(A->getChar(i-1));
            minScore += scoreM->ScoreOfGapExtend(A->getChar(i-1));
            itA = outA.insert(itA, A->getValue(i-1));
            itB = outB.insert(itB, -1);
            i--;
          }

        while (j> 0)
          {
            //            std::cout << "Overlapp B"<< std::endl;
            if (j == 1)
              minScore += scoreM->ScoreOfGapOpen(B->getChar(j-1));
            minScore += scoreM->ScoreOfGapExtend(B->getChar(j-1));
            itB = outB.insert(itB, B->getValue(j-1));
            itA = outA.insert(itA, -1);
            j--;
          }

        bpp::Sequence* seqA = new bpp::Sequence(A->getName(), outA, A->getAlphabet());
        bpp::Sequence* seqB = new bpp::Sequence(B->getName(), outB, B->getAlphabet());

        //Save Score to Container
        out.set_Score(minScore);

        //Add Alinged Sequences to Container
        out.addSequence(seqA);
        out.addSequence(seqB);

        delete seqA;
        delete seqB;

        return out;

      }
    void AlignmentFactory::ForwardRecursionNMW(const bpp::Sequence* A,
        const bpp::Sequence* B, ScoringModel* scoreM, ScoreMatrix* D)
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
    Alignment AlignmentFactory::BacktrackingNMW(const bpp::Sequence* A,
        const bpp::Sequence* B, ScoringModel* scoreM, const ScoreMatrix* D,
        uint &i, uint &j)
      {

        std::vector<int> outA;
        std::vector<int>::iterator itA = outA.begin();
        std::vector<int> outB;
        std::vector<int>::iterator itB = outB.begin();

        Alignment out(A->getAlphabet());

        double minScore = double(0);

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
                //std::cout << "Left"<< std::endl;
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
                //std::cout << "Up"<< std::endl;
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
                //std::cout << "Diag"<< std::endl;
                itA = outA.insert(itA, A->getValue(i-1));
                itB = outB.insert(itB, B->getValue(j-1));
                i--;
                j--;
              }

            else
              {
                std::cout << "42!!! Lost in D"<< std::endl;
              }

          }

        while (i > 0)
          {
            //std::cout << "Overlapp A"<< std::endl;
            minScore += scoreM->ScoreOfGapOpen(A->getChar(i-1));
            itA = outA.insert(itA, A->getValue(i-1));
            itB = outB.insert(itB, -1);
            i--;
          }

        while (j > 0)
          {
            //std::cout << "Overlapp B"<< std::endl;
            minScore += scoreM->ScoreOfGapOpen(B->getChar(j-1));
            itB = outB.insert(itB, B->getValue(j-1));
            itA = outA.insert(itA, -1);
            j--;
          }

        bpp::Sequence* seqA = new bpp::Sequence(A->getName(), outA, A->getAlphabet());
        bpp::Sequence* seqB = new bpp::Sequence(B->getName(), outB, B->getAlphabet());

        //Save Score to Container
        out.set_Score(minScore);

        //Add Alinged Sequences to Container
        out.addSequence(seqA);
        out.addSequence(seqB);

        delete seqA;
        delete seqB;

        return out;

      }
    double AlignmentFactory::ForwardRecursionSmithWaterman(
        const bpp::Sequence* A, const bpp::Sequence* B, ScoringModel* scoreM,
        ScoreMatrix* D, uint &bi, uint &bj)
      {
        double gapOpenP;
        double gapOpenQ;
        double diagScore;

        double bestScore = 0;

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
              if (scoreM->BestOfThree(diagScore, gapOpenP, gapOpenQ) > 0)
                D->at(i).at(j) = scoreM->BestOfThree(diagScore, gapOpenP,
                    gapOpenQ);
              else
                D->at(i).at(j) = 0;

              //Save best Score
              if (scoreM->BestOfTwo(bestScore, D->at(i).at(j)) != bestScore)
                {
                  bestScore = D->at(i).at(j);
                  bi = i;
                  bj = j;
                }
            }
        return bestScore;

      }

    Alignment AlignmentFactory::BacktrackingSmithWaterman(
        const bpp::Sequence* A, const bpp::Sequence* B, ScoringModel* scoreM,
        const ScoreMatrix* D, uint &i, uint &j)
      {
        std::vector<int> outA;
        std::vector<int>::iterator itA = outA.begin();
        std::vector<int> outB;
        std::vector<int>::iterator itB = outB.begin();

        Alignment out(A->getAlphabet());

        double minScore = double(0);

        double ScoreD;
        double ScoreDiag;
        double ScoreUpD;
        double ScoreLeftD;

        while ( (i>0) && (j>0))
          {
            ScoreD = D->at(i).at(j);

            if (ScoreD == 0)
              break;

            ScoreDiag = D->at(i-1).at(j-1);
            ScoreLeftD = D->at(i).at(j-1);
            ScoreUpD = D->at(i-1).at(j);

            //Lonely Gap in A
            if (ScoreD == ScoreLeftD + scoreM->ScoreOfGapOpen(A->getChar(i -1))
                + scoreM->ScoreOfGapExtend(A->getChar(i-1)))
              {
                //std::cout << "Left"<< std::endl;
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
                //std::cout << "Up"<< std::endl;
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
                //std::cout << "Diag"<< std::endl;
                itA = outA.insert(itA, A->getValue(i-1));
                itB = outB.insert(itB, B->getValue(j-1));
                i--;
                j--;
              }

            else
              {
                std::cout << "42!!! Lost in D"<< std::endl;
              }

          }

        bpp::Sequence* seqA = new bpp::Sequence(A->getName(), outA, A->getAlphabet());
        bpp::Sequence* seqB = new bpp::Sequence(B->getName(), outB, B->getAlphabet());

        //Save Score to Container
        out.set_Score(minScore);

        //Add Alinged Sequences to Container
        out.addSequence(seqA);
        out.addSequence(seqB);

        delete seqA;
        delete seqB;

        return out;
      }

    double AlignmentFactory::ForwardRecursionSmithWatermanAffin(
        const bpp::Sequence* A, const bpp::Sequence* B, ScoringModel* scoreM,
        ScoreMatrix* D, ScoreMatrix* P, ScoreMatrix* Q, uint &bi, uint &bj)
      {
        double gapOpenP = 0;
        double gapExtendP = 0;
        double gapOpenQ = 0;
        double gapExtendQ = 0;
        double diagScore = 0;

        double bestScore = 0;

        double gapExA;
        double gapExB;

        //Forward
        for (uint i=1; i<A->size()+1; i++)
          for (uint j=1; j<B->size()+1; j++)
            {

              //Score of Match eg. Mismatch
              diagScore = D->at(i-1).at(j-1) + scoreM->ScoreOf(A->getChar(i-1), B->getChar(j
                  -1));

              gapExB = scoreM->ScoreOfGapExtend(B->getChar(j-1));
              gapExA = scoreM->ScoreOfGapExtend(A->getChar(i-1));

              //Score of Open Gap in A
              gapOpenP = D->at(i-1).at(j);
              gapOpenP += scoreM->ScoreOfGapOpen(B->getChar(j-1));
              gapOpenP += gapExB;

              //Score of continue Gap in A
              gapExtendP = P->at(i-1).at(j);
              gapExtendP +=+gapExB;

              //Score of Open Gap in B
              gapOpenQ = D->at(i).at(j-1);
              gapOpenQ += scoreM->ScoreOfGapOpen(A->getChar(i-1));
              gapOpenQ += gapExA;

              //Score of continue Gap in B
              gapExtendQ = Q->at(i).at(j-1);
              gapExtendQ += gapExA;

              //Set Helper Matrices Values
              P->at(i).at(j) = scoreM->BestOfTwo(gapOpenP, gapExtendP);
              Q->at(i).at(j) = scoreM->BestOfTwo(gapOpenQ, gapExtendQ);

              //Set Distance Matrix Value
              D->at(i).at(j) = scoreM->BestOfThree(diagScore, P->at(i).at(j), Q->at(i).at(j) );

              //Check for 0
              if (D->at(i).at(j) < 0)
                D->at(i).at(j) = 0;

              //Save best Score
              if (scoreM->BestOfTwo(bestScore, D->at(i).at(j)) != bestScore)
                {
                  bestScore = D->at(i).at(j);
                  bi = i;
                  bj = j;
                }

            }
        return bestScore;
      }

    Alignment AlignmentFactory::BacktrackingSmithWatermanAffin(
        const bpp::Sequence* A, const bpp::Sequence* B, ScoringModel* scoreM,
        const ScoreMatrix* D, const ScoreMatrix* P, const ScoreMatrix* Q,
        uint &i, uint &j)
      {
        Alignment out(A->getAlphabet());

        std::vector<int> outA;
        std::vector<int>::iterator itA = outA.begin();
        std::vector<int> outB;
        std::vector<int>::iterator itB = outB.begin();

        double minScore = 0;

        while ( (i>0) && (j>0))
          {
            //                        std::cout << "Aktuelle Score: " << minScore << std::endl;

            if (D->at(i).at(j) == 0)
              break;

            //Change to Q
            if (scoreM->BestOfTwo(D->at(i).at(j), Q->at(i).at(j) ) == Q->at(i).at(j))
              {

                //Lonely Gap in A
                if (Q->at(i).at(j) == D->at(i).at(j-1) + scoreM->ScoreOfGapOpen(A->getChar(i -1))
                    + scoreM->ScoreOfGapExtend(A->getChar(i-1)))
                  {
                    minScore +=scoreM->ScoreOfGapOpen(A->getChar(i-1));
                    //                                        std::cout << "Left"<< std::endl;
                    //                    std::cout << "Score + " << scoreM->ScoreOfGapOpen(A->getChar(i-1))
                    //                    << std::endl;
                  }
                //Continued Gap in A

                else if (Q->at(i).at(j) == Q->at(i).at(j-1) + scoreM->ScoreOfGapExtend(A->getChar(i-1)))
                  {
                    //                                        std::cout << "Left cont "<< std::endl;
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
                    //                                        std::cout << "Up"<< std::endl;
                    minScore +=scoreM->ScoreOfGapOpen(B->getChar(j-1));
                    //                    std::cout << "Score + " << scoreM->ScoreOfGapOpen(B->getChar(j-1))
                    //                    << std::endl;
                  }
                //Continued Gap in B

                else if (P->at(i).at(j) == P->at(i-1).at(j) + scoreM->ScoreOfGapExtend(B->getChar(j -1)))
                  {
                    //                                        std::cout << "Up cont"<< std::endl;
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
                //                                                std::cout << "Diag"<< std::endl;
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

        bpp::Sequence* seqA = new bpp::Sequence(A->getName(), outA, A->getAlphabet());
        bpp::Sequence* seqB = new bpp::Sequence(B->getName(), outB, B->getAlphabet());

        //Save Score to Container
        out.set_Score(minScore);

        //Add Alinged Sequences to Container
        out.addSequence(seqA);
        out.addSequence(seqB);

        delete seqA;
        delete seqB;

        return out;
      }

    Alignment AlignmentFactory::NeedlemanWunschAlignment(
        const bpp::Sequence* inA, const bpp::Sequence* inB,
        ScoringModel* scoreM, bool verbose)
      {
        ScoreMatrix D =
            matrix->InitializeScoreMatrixDistances(inA, inB, scoreM);

        //Forward Iteration
        ForwardRecursionNMW(inA, inB, scoreM, &D);

        uint i = D.size()-1;
        uint j = D.at(0).size()-1;

        return BacktrackingNMW(inA, inB, scoreM, &D, i, j);

      }

    Alignment AlignmentFactory::GotohAlignment(const bpp::Sequence* inA,
        const bpp::Sequence* inB, ScoringModel* scoreM, bool verbose)
      {

        ScoreMatrix D = matrix->InitScoreMatrixWith(inA, inB, 0);
        ScoreMatrix P = matrix->InitScoreMatrixWith(inA, inB, 0);
        ScoreMatrix Q = matrix->InitScoreMatrixWith(inA, inB, 0);

        //Forward Iteration
        ForwardRecursionGotoh(inA, inB, scoreM, &D, &P, &Q);

        uint i = D.size()-1;
        uint j = D.at(0).size()-1;

        return BacktrackingGotohGlobal(inA, inB, scoreM, &D, &P, &Q, i, j);

      }
    Alignment AlignmentFactory::SmithWaterman(const bpp::Sequence* inA,
        const bpp::Sequence* inB, ScoringModel* scoreM, bool verbose)
      {
        ScoreMatrix D =
            matrix->InitializeScoreMatrixDistances(inA, inB, scoreM);

        uint i = 0;
        uint j = 0;

        ForwardRecursionSmithWaterman(inA, inB, scoreM, &D, i, j);

        return BacktrackingSmithWaterman(inA, inB, scoreM, &D, i, j);

      }
    Alignment AlignmentFactory::SmithWatermanAffin(const bpp::Sequence* inA,
        const bpp::Sequence* inB, ScoringModel* scoreM, bool verbose)
      {
        ScoreMatrix D = matrix->InitScoreMatrixWith(inA, inB, 0);
        ScoreMatrix P = matrix->InitScoreMatrixWith(inA, inB, 0);
        ScoreMatrix Q = matrix->InitScoreMatrixWith(inA, inB, 0);

        uint i = 0;
        uint j = 0;

        ForwardRecursionSmithWatermanAffin(inA, inB, scoreM, &D, &P, &Q, i, j);

        return BacktrackingSmithWatermanAffin(inA, inB, scoreM, &D, &P, &Q, i,
            j);

      }

  }
