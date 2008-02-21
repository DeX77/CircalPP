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
    void PseudoCircularAlignmentFactory::ForwardRecursionGotoh(
        const bpp::Sequence* A, const PseudoRotatedSequence* B,
        const ScoringModel* scoreM, const int &delta, ScoreMatrix3D* D,
        ScoreMatrix3D* P, ScoreMatrix3D* Q)
      {

        double gapOpenP;
        double gapExtendP;
        double gapOpenQ;
        double gapExtendQ;
        double diagScore;
        double maxScore = 0;
        double maxScoreJ = 0;

        //B is doubled=pseudorotated so correct size is B/2
        int slaps = B->size()/2;
        //catch (stupid User Exception)
        if (delta != 0)
          slaps /= delta;

        //Forward
        //#ifdef _OPENMP            
        //#pragma omp parallel for
        //#endif           
        for (uint i=1; i<A->size()+1; i++)
          for (uint j=1; j<B->size()+1; j++)
            {
              if (j % delta == 1)
                {

                  //Score of Match eg. Mismatch
                  diagScore = D->at(i-1).at(j-1).at(0) + scoreM->ScoreOf(A->getChar(i-1),
                      B->getChar(j -1));

                  //Score of Open Gap in A
                  gapOpenP = D->at(i-1).at(j).at(0) + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                      + scoreM->ScoreOfGapExtend(B->getChar(j-1));
                  //Score of continue Gap in A
                  gapExtendP = P->at(i-1).at(j).at(0)+ scoreM->ScoreOfGapExtend(B->getChar(j-1));

                  //Score of Open Gap in B
                  gapOpenQ = D->at(i).at(j-1).at(0) + scoreM->ScoreOfGapOpen(A->getChar(i-1))
                      + scoreM->ScoreOfGapExtend(A->getChar(i-1));
                  //Score of continue Gap in B
                  gapExtendQ = Q->at(i).at(j-1).at(0) + scoreM->ScoreOfGapExtend(A->getChar(i-1));

                  //Set Helper Matrices Values
                  P->at(i).at(j).at(0) = scoreM->BestOfTwo(gapOpenP, gapExtendP);
                  Q->at(i).at(j).at(0) = scoreM->BestOfTwo(gapOpenQ, gapExtendQ);

                  //Set Distance Matrix Value
                  D->at(i).at(j).at(0) = scoreM->BestOfThree(diagScore, P->at(i).at(j).at(0), Q->at(i).at(j).at(0) );

                  for (int k=1; k<slaps; k++)
                    {
                      //Score of Match eg. Mismatch
                      diagScore = D->at(i-1).at(j-1).at(k-1) + scoreM->ScoreOf(A->getChar(i-1),
                          B->getChar(j -1));

                      //Score of Open Gap in A
                      gapOpenP = D->at(i-1).at(j).at(k) + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                          + scoreM->ScoreOfGapExtend(B->getChar(j-1));
                      //Score of continue Gap in A
                      gapExtendP = P->at(i-1).at(j).at(k)+ scoreM->ScoreOfGapExtend(B->getChar(j-1));

                      //Score of Open Gap in B
                      gapOpenQ = D->at(i).at(j-1).at(k-1)
                          + scoreM->ScoreOfGapOpen(A->getChar(i-1))
                          + scoreM->ScoreOfGapExtend(A->getChar(i-1));
                      //Score of continue Gap in B
                      gapExtendQ = Q->at(i).at(j-1).at(k-1) + scoreM->ScoreOfGapExtend(A->getChar(i
                          -1));

                      //Set Helper Matrices Values
                      P->at(i).at(j).at(k) = scoreM->BestOfTwo(gapOpenP, gapExtendP);
                      Q->at(i).at(j).at(k) = scoreM->BestOfTwo(gapOpenQ, gapExtendQ);

                      //Set Distance Matrix Value
                      D->at(i).at(j).at(k) = scoreM->BestOfThree(diagScore, P->at(i).at(j).at(k), Q->at(i).at(j).at(k) );

                      //Save actual best Score
                      if (scoreM->BestOfTwo(maxScore, D->at(i).at(j).at(k)) != maxScore)
                        {
                          maxScore = D->at(i).at(j).at(k);
                          maxScoreJ = j;
                        }
                    }
                }
              else
                for (int k=0; k<slaps; k++)
                  {

                    //Score of Match eg. Mismatch
                    diagScore = D->at(i-1).at(j-1).at(k) + scoreM->ScoreOf(A->getChar(i-1),
                        B->getChar(j -1));

                    //Score of Open Gap in A
                    gapOpenP = D->at(i-1).at(j).at(k) + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                        + scoreM->ScoreOfGapExtend(B->getChar(j-1));
                    //Score of continue Gap in A
                    gapExtendP = P->at(i-1).at(j).at(k)+ scoreM->ScoreOfGapExtend(B->getChar(j-1));

                    //Score of Open Gap in B
                    gapOpenQ = D->at(i).at(j-1).at(k) + scoreM->ScoreOfGapOpen(A->getChar(i-1))
                        + scoreM->ScoreOfGapExtend(A->getChar(i-1));
                    //Score of continue Gap in B
                    gapExtendQ = Q->at(i).at(j-1).at(k) + scoreM->ScoreOfGapExtend(A->getChar(i-1));

                    //Set Helper Matrices Values
                    P->at(i).at(j).at(k) = scoreM->BestOfTwo(gapOpenP, gapExtendP);
                    Q->at(i).at(j).at(k) = scoreM->BestOfTwo(gapOpenQ, gapExtendQ);

                    //Set Distance Matrix Value
                    D->at(i).at(j).at(k) = scoreM->BestOfThree(diagScore, P->at(i).at(j).at(k), Q->at(i).at(j).at(k) );

                    //Save actual best Score
                    if (scoreM->BestOfTwo(maxScore, D->at(i).at(j).at(k)) != maxScore)
                      {
                        maxScore = D->at(i).at(j).at(k);
                        maxScoreJ = j;
                      }

                  }

            }
        cout << "Actual best Score: " << maxScore << "@ j=" << maxScoreJ
            << std::endl;

      }
    Alignment PseudoCircularAlignmentFactory::BacktrackingGotohLocal(
        const bpp::Sequence* A, const PseudoRotatedSequence* B,
        const ScoringModel* scoreM, const int &delta, const ScoreMatrix3D* D,
        const ScoreMatrix3D* P, const ScoreMatrix3D* Q, uint &i, uint &j)
      {

        //Debug purpose
        stringstream cout;

        Alignment out(A->getAlphabet());

        vector<int> outA;
        vector<int>::iterator itA = outA.begin();
        vector<int> outB;
        vector<int>::iterator itB = outB.begin();

        double minScore = 0;

        //B is doubled=pseudorotated so correct size is B/2
        int slaps = B->size()/2;

        if (delta != 0)
          slaps /= delta;

        while ( (i>0) && (j>0))
          {

            int gk = 0;
            double bestScore = 0;

            //Search maximum score within k
            for (int k=0; k<slaps; k++)
              {
                //Search in D
                if (scoreM->BestOfTwo(D->at(i).at(j).at(k), bestScore) != bestScore)
                  {
                    bestScore = D->at(i).at(j).at(k);
                    gk = k;
                  }

                //.. or Q
                if (scoreM->BestOfTwo(Q->at(i).at(j).at(k), bestScore) != bestScore)
                  {
                    bestScore = Q->at(i).at(j).at(k);
                    gk = k;
                  }
                //.. or P
                if (scoreM->BestOfTwo(P->at(i).at(j).at(k), bestScore) != bestScore)
                  {
                    bestScore = P->at(i).at(j).at(k);
                    gk = k;
                  }

              }

            //Change to Q
            if (scoreM->BestOfTwo(D->at(i).at(j).at(gk), Q->at(i).at(j).at(gk) ) == Q->at(i).at(j).at(gk))
              {

                //Lonely Gap in A
                if (Q->at(i).at(j).at(gk) == D->at(i).at(j-1).at(gk) + scoreM->ScoreOfGapOpen(A->getChar(i -1))
                    + scoreM->ScoreOfGapExtend(A->getChar(i-1)))
                  {
                    minScore +=scoreM->ScoreOfGapOpen(A->getChar(i-1));
                    cout << "Left"<< endl;
                    cout << "Score + "
                        << scoreM->ScoreOfGapOpen(A->getChar(i-1)) << std::endl;
                  }
                //Continued Gap in A
                else if (Q->at(i).at(j).at(gk) == Q->at(i).at(j-1).at(gk) + scoreM->ScoreOfGapExtend(A->getChar(i-1)))
                  {
                    cout << "Left cont "<< endl;
                  }
                else
                  {
                    cout << "Changed to Q but not possible?" << std::endl;
                  }
                itA = outA.insert(itA, -1);
                itB = outB.insert(itB, B->getValue(j-1));
                j--;
                minScore +=scoreM->ScoreOfGapExtend(A->getChar(i-1));
                cout << "Score + " << scoreM->ScoreOfGapExtend(A->getChar(i-1))
                    << std::endl;
                continue;
              }

            //Change to P
            if (scoreM->BestOfTwo(D->at(i).at(j).at(gk) , P->at(i).at(j).at(gk) ) == P->at(i).at(j).at(gk))
              {
                //Lonely Gap in B
                if (P->at(i).at(j).at(gk) == D->at(i-1).at(j).at(gk) + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                    + scoreM->ScoreOfGapExtend(B->getChar(j-1)))
                  {
                    cout << "Up"<< endl;
                    minScore +=scoreM->ScoreOfGapOpen(B->getChar(j-1));
                    cout << "Score + "
                        << scoreM->ScoreOfGapOpen(B->getChar(j-1)) << std::endl;
                  }
                //Continued Gap in B
                else if (P->at(i).at(j).at(gk) == P->at(i-1).at(j).at(gk) + scoreM->ScoreOfGapExtend(B->getChar(j -1)))
                  {
                    cout << "Up cont"<< endl;
                  }
                else
                  {
                    //                    cout << "Changed to P but not possible?" << std::endl;
                  }
                minScore +=scoreM->ScoreOfGapExtend(B->getChar(j-1));
                cout << "Score + " << scoreM->ScoreOfGapExtend(B->getChar(j-1))
                    << std::endl;
                itA = outA.insert(itA, A->getValue(i-1));
                itB = outB.insert(itB, -1);
                i--;
                continue;

              }

            //Diagonal
            else if (D->at(i).at(j).at(gk) == D->at(i-1).at(j-1).at(gk) + scoreM->ScoreOf(A->getChar(i-1), B->getChar(j-1)))
              {
                minScore += scoreM->ScoreOf(A->getChar(i-1), B->getChar(j-1));
                cout << "Diag"<< endl;
                cout << "Score + " << scoreM->ScoreOf(A->getChar(i-1),
                    B->getChar(j-1)) << std::endl;
                itA = outA.insert(itA, A->getValue(i-1));
                itB = outB.insert(itB, B->getValue(j-1));
                i--;
                j--;
                continue;
              }

            else
              {

                cout << "42!!!"<< std::endl;
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
    void PseudoCircularAlignmentFactory::ForwardRecursionNMW(
        const bpp::Sequence* A, const PseudoRotatedSequence* B,
        const ScoringModel* scoreM, const int &delta, ScoreMatrix3D* D)
      {
        double gapOpenP;
        double gapOpenQ;
        double diagScore;

        //B is doubled=pseudorotated so correct size is B/2
        int slaps = B->size()/2;

        if (delta != 0)
          slaps /= delta;

        //Forward
        for (uint i=1; i<A->size()+1; i++)
          for (uint j=1; j<B->size()+1; j++)
            {
              if (j % delta ==1)
                {
                  //Score of Match eg. Mismatch
                  diagScore = D->at(i-1).at(j-1).at(0) + scoreM->ScoreOf(A->getChar(i-1),
                      B->getChar(j -1));

                  //Score of Open Gap in A
                  gapOpenP = D->at(i-1).at(j).at(0) + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                      + scoreM->ScoreOfGapExtend(B->getChar(j-1));

                  //Score of Open Gap in B
                  gapOpenQ = D->at(i).at(j-1).at(0) + scoreM->ScoreOfGapOpen(A->getChar(i-1))
                      + scoreM->ScoreOfGapExtend(A->getChar(i-1));

                  //Set Distance Matrix Value
                  D->at(i).at(j).at(0) = scoreM->BestOfThree(diagScore, gapOpenP,
                      gapOpenQ);

                  for (int k=1; k<slaps; k++)
                    {
                      //Score of Match eg. Mismatch
                      //!! Changed from Original
                      diagScore = D->at(i-1).at(j-1).at(k-1) + scoreM->ScoreOf(A->getChar(i-1),
                          B->getChar(j -1));

                      //Score of Open Gap in A
                      gapOpenP = D->at(i-1).at(j).at(k) + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                          + scoreM->ScoreOfGapExtend(B->getChar(j-1));

                      //Score of Open Gap in B
                      //!! Changed from Original
                      gapOpenQ = D->at(i).at(j-1).at(k-1)
                          + scoreM->ScoreOfGapOpen(A->getChar(i-1))
                          + scoreM->ScoreOfGapExtend(A->getChar(i-1));

                      //Set Distance Matrix Value
                      D->at(i).at(j).at(k) = scoreM->BestOfThree(diagScore, gapOpenP,
                          gapOpenQ);

                    }

                }
              else
                for (int k=0; k<slaps; k++)
                  {
                    //Score of Match eg. Mismatch
                    diagScore = D->at(i-1).at(j-1).at(k) + scoreM->ScoreOf(A->getChar(i-1),
                        B->getChar(j -1));

                    //Score of Open Gap in A
                    gapOpenP = D->at(i-1).at(j).at(k) + scoreM->ScoreOfGapOpen(B->getChar(j-1))
                        + scoreM->ScoreOfGapExtend(B->getChar(j-1));

                    //Score of Open Gap in B
                    gapOpenQ = D->at(i).at(j-1).at(k) + scoreM->ScoreOfGapOpen(A->getChar(i-1))
                        + scoreM->ScoreOfGapExtend(A->getChar(i-1));

                    //Set Distance Matrix Value
                    D->at(i).at(j).at(k) = scoreM->BestOfThree(diagScore, gapOpenP,
                        gapOpenQ);

                  }

            }
      }
    Alignment PseudoCircularAlignmentFactory::BacktrackingNMW(
        const bpp::Sequence* A, const PseudoRotatedSequence* B,
        const ScoringModel* scoreM, const int &delta, const ScoreMatrix3D* D,
        uint &i, uint &j)
      {

        vector<int> outA;
        vector<int>::iterator itA = outA.begin();
        vector<int> outB;
        vector<int>::iterator itB = outB.begin();

        Alignment out(A->getAlphabet());

        double minScore = double(0);

        //B is doubled=pseudorotated so correct size is B/2
        int slaps = B->size()/2;
        if (delta != 0)
          slaps /= delta;

        double ScoreD;
        double ScoreDiag;
        double ScoreUpD;
        double ScoreLeftD;

        while ( (i>0) && (j>0))
          {
            int k = 0;
            ScoreD = D->at(i).at(j).at(k);

            //Search maximum score within k
            for (k=1; k<slaps; k++)
              {
                if (scoreM->BestOfTwo(D->at(i).at(j).at(k), ScoreD) != ScoreD)
                  ScoreD = D->at(i).at(j).at(k);
              }

            ScoreDiag = D->at(i-1).at(j-1).at(k);
            ScoreLeftD = D->at(i).at(j-1).at(k);
            ScoreUpD = D->at(i-1).at(j).at(k);

            //Lonely Gap in A
            if (ScoreD == ScoreLeftD + scoreM->ScoreOfGapOpen(A->getChar(i -1)))
              {
                //cout << "Left"<< endl;
                minScore += scoreM->ScoreOfGapOpen(A->getChar(i-1));
                itA = outA.insert(itA, -1);
                itB = outB.insert(itB, B->getValue(j-1));
                j--;
              }

            //Lonely Gap in B
            else if (ScoreD == ScoreUpD
                + scoreM->ScoreOfGapOpen(B->getChar(j-1)))
              {
                //cout << "Up"<< endl;
                minScore += scoreM->ScoreOfGapOpen(B->getChar(j-1));
                itA = outA.insert(itA, A->getValue(i-1));
                itB = outB.insert(itB, -1);
                i--;
              }
            //Diagonal
            else if (ScoreD == ScoreDiag + scoreM->ScoreOf(A->getChar(i-1),
                B->getChar(j-1)))
              {
                minScore += scoreM->ScoreOf(A->getChar(i-1), B->getChar(j-1));
                //cout << "Diag"<< endl;
                itA = outA.insert(itA, A->getValue(i-1));
                itB = outB.insert(itB, B->getValue(j-1));
                i--;
                j--;
              }

            else
              {
                cout << "42!!! Lost in D"<< endl;
              }

          }

        while (i > 0)
          {
            //cout << "Overlapp A"<< endl;
            minScore += scoreM->ScoreOfGapOpen(A->getChar(i-1));
            itA = outA.insert(itA, A->getValue(i-1));
            itB = outB.insert(itB, -1);
            i--;
          }

        while (j > 0)
          {
            //cout << "Overlapp B"<< endl;
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

    Alignment PseudoCircularAlignmentFactory::NeedlemanWunschAlignment(
        const bpp::Sequence* inA, const bpp::Sequence* inB,
        const ScoringModel* scoreM, const int &delta)
      {
        PseudoRotatedSequence* B = new PseudoRotatedSequence(inB);

        ScoreMatrix D = matrix->InitializeScoreMatrixDistances(inA, B, scoreM);

        //Forward Iteration
        AlignmentFactory::ForwardRecursionNMW(inA, B, scoreM, &D);

        //Search Starting Node
        uint i = D.size()-1;
        uint j = D.at(0).size()-1;
        i = matrix->SearchBestInColumn(&D, scoreM, i, j);

        uint horizontalStart = i;
        Alignment temp = AlignmentFactory::BacktrackingNMW(inA, B, scoreM, &D,
            i, j);
        uint horizontalEnd = i;

        //Check for length of actual Alignment
        if ((horizontalStart-horizontalEnd) < inB->size())
          {
            //Seems to be ok
            return temp;
          }
        ScoreMatrix3D kD = matrix->InitScoreMatrix3DWith(inA, B, delta, 0);

        ForwardRecursionNMW(inA, B, scoreM, delta, &kD);

        return BacktrackingNMW(inA, B, scoreM, delta, &kD, i, j);

      }

    Alignment PseudoCircularAlignmentFactory::GotohAlignment(
        const bpp::Sequence* inA, const bpp::Sequence* inB,
        const ScoringModel* scoreM, const int &delta)
      {
        PseudoRotatedSequence B(inB);

        ScoreMatrix D = matrix->InitScoreMatrixWith(inA, &B, 0);
        ScoreMatrix P = matrix->InitScoreMatrixWith(inA, &B, 0);
        ScoreMatrix Q = matrix->InitScoreMatrixWith(inA, &B, 0);

        //Forward Iteration
        AlignmentFactory::ForwardRecursionGotoh(inA, &B, scoreM, &D, &P, &Q);

        //Search Starting Node
        uint i = D.size()-1;
        uint j = D.at(0).size()-1;

        //j = matrix->SearchBestInColumn(&D, scoreM, i, j);
        uint horizontalStart = j;
        Alignment temp = AlignmentFactory::BacktrackingGotohGlocal(inA, &B,
            scoreM, &D, &P, &Q, i, j);
        uint horizontalEnd = j;

        //Check for length of actual Alignment
        if ((horizontalStart-horizontalEnd) <= inB->size())
          {
            //Seems to be ok
            std::cout << "Optimum sofort gefunden!" << std::endl;
            return temp;
          }

        //The hard way
        ScoreMatrix3D kD = matrix->InitScoreMatrix3DWith(inA, &B, delta, 0);
        ScoreMatrix3D kP = matrix->InitScoreMatrix3DWith(inA, &B, delta, 0);
        ScoreMatrix3D kQ = matrix->InitScoreMatrix3DWith(inA, &B, delta, 0);

        //Forward Iteration
        ForwardRecursionGotoh(inA, &B, scoreM, delta, &kD, &kP, &kQ);

        i = kD.size()-1;
        j = kD.at(0).size()-1;
        uint k= kD.at(0).at(0).size()-1;

        j = matrix->SearchBestInColumn3D(&kD, scoreM, i, j, k);
        std::cout << "Best Score: " << kD.at(i).at(j).at(k) << " gefunden in j=" << j << std::endl;

        return BacktrackingGotohLocal(inA, &B, scoreM, delta, &kD, &kP, &kQ, i,
            j);

      }
  }
