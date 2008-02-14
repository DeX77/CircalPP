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

#include "Output.h"
#include "MatrixHelper.h"
#include "Alignment.h"
#include "VertebrateMitochondrialGenomeAlphabet.h"
#include "PseudoRotatedSequence.h"
#include <Seq/SymbolList.h>

#include <Seq/sequences>
#include <iostream>
#include <sstream>

namespace Circal
  {
    Output::Output()
      {
      }

    Output::~Output()
      {
        
      }

    std::string Output::ScoreMatrixPrettyPrint(const bpp::Sequence* A,
        const bpp::Sequence* B, const ScoreMatrix &D)
      {
        std::stringstream out;

        out << "\t";
        for (uint k=0; k<B->size(); k++)
          out << "\t"<< B->getChar(k);
        out << std::endl;
        for (uint i=0; i<D.size(); i++)
          {
            if (i != 0)
              out << A->getChar(i-1);
            out << "\t";
            for (uint j=0; j<D.at(i).size(); j++)
              out << D.at(i).at(j)<< "\t";
            out << std::endl;
          }
        return out.str();
      }

    std::string Output::AlignmentPrettyPrint(const Alignment* aln)
      {
        std::stringstream out;
        for (uint i=0; i<aln->getNumberOfSequences(); i++)
          out
              << SequencePrettyPrint(const_cast<bpp::Sequence*>(aln->getSequence(i)))
              << std::endl;
        if (aln->getNumberOfSequences())
          {
            out << "---------------------------------------------------"
                << std::endl;
            out << "Score: " << aln->get_Score() << std::endl;
          }
        return out.str();
      }

    std::string Output::AdjacenceMatrixPrettyPrint(const BoolMatrix &G)
      {
        std::stringstream out;
        for (uint i=0; i<G.size(); i++)
          {
            out << "Aln"<< i << "\t\t"<< std::flush;
            for (uint j=0; j<G[0].size(); j++)
              out << " "<< G[i][j]<< std::flush;
            out << std::endl;
          }
        return out.str();
      }

    std::string Output::SequencePrettyPrint(const bpp::Sequence* A)
      {
        std::stringstream out;
        out << "Sequenz " << A->getName() << " Laenge: " << A->size() << " :" << std::endl;

        for (uint i=0; i< A->size(); i++)
          out << " "<< A->getChar(i);
        return out.str();

      }
    //Shamelessly stolen from MARNA 
    std::string Output::TCoffeeLibFormat(const Alignment* aln)
      {
        std::stringstream out;
        int p1 = 1;
        int p2 = 1;

        int index1 = 1;
        int index2 = 2;
        out << aln->get_origSize() << endl;


        out << aln->getSequence(0)->getName() << " ";
        out << aln->getSequence(0)->size() << " ";
        out << PseudoRotatedSequence(aln->getSequence(0)).toString();
        out << std::endl;

        //Followed by a line with SeqName, Nr of Sited, Sequence for each Sequence        
        for (uint i=1; i < ((aln->get_origSize()-1)*2); i+=2)
          {
            out << aln->getSequence(i)->getName() << " ";
            out << aln->getSequence(i)->size() << " ";
            out << PseudoRotatedSequence(aln->getSequence(i)).toString();
            out << std::endl;
          }

        //As we save pairs of Sequences we have to iterate by 2
        for (uint j=0; j < aln->getNumberOfSequences(); j+=2)
          {
            p1 = 1;
            p2 = 1;

            if (index2 == aln->get_origSize())
              {
                index1++;
                index2 = index1+1;
              }
            out << "#";
            out << index1;
            out << " ";
            out << index2;
            out << std::endl;
            index2++;

            for (uint i=0; i < aln->getSequence(j)->size(); i++)
              {

                if (aln->getSequence(j)->getValue(i) != -1)
                  ++p1;
                if (aln->getSequence(j+1)->getValue(i) != -1)
                  ++p2;
                if ( (aln->getSequence(j)->getValue(i) != -1) && (aln->getSequence(j+1)->getValue(i) != -1))
                  //TODO Weight is hardcoded here
                  out << p1-1 << " " << p2-1 << " " << 100 << std::endl;
              }
          }
        //Tail of T_Coffee_Lib format
        out << "CPU 0" << std::endl;
        out << "! SEQ_1_TO_N" << std::endl;

        return out.str();
      }

  }
