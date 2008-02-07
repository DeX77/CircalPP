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

#include "MultipleAlignment.h"
#include "Output.h"
#include "MatrixHelper.h"

namespace Circal
  {
    MultipleAlignment::MultipleAlignment(const bpp::Alphabet* alpha) :
      VectorSequenceContainer(alpha), Alignment(alpha)
      {
      }

    MultipleAlignment::~MultipleAlignment()
      {
      }

    MultipleAlignment::MultipleAlignment(const VectorSequenceContainer* input,
        const ScoringModel* scoreM) :
      VectorSequenceContainer(input->getAlphabet()),
          Alignment(input->getAlphabet())
      {
      }

    void MultipleAlignment::GotohalignMultiple(
        const VectorSequenceContainer* input, const ScoringModel* scoreM)
      {
        set_origSize(input->getNumberOfSequences());

        for (uint u=0; u<input->getNumberOfSequences(); u++)
          {

            for (uint k=u+1; k<input->getNumberOfSequences(); k++)
              {

                Alignment* temp = GotohAlignment(input->getSequence(u),
                    input->getSequence(k), scoreM);
                std:cerr << prettyPrint->AlignmentPrettyPrint(temp);

                bpp::SequenceContainerTools::append(*this, *temp);
              }

          }

      }
    void MultipleAlignment::NMWalignMultiple(
        const VectorSequenceContainer* input, const ScoringModel* scoreM)
      {

        std::cout << input->getNumberOfSequences()
            *input->getNumberOfSequences() << " to go" << std::endl;

          {

            for (uint u=1; u<input->getNumberOfSequences(); u++)
              {
                for (uint k=u+1; k<input->getNumberOfSequences(); k++)
                  {

                    Alignment* temp = NeedlemanWunschAlignment(
                        input->getSequence(u), input->getSequence(k), scoreM);
                    
                    prettyPrint->AlignmentPrettyPrint(temp);

                    bpp::SequenceContainerTools::append(*this, *temp);
                    std::cout
                        << "******************************** Done with Nr. "
                        << u << "," << k << std::endl;
                  }

              }
          }
      }

    MultipleAlignment* MultipleAlignment::TopologicSortedAlignment(
        BoolMatrix G, const Alignment* pairWiseAlignments)
      {

        MultipleAlignment* out = new MultipleAlignment(pairWiseAlignments->getAlphabet());

        //Search connected Compounds.. those are represented by 1 in the same Column of the Matrix and not divided by 0

          {
            for (uint i=0; i<G.size(); i++)
              {

                for (uint j=0; j< G[i].size(); j++)
                  {
                    //Check for match or not mismatch := !gap
                    if (pairWiseAlignments->getSequence(i)->getValue(j) != -1)
                      G[i/2][j] = true;
                    else
                      G[i/2][j] = false;
                  }
              }
          }
        return out;
      }
  }
