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

#include "MultipleAlignmentFactory.h"
#include "Output.h"
#include "MatrixHelper.h"

namespace Circal
  {
    MultipleAlignmentFactory::MultipleAlignmentFactory()
      {
      }

    MultipleAlignmentFactory::~MultipleAlignmentFactory()
      {
      }

    Alignment MultipleAlignmentFactory::GotohalignMultiple(
        const bpp::VectorSequenceContainer* input, ScoringModel* scoreM,
        int &delta, bool verbose)
      {

        Alignment out(input->getAlphabet());

#ifdef _OPENMP            
#pragma omp parallel for
#endif    
        for (uint u=0; u<input->getNumberOfSequences(); u++)
          {

            for (uint k=u+1; k<input->getNumberOfSequences(); k++)
              {
                if (verbose)
                  std::cout << "*" << std::flush;
                Alignment temp = GotohAlignment(input->getSequence(u),
                    input->getSequence(k), scoreM, verbose);
                out.addSequence(temp.getSequence(0));
                out.addSequence(temp.getSequence(1));
              }
            if (verbose)
              std::cout << std::endl;

          }
        return out;

      }
    Alignment MultipleAlignmentFactory::NMWalignMultiple(
        const bpp::VectorSequenceContainer* input, ScoringModel* scoreM,
        int &delta, bool verbose)
      {
        Alignment out(input->getAlphabet());

#ifdef _OPENMP            
#pragma omp parallel for
#endif    
        for (uint u=0; u<input->getNumberOfSequences(); u++)
          {

            for (uint k=u+1; k<input->getNumberOfSequences(); k++)
              {
                if (verbose)
                  std::cout << "*" << std::flush;
                Alignment temp = NeedlemanWunschAlignment(
                    input->getSequence(u), input->getSequence(k), scoreM);

                out.addSequence(temp.getSequence(0));
                out.addSequence(temp.getSequence(1));
              }
            if (verbose)
              std::cout << std::endl;

          }
        return out;

      }

  }
