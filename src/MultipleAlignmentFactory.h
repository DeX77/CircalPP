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

#ifndef MULTIPLEALIGNMENT_H_
#define MULTIPLEALIGNMENT_H_

#include "AlignmentFactory.h"
#include "Alignment.h"

namespace Circal
  {

    class MultipleAlignmentFactory : public virtual AlignmentFactory
      {

  public:
      explicit MultipleAlignmentFactory();
      virtual ~MultipleAlignmentFactory();

      void
          GotohalignMultiple(Alignment* out,
              const bpp::VectorSequenceContainer* input,
              const ScoringModel* scoreM);
      void
          NMWalignMultiple(Alignment* out,
              const bpp::VectorSequenceContainer* input,
              const ScoringModel* scoreM);

      };
  }
#endif /*MULTIPLEALIGNMENT_H_*/
