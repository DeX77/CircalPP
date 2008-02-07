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

#ifndef MULTIPLECIRCULARALIGNMENT_H_
#define MULTIPLECIRCULARALIGNMENT_H_

#include "CircularAlignment.h"
#include "MultipleAlignment.h"

namespace Circal
  {
    class MultipleCircularAlignment : public CircularAlignment,
      public MultipleAlignment
      {
    public:
        explicit MultipleCircularAlignment(const bpp::Alphabet* alpha);
        MultipleCircularAlignment(const VectorSequenceContainer* input,
            const ScoringModel* scoreM);
        virtual ~MultipleCircularAlignment();
      };
  }

#endif /*MULTIPLECIRCULARALIGNMENT_H_*/
