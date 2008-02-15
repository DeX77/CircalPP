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

#ifndef PSEUDOCIRCULARALIGNMENT_H_
#define PSEUDOCIRCULARALIGNMENT_H_

#include "AlignmentFactory.h"

namespace Circal
  {
    class PseudoCircularAlignmentFactory : public virtual AlignmentFactory
      {

  public:
      PseudoCircularAlignmentFactory();
      virtual ~PseudoCircularAlignmentFactory();

      void NeedlemanWunschAlignment(Alignment* out, const bpp::Sequence* A,
          const bpp::Sequence* B, const ScoringModel* scoreM);

      void GotohAlignment(Alignment* out, const bpp::Sequence* A,
          const bpp::Sequence* B, const ScoringModel* scoreM);

      };
  }

#endif /*PSEUDOCIRCULARALIGNMENT_H_*/
