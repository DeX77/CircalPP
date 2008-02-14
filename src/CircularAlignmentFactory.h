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

#ifndef CIRCULARALIGNMENT_H_
#define CIRCULARALIGNMENT_H_

#include "AlignmentFactory.h"
namespace Circal
  {

    class CircularAlignmentFactory : public virtual AlignmentFactory
      {

  public:
      explicit CircularAlignmentFactory();
      virtual ~CircularAlignmentFactory();

      virtual Alignment* NeedlemanWunschAlignment(const bpp::Sequence* A,
          const bpp::Sequence* B, const ScoringModel* scoreM);

      virtual Alignment* GotohAlignment(const bpp::Sequence* A,
          const bpp::Sequence* B, const ScoringModel* scoreM);
      };
  }
#endif /*CIRCULARALIGNMENT_H_*/
