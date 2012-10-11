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

#ifndef WHITESPACEFASTA_H_
#define WHITESPACEFASTA_H_

#include <sstream>
#include <Bpp/Seq/Io/Fasta.h>


namespace Circal
  {

    class WhitespaceFasta : public bpp::Fasta
      {
  public:
      WhitespaceFasta();
      virtual ~WhitespaceFasta();
      void
          appendFromStream(istream & input, bpp::VectorSequenceContainer & vsc) const
              throw (bpp::Exception);
      };

  }

#endif /*WHITESPACEFASTA_H_*/
