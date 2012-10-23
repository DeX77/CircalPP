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

#include "SequenceProxy.h"

namespace Circal
  {
    SequenceProxy::SequenceProxy(const std::string& name,
        const std::vector<int>& sequence, const bpp::Alphabet* alpha)
            throw (bpp::BadIntException) :
        bpp::BasicSequence(name, sequence, alpha)
      {
        // TODO Auto-generated constructor stub

      }
    SequenceProxy::SequenceProxy(const bpp::Sequence& sequence) :
        bpp::BasicSequence(sequence)
      {

      }

    SequenceProxy::SequenceProxy(const std::string& name,
        const std::string& sequence, const bpp::Alphabet* alpha) :
        bpp::BasicSequence(name, sequence, alpha)
      {

      }
    SequenceProxy::~SequenceProxy()
      {
        // TODO Auto-generated destructor stub
      }

  }
