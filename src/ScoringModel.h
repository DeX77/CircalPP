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

#ifndef SCORINGMODEL_H_
#define SCORINGMODEL_H_

#include <map>
#include <string>

namespace Circal
  {
    class AlignmentSymbol;
    typedef std::map<const std::string, AlignmentSymbol> ModelValues;

    class ScoringModel
      {

      ModelValues model;
  public:

      ScoringModel();
      explicit ScoringModel(const std::string & path);
      explicit ScoringModel(std::istream & input);
      virtual ~ScoringModel();
      void read(const std::string & path);
      void read(std::istream & input);

      double ScoreOf(const std::string &A, const std::string &B) const;
      double ScoreOfGapOpen(const std::string A) const;
      double ScoreOfGapExtend(const std::string A) const;
      int GetModelSize(void) const;

      const double &BestOfTwo(const double &A, const double &B) const;
      const double &BestOfThree(double &A, double &B, double &C) const;
      ModelValues::const_iterator itBegin(void) const;
      ModelValues::const_iterator itEnd(void) const;

      };
  }

#endif /*SCORINGMODEL_H_*/
