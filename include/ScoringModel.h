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

#include "AlignmentSymbol.h"
//Try tu include unordered set
#include <tr1/unordered_map>

//Didn't work? Ok we use map instead
#ifndef _TR1_UNORDERED_MP
#include <map>
#endif

#include <string>

namespace Circal
  {    

#ifdef _TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<std::string, AlignmentSymbol> ModelValues;
#else
    typedef std::map<const std::string, AlignmentSymbol> ModelValues;

#endif
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

      double ScoreOf(const std::string &A, const std::string &B);
      double ScoreOfGapOpen(const std::string A) ;
      double ScoreOfGapExtend(const std::string A) ;
      int GetModelSize(void) const;

      double BestOfTwo(const double A, const double B) const;
      double BestOfThree(const double A, const double B, const double C) const;

      void NormalizeScores(const double &maxValue);

      ModelValues::const_iterator constitStart(void) const;
      ModelValues::const_iterator constitEnd(void) const;

      ModelValues::iterator itStart(void);
      ModelValues::iterator itEnd(void);

      };
  }

#endif /*SCORINGMODEL_H_*/
