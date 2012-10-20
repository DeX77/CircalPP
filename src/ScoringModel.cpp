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

#include "ScoringModel.h"
#include "AlignmentSymbol.h"

#include <Bpp/Text/TextTools.h>
#include <algorithm>
#include <sstream>

namespace Circal
  {
    ScoringModel::ScoringModel() :
        model()
      {
      }
    ScoringModel::ScoringModel(const std::string & path) :
        model()
      {
        read(path);
      }

    ScoringModel::ScoringModel(std::istream & input) :
        model()
      {
        read(input);
      }

    ScoringModel::~ScoringModel()
      {
      }
    void ScoringModel::read(std::istream & input)
      {
        if (!input)
          {
            throw bpp::IOException(
                "ScoringModel::ScoringModel: fail to open file");
          }

        std::string temp;
        std::string type;

        double gapStart = 0;
        double gapExtend = 0;
        double match = 0;
        double missmatch = 0;
        double negativeMatch = 0;

        double maxValue = 1;

        uint Size = 1;

        // Main loop : for all file lines
        while (!input.eof())
          {
            getline(input, temp, '\n'); // Copy current line in temporary string

            if (temp[0] != '#')
              {
                std::istringstream istring(temp);
                istring >> type >> gapStart >> gapExtend >> match
                    >> negativeMatch >> missmatch;
                istring >> temp; //Cut out the ":" char

                if (match > maxValue)
                  maxValue = match;

                if (negativeMatch > maxValue)
                  maxValue = negativeMatch;

                while (istring)
                  {
                    AlignmentSymbol sc;
                    istring >> temp;
                    if (!temp.empty())
                      {
                        sc.symbol = bpp::TextTools::toUpper(temp);
                        if (sc.symbol.length() > Size)
                          Size = sc.symbol.length();
                        sc.type = type;
                        //Changed gapOpen & gapExtend to negative Values to correspond with Maximize
                        sc.gapOpen = -gapStart;
                        sc.gapExtend = -gapExtend;
                        sc.match = match;
                        sc.negativeMatch = negativeMatch;
                        sc.missmatch = -missmatch;

                        model.insert(
                            std::pair<const std::string, AlignmentSymbol>(
                                sc.symbol, sc));
                        //Also add negative Symbol Version
                        sc.symbol = "-" + sc.symbol;
                        model.insert(
                            std::pair<const std::string, AlignmentSymbol>(
                                sc.symbol, sc));
                      }
                  }

              }
          }
        NormalizeScores(maxValue);

      }

    void ScoringModel::read(const std::string &path)
      {
        std::ifstream input(path.c_str(), std::ios::in);
        read(input);
        input.close();
      }

    double ScoringModel::ScoreOf(const std::string &A, const std::string &B)
      {

        double outA = 0;
        double outB = 0;

        //Check if two symbols are completly identical
        if (A.compare(B) == 0)
          {
            return model[A].match;
          }
        else
        //Check for negative values
        if ((A.compare("-" + B) == 0) || (B.compare("-" + A) == 0))
          {
            return model[A].negativeMatch;
          }
        else
          {
            outA = model[A].missmatch;
            outB = model[B].missmatch;
            //Check wich one is worse
            if (BestOfTwo(outA, outB) != outA)
              return outA;
            else
              return outB;
          }
      }

    double ScoringModel::ScoreOfGapOpen(const std::string A)
      {
        return model[A].gapOpen;
      }
    double ScoringModel::ScoreOfGapExtend(const std::string A)
      {
        return model[A].gapExtend;
      }
    int ScoringModel::GetModelSize(void) const
      {
        return model.size();
      }

    double ScoringModel::BestOfTwo(const double A, const double B) const
      {

        return std::max(A, B);
      }

    double ScoringModel::BestOfThree(const double A, const double B,
        const double C) const
      {

        return BestOfTwo(A, BestOfTwo(B, C));

      }

    ModelValues::const_iterator ScoringModel::constitStart(void) const
      {
        return model.begin();
      }

    ModelValues::const_iterator ScoringModel::constitEnd(void) const
      {
        return model.end();
      }

    ModelValues::iterator ScoringModel::itStart(void)
      {
        return model.begin();
      }
    ModelValues::iterator ScoringModel::itEnd(void)
      {
        return model.end();
      }

    void ScoringModel::NormalizeScores(const double &maxValue)
      {
        for (ModelValues::iterator foo = itStart(); foo != itEnd(); foo++)
          {
            foo->second.match /= maxValue;
            foo->second.negativeMatch /= maxValue;

          }
      }
  }
