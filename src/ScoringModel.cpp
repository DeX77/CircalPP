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

namespace Circal
  {
    ScoringModel::ScoringModel()
      {
      }
    ScoringModel::ScoringModel(const std::string & path)
      {
        read(path);
      }

    ScoringModel::ScoringModel(std::istream & input)
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
            throw bpp::IOException ("ScoringModel::ScoringModel: fail to open file");
          }

        string temp;
        string type;

        double gapStart;
        double gapExtend;
        uint Size = 1;

        // Main loop : for all file lines
        while (!input.eof())
          {
            getline(input, temp, '\n'); // Copy current line in temporary string

            if (temp[0] != '#')
              {
                istringstream istring(temp);
                istring >> type >> gapStart >> gapExtend;
                istring >> temp; //Cut out the ":" char

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
                        model[sc.symbol] = sc;
                        //Also add negative Symbol Version
                        sc.symbol = "-" + sc.symbol;
                        model[sc.symbol] = sc;
                      }
                  }

              }
          }
      }

    void ScoringModel::read(const std::string &path)
      {
        ifstream input(path.c_str(), ios::in);
        read(input);
        input.close();
      }

    double ScoringModel::ScoreOf(const std::string &A, const std::string &B) const
      {
        //Check if two symbols are completly identical
        if (A.compare(B) == 0)
          return 20;
        else
          return -10;
      }

    double ScoringModel::ScoreOfGapOpen(const std::string A) const
      {
        ModelValues::const_iterator foo = model.find(A);
        return foo->second.gapOpen;
      }
    double ScoringModel::ScoreOfGapExtend(const std::string A) const
      {
        ModelValues::const_iterator foo = model.find(A);
        return foo->second.gapExtend;
      }
    int ScoringModel::GetModelSize(void) const
      {
        return model.size();
      }

    const double &ScoringModel::BestOfTwo(const double &A, const double &B) const
      {
        //If we minimize
        //return std::min(A, B);
        //If we maximize
        return std::max(A, B);
      }

    const double &ScoringModel::BestOfThree(double &A, double &B, double &C) const
      {

        return BestOfTwo(A, BestOfTwo(B, C) );

      }
    ModelValues::const_iterator ScoringModel::itBegin(void) const
      {
        return model.begin();
      }
    ModelValues::const_iterator ScoringModel::itEnd(void) const
      {
        return model.end();
      }
  }
