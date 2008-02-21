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

#include "ExtendedFasta.h"

namespace Circal
  {
    ExtendedFasta::ExtendedFasta()
      {
      }

    ExtendedFasta::~ExtendedFasta()
      {
      }

    //    void ExtendedFasta::appendFromStream(istream & input,
    //        bpp::VectorSequenceContainer & vsc) const throw(bpp::Exception)
    //      {
    //        if (!input)
    //          {
    //            throw bpp::IOException ("Fasta::read: fail to open file");
    //          }
    //
    //        string temp, name, base = ""; // Initialization
    //        vector<int> sequence;
    //
    //        // Main loop : for all file lines
    //        while (!input.eof())
    //          {
    //            getline(input, temp, '\n'); // Copy current line in temporary string
    //
    //
    //            // If first character is >
    //            if (temp[0] == '>')
    //              {
    //                // If a name and a sequence were found
    //                if ((name != "") && (!sequence.empty()))
    //                  {
    //                    // New sequence creation, and addition in existing VectorSequenceContainer
    //                    bpp::Sequence* seq = new bpp::Sequence(name, sequence, vsc.getAlphabet());
    //                    vsc.addSequence(*seq);
    //                    delete seq;
    //                  }
    //                // Sequence name isolation
    //                name = temp;
    //                name.erase(name.begin()); // Character > deletion
    //              }
    //            else if (!temp.empty())// Sequence isolation
    //              {
    //                sequence.clear();
    //                stringstream seqstream(temp);
    //
    //                while (seqstream)
    //                  {
    //                    seqstream >> base;
    //                    base.resize(
    //                        bpp::AlphabetTools::getAlphabetCodingSize(vsc.getAlphabet()),
    //                        ' ');
    //                    sequence.push_back(vsc.getAlphabet()->charToInt(base));
    //                  }
    //                //Seems like the last added is double.. remove it
    //                sequence.pop_back();
    //              }
    //          }
    //
    //        // Addition of the last sequence in file
    //        if ((name != "") && (!sequence.empty()))
    //          {
    //            bpp::Sequence * seq = new bpp::Sequence(name, sequence, vsc.getAlphabet());
    //            vsc.addSequence(*seq);
    //            delete seq;
    //          }
    //      }
    void ExtendedFasta::appendFromStream(istream & input,
        bpp::VectorSequenceContainer & vsc) const throw (bpp::Exception)
      {
        if (!input)
          {
            throw bpp::IOException ("Fasta::read: fail to open file");
          }

        string temp, name, sequence = ""; // Initialization

        // Main loop : for all file lines
        while (!input.eof())
          {
            getline(input, temp, '\n'); // Copy current line in temporary string

            // If first character is >
            if (temp[0] == '>')
              {
                // If a name and a sequence were founded
                if ((name != "") && (sequence != ""))
                  {
                    // New sequence creation, and addition in existing VectorSequenceContainer
                    bpp::Sequence * seq = new bpp::Sequence(name, sequence, vsc.getAlphabet());
                    vsc.addSequence(*seq);
                    delete seq;
                    //name = ""; no need for that, no?
                    sequence = "";
                  }
                // Sequence name isolation
                name = temp;
                name.erase(name.begin()); // Character > deletion
              }
            else
              sequence += temp; // Sequence isolation
          }

        // Addition of the last sequence in file
        if ((name != "") && (sequence != ""))
          {
            bpp::Sequence * seq = new bpp::Sequence(name, sequence, vsc.getAlphabet());
            vsc.addSequence(*seq);
            delete seq;
          }
      }
  }
