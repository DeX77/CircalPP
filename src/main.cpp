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

#ifdef DHAVE_CONFIG_H
#include "config.h"
#endif

#include "MultipleCircularAlignment.h"
#include "MultipleMetaCircularAlignment.h"
#include "ScoringModel.h"
#include "ExtendedFasta.h"
#include "Output.h"
#include "VertebrateMitochondrialGenomeAlphabet.h"
#include <Seq/alphabets>
#include <iostream>
#include <fstream>

int main(int args, char* argv[])
  {
    //Check for correct call
    if (args < 4)
      {
        std::cerr<< "Usage: Circal++ Input.fasta  XX.score Output.fasta"
            << std::endl;;
        exit(-1);
      }
    std::string seqFilename = argv[1];
    std::string scoreFilename = argv[2];
    std::string outFilename = argv[3];

    //Check for Sequence File Existence
    if (!bpp::FileTools::fileExists(seqFilename))
      {
        std::cerr << "Can't open Sequence File "<< seqFilename << std::endl;
        exit(-1);
      }

    //Check for Scoring File Existence
    if (!bpp::FileTools::fileExists(scoreFilename))
      {
        std::cerr << "Can't open Scoring File "<< scoreFilename << std::endl;
        exit(-1);
      }

    //Initialize FASTA Reader
    //Circal::ExtendedFasta seqReader;
    bpp::Fasta seqReader;

    bpp::Fasta seqWriter;

    //Init Scoring Model
    Circal::ScoringModel* scoreM = new Circal::ScoringModel(scoreFilename);

    //Use Self defined Alphabet
    //    Circal::VertebrateMitochondrialGenomeAlphabet * alpha =
    //        new Circal::VertebrateMitochondrialGenomeAlphabet(scoreM);

    bpp::DNA* alpha = new bpp::DNA();

    //Container for read sequences
    bpp::VectorSequenceContainer* sequences;

    //Try parsing the file for FASTA
    try
      {
        sequences = dynamic_cast<bpp::VectorSequenceContainer*>
        (seqReader.read(seqFilename,alpha));
      }
    catch(bpp::Exception &e)
      {
        std::cerr << "Invalid File Format!" << std::endl;
        std::cerr << e.what() << std::endl;
        //Clean up
        delete alpha;
        exit(-1);
      }

    Circal::MultipleMetaCircularAlignment* circal;

    //Create Multiple Circular Alignment
    std::cerr << "Meta circuläres Alignment" << std::endl;
    circal = new Circal::MultipleMetaCircularAlignment(sequences,scoreM);

    Circal::Output* foo = new Circal::Output();
    seqWriter.write(outFilename, *dynamic_cast<bpp::SequenceContainer*>(circal));

    std::cout << foo->TCoffeeLibFormat(circal);

    //Clean up
    delete alpha;
    delete sequences;
    return 0;
  }

