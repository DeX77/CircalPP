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

#include "MultipleCircularAlignmentFactory.h"
#include "MultiplePseudoCircularAlignmentFactory.h"
#include "ScoringModel.h"
#include "ExtendedFasta.h"
#include "Output.h"
#include "RandomSequence.h"
#include "VertebrateMitochondrialGenomeAlphabet.h"
#include <Seq/alphabets>
#include <iostream>
#include <fstream>
#include <sstream>

int main(int args, char* argv[])
  {

    //Check for correct call
    if (args < 5)
      {
        std::cerr
            << "Usage: Circal++ DeltaValue Input.fasta  XX.score Output.fasta"
            << std::endl;;
        exit(-1);
      }

    std::stringstream parser(argv[1]);
    int delta;
    parser >> delta;

    std::string seqFilename = argv[2];
    std::string scoreFilename = argv[3];
    std::string outFilename = argv[4];

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

    Circal::Output* out = new Circal::Output();

    std::string gotoh = "GOTOH";
    double Alignments = 0;
    double AlignmentsGood = 0;
    double AlignmentsRev = 0;

    Circal::CircularAlignmentFactory* cicAln =
        new Circal::CircularAlignmentFactory();
    Circal::PseudoCircularAlignmentFactory* pseudoCircAln =
        new Circal::PseudoCircularAlignmentFactory();

    //#ifdef _OPENMP            
    //#pragma omp parallel for
    //#endif      
    for (uint u=0; u<sequences->getNumberOfSequences(); u++)
      {

        for (uint k=u+1; k<sequences->getNumberOfSequences(); k++)
          {

            Circal::Alignment* regCircular = cicAln->GotohAlignment(
                sequences->getSequence(u), sequences->getSequence(k), scoreM);
            Circal::Alignment* pseudoCircular= pseudoCircAln->GotohAlignment(
                sequences->getSequence(u), sequences->getSequence(k), scoreM,
                delta);
            if (regCircular->get_Score() != pseudoCircular->get_Score())
              {

              }
            else
              AlignmentsGood++;
            Alignments++;
            std::cout << "Optimum: " << std::endl;
            std::cout << out->AlignmentPrettyPrint(regCircular) << std::endl;
            std::cout << "Heuristik: " << std::endl;
            std::cout << out->AlignmentPrettyPrint(pseudoCircular) << std::endl;
            std::cout << "Insgesamt: " << Alignments << " davon schlecht: "
                << (Alignments-AlignmentsGood) << " Rate: " << (AlignmentsGood
                /Alignments) << std::endl;
          }
      }

    /*
     Circal::MultiplePseudoCircularAlignment* circal;
     //Create Multiple Circular Alignment
     std::cerr << "Pseudo circuläres Alignment" << std::endl;
     circal = new Circal::MultiplePseudoCircularAlignment(sequences,scoreM);

     
     Circal::Output* foo = new Circal::Output();
     seqWriter.write(outFilename, *dynamic_cast<bpp::SequenceContainer*>(circal));

     std::cout << foo->TCoffeeLibFormat(circal);
     */

    //Clean up
    delete alpha;
    delete sequences;
    delete out;
    delete scoreM;
    delete cicAln;
    delete pseudoCircAln;

    return 0;

    /*
     //Check for correct call
     if (args < 4)
     {
     std::cerr<< "Usage: Circal++ XX.score Runs maxSize" << std::endl;;
     exit(-1);
     }
     std::string scoreFilename = argv[1];
     std::stringstream parser1(argv[2]);
     std::stringstream parser2(argv[3]);
     uint runs;
     parser1 >> runs;

     uint maxSize;
     parser2 >> maxSize;

     //Check for Scoring File Existence
     if (!bpp::FileTools::fileExists(scoreFilename))
     {
     std::cerr << "Can't open Scoring File "<< scoreFilename << std::endl;
     exit(-1);
     }
     //Init Scoring Model
     Circal::ScoringModel* scoreM = new Circal::ScoringModel(scoreFilename);
     bpp::DNA* alpha = new bpp::DNA();

     Circal::Output* out = new Circal::Output();

     double Alignments = 0;
     double AlignmentsGood = 0;
     double GlobalAlignments = 0;
     double GlobalAlignmentsGood = 0;
     std::string gotoh = "GOTOH";

     #ifdef _OPENMP            
     #pragma omp parallel for
     #endif 
     for (uint j=1; j < maxSize; j++)
     {
     for (uint k=j+1; k < maxSize; k++)
     {
     Alignments = 0;
     AlignmentsGood = 0;
     for (uint i=0; i < runs; i++)
     {
     Circal::RandomSequence* test1 = new Circal::RandomSequence(j,alpha);
     Circal::RandomSequence* test2 = new Circal::RandomSequence(k,alpha);
     Circal::CircularAlignment* cicAln =
     new Circal::CircularAlignment(gotoh, test1, test2, scoreM);
     Circal::PseudoCircularAlignment* pseudoCircAln =
     new Circal::PseudoCircularAlignment(gotoh, test1, test2, scoreM);
     Circal::PseudoCircularAlignment* pseudoCircAlnRev =
     new Circal::PseudoCircularAlignment(gotoh, test2, test1, scoreM);
     if (cicAln->get_Score()> pseudoCircAln->get_Score())
     {
     std::cout << "Optimum: " << std::endl;
     std::cout << out->AlignmentPrettyPrint(cicAln) << std::endl;
     std::cout << "Heuristik: " << std::endl;
     std::cout << out->AlignmentPrettyPrint(pseudoCircAln)
     << std::endl;
     }
     else if (cicAln->get_Score()> pseudoCircAlnRev->get_Score())
     {
     std::cout << "Optimum: " << std::endl;
     std::cout << out->AlignmentPrettyPrint(cicAln) << std::endl;
     std::cout << "Heuristik: " << std::endl;
     std::cout << out->AlignmentPrettyPrint(pseudoCircAlnRev)
     << std::endl;
     }
     else
     AlignmentsGood++;
     Alignments++;

     }
     std::cout << "Laenge 1= " << j << " Laenge 2= " << k
     << " schlecht: " << (Alignments-AlignmentsGood) << " Rate: "
     << (AlignmentsGood /Alignments) << std::endl;
     GlobalAlignmentsGood += AlignmentsGood;
     GlobalAlignments += Alignments;
     }
     std::cout << "Isgesamt: " << GlobalAlignments << " davon schlecht: "
     << (GlobalAlignments-GlobalAlignmentsGood) << " Rate: "
     << (GlobalAlignmentsGood/GlobalAlignments) << std::endl;

     }


     //Clean up
     delete alpha;
     delete scoreM;
     */
  }

