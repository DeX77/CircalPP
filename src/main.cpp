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
#include "CorrectedFasta.h"
#include "WhitespaceFasta.h"
#include "Output.h"
#include "RandomSequence.h"
#include "VertebrateMitochondrialGenomeAlphabet.h"
#include <Seq/alphabets>
#include <NumCalc/RandomTools.h>
#include <iostream>
#include <fstream>
#include <sstream>

std::string usage()
  {
    std::stringstream out;

    out<< std::endl;
    out<<"usage: CircalPP [options]" << std::endl;
    out<< std::endl;
    out<<"options:" << std::endl;
    out<<" -v  verbose mode" << std::endl;
    out<<" -d  <integer>  delta value" << std::endl;
    out<<" -D  input is dna" << std::endl;
    out<<" -R  input is rna" << std::endl;
    out<<" -G  alphabet is build from scoring sheme" << std::endl;
    out<<" -S  <filename> file including the scoring sheme" << std::endl;
    out<<" -I  <filename> read data from fasta file instead of stdin"
        << std::endl;
    out<<" -F  <filename> write all pairwise alignments as fasta to <filename>"
        << std::endl;
    out<<" -O <filename> write output to <filename> instead of stdout"
        << std::endl;
    return out.str();
  }
void doAllignment(const bpp::Alphabet* alpha, const std::string &seqFilename,
    const Circal::ScoringModel &scoreM, const std::string &outFilename,
    const std::string &resultFilename, bool outF, bool resultF, int &delta,
    bool verbose)
  {

    if (verbose)
      std::cout << "AlphaTyp: " << alpha->getAlphabetType() << std::endl;

    //Container for read sequences
    bpp::VectorSequenceContainer sequences(alpha);

    //Initialize FASTA Reader
    if (alpha->getAlphabetType() == "Vertebrate Mitochondrial Genome")
      {
        if (verbose)
          std::cout << "Using WhitespaceFasta Reader" << std::endl;
        Circal::WhitespaceFasta seqReader;
        //Check for Sequence File Existence
        if (!bpp::FileTools::fileExists(seqFilename))
          {
            std::cerr << "Can't open Sequence File "<< seqFilename
                << " using stdin instead" << std::endl;
            seqReader.read(std::cin, sequences);
          }
        else
          //Try parsing the file for FASTA
          try
            {
              seqReader.read(seqFilename,sequences);
            }
          catch(bpp::Exception &e)
            {
              std::cerr << "Invalid File Format!" << std::endl;
              std::cerr << e.what() << std::endl;
              exit(-1);
            }
      }
    else
      {
        if (verbose)
          std::cout << "Using CorrectedFasta Reader" << std::endl;

        Circal::WhitespaceFasta seqReader;
        //Check for Sequence File Existence
        if (!bpp::FileTools::fileExists(seqFilename))
          {
            std::cerr << "Can't open Sequence File "<< seqFilename
                << " using stdin instead" << std::endl;
            seqReader.read(std::cin, sequences);
          }
        else
          //Try parsing the file for FASTA
          try
            {
              seqReader.read(seqFilename,sequences);
            }
          catch(bpp::Exception &e)
            {
              std::cerr << "Invalid File Format!" << std::endl;
              std::cerr << e.what() << std::endl;
              exit(-1);
            }
      }

    Circal::CorrectedFasta seqWriter;

    //Create Multiple Circular Alignment
    Circal::MultiplePseudoCircularAlignmentFactory circal;
    Circal::Alignment multi = circal.GotohalignMultiple(&sequences, &scoreM,
        delta, verbose);

    if (outF)
      seqWriter.write(outFilename,
          *dynamic_cast<bpp::SequenceContainer*>(&multi));

    Circal::Output prettyPrint;

    if (resultF)
      {
        ofstream foo(resultFilename.c_str());
        foo << prettyPrint.TCoffeeLibFormat(&multi, &sequences);
      }
    else
      std::cout << prettyPrint.TCoffeeLibFormat(&multi, &sequences);
  }

int main(int args, char* argv[])
  {
    bool verbose = false;
    bool outF = false;
    bool resultF = false;
    bool alphaSet = false;
    bool scoreSet = false;

    std::string alphaTyp;
    std::string seqFilename;
    std::string scoreFilename;
    std::string outFilename;
    std::string resultFilename;
    int delta = 1;

    for (int i=1; i<args; i++)
      {
        //        std::cout << argv[i] << std::endl;

        if (argv[i][0]=='-')
          switch (argv[i][1])
            {
          case 'v':
            verbose = true;
            break;
          case 'd':
            {
              std::stringstream parser(argv[i+1]);
              parser >> delta;
              if (verbose)
                std::cout << "Delta: " << delta << std::endl;
              break;
            }
          case 'D':
            {
              alphaTyp = "DNA";
              alphaSet = true;
              if (verbose)
                std::cout << "Using DNA alphabet" << std::endl;
              break;
            }
          case 'R':
            {

              alphaTyp = "RNA";

              alphaSet = true;
              if (verbose)
                std::cout << "Using RNA alphabet" << std::endl;
              break;
            }
          case 'G':
            {
              //Use Self defined Alphabet
              alphaTyp = "SELF";

              if (verbose)
                std::cout << "Using self defined alphabet" << std::endl;
              alphaSet = true;
              break;
            }
          case 'S':
            {
              scoreFilename = argv[i+1];
              if (verbose)
                std::cout << "ScoreFile: " << scoreFilename << std::endl;
              //Check for Scoring File Existence
              if (!bpp::FileTools::fileExists(scoreFilename))
                {
                  std::cerr << "Can't open Scoring File "<< scoreFilename
                      << std::endl;
                  exit(-1);
                }
              scoreSet = true;
              break;
            }
          case 'F':
            {
              outFilename = argv[i+1];
              if (verbose)
                std::cout << "Output Fasta-File: " << outFilename << std::endl;
              outF = true;
              break;
            }
          case 'O':
            {
              resultFilename = argv[i+1];
              if (verbose)
                std::cout << "Output File: " << resultFilename << std::endl;
              resultF = true;
              break;
            }
          case 'I':
            {
              seqFilename = argv[i+1];
              if (verbose)
                std::cout << "InputFile: " << seqFilename << std::endl;
              break;
            }

          default:
            std: cerr << usage();
            }
      }

    //Check for Alphabet and Scoring Sheme set
    if (!alphaSet | !scoreSet)
      {
        std::cerr << usage();
        exit(-1);
      }

    if ((alphaTyp == "SELF") && (scoreFilename == ""))
      {
        std::cerr << usage();
        exit(-1);
      }

    //Init Alphabet
    if (alphaTyp == "DNA")
      {
        bpp::DNA* alpha = new bpp::DNA();

        //Init Scoring Model
        Circal::ScoringModel scoreM(scoreFilename);

        doAllignment(alpha, seqFilename, scoreM, outFilename, resultFilename,
            outF, resultF, delta, verbose);
        delete alpha;
      }
    else if (alphaTyp == "RNA")
      {
        bpp::RNA* alpha = new bpp::RNA();

        //Init Scoring Model
        Circal::ScoringModel scoreM(scoreFilename);

        doAllignment(alpha, seqFilename, scoreM, outFilename, resultFilename,
            outF, resultF, delta, verbose);
        delete alpha;
      }
    else
      {
        //Init Scoring Model
        Circal::ScoringModel scoreM(scoreFilename);

        Circal::VertebrateMitochondrialGenomeAlphabet* alpha =
            new Circal::VertebrateMitochondrialGenomeAlphabet(&scoreM);

        doAllignment(alpha, seqFilename, scoreM, outFilename, resultFilename,
            outF, resultF, delta, verbose);
        delete alpha;
      }

    //WTF!=?? Why do I have to delete this!?
    delete bpp::RandomTools::DEFAULT_GENERATOR;
    return 0;

  }

