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

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

//Includes from Bio++ 
#include <Seq/AbstractSequenceContainer.h>
#include <Seq/ioseq>
#include <Utils/FileTools.h>

namespace bpp
  {
    class Sequence;
  }

namespace Circal
  {

    class ScoringModel;
    class MatrixHelper;
    class Output;
    typedef std::vector< std::vector<double> > ScoreMatrix;

    class Alignment : public bpp::AlignedSequenceContainer
      {
    protected:
        int Score;
        Output* prettyPrint;
        MatrixHelper* matrix;
        uint origSize;
        
    public:
        explicit Alignment(const bpp::Alphabet* alpha);
        virtual ~Alignment();
        uint get_origSize() const;
        int get_Score() const;
        void set_score(const int &s);
        void set_origSize(const uint &orig);
        void addSequence(const bpp::Sequence & sequence, bool checkNames)
            throw(bpp::Exception);
        void addSequence(const bpp::Sequence & sequence)
            throw(bpp::Exception);

        //Needleman Wunsch Alignment
        virtual Alignment* NeedlemanWunschAlignment(const bpp::Sequence* inA,
            const bpp::Sequence* inB, const ScoringModel* scoreM);
        virtual void
            ForwardRecursionNMW(const bpp::Sequence* A,
                const bpp::Sequence* B, const ScoringModel* scoreM,
                ScoreMatrix* D);
        virtual Alignment* BacktrackingNMW(const bpp::Sequence* outA,
            const bpp::Sequence* outB, const ScoringModel* scoreM,
            const ScoreMatrix* D);

        //Gotoh Alignment
        virtual Alignment* GotohAlignment(const bpp::Sequence* inA,
            const bpp::Sequence* inB, const ScoringModel* scoreM);
        virtual void ForwardRecursionGotoh(const bpp::Sequence* A,
            const bpp::Sequence* B, const ScoringModel* scoreM,
            ScoreMatrix* D, ScoreMatrix* P, ScoreMatrix* Q);
        virtual Alignment* BacktrackingGotoh(const bpp::Sequence* outA,
            const bpp::Sequence* outB, const ScoringModel* scoreM,
            const ScoreMatrix* D, const ScoreMatrix* P, const ScoreMatrix* Q);

      };
  }
#endif /*ALIGNMENT_H_*/
