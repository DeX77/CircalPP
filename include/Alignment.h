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

#include "SequenceProxy.h"

//Includes from Bio++ 
#include <Bpp/Seq/Container/AbstractSequenceContainer.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Io/IoSequence.h>
#include <Bpp/Io/FileTools.h>

namespace bpp
  {
    class BasicSequence;
  }

namespace Circal
  {

    class Alignment: public bpp::AlignedSequenceContainer
      {
      uint offsetA;
      uint offsetB;
      double Score;

    public:
      Alignment(const bpp::Alphabet* alpha);

      virtual ~Alignment();

      uint get_offsetA();
      void set_offsetA(const uint &orig);

      uint get_offsetB();
      void set_offsetB(const uint &orig);

      double get_Score();
      void set_Score(const double &s);

      void addSequence(const SequenceProxy sequence);
      //void addSequence(const bpp::Sequence & sequence) throw(bpp::Exception);

      void deleteSequence(const std::string& name)
          throw (bpp::SequenceNotFoundException)
        {
          bpp::AlignedSequenceContainer::deleteSequence(name);
        }
      ;
      void deleteSequence(unsigned int sequenceIndex)
          throw (bpp::IndexOutOfBoundsException)
        {
          bpp::AlignedSequenceContainer::deleteSequence(sequenceIndex);
        }
      ;
      std::string toString(unsigned int sequenceIndex) const
          throw (bpp::IndexOutOfBoundsException)
        {
          return bpp::AlignedSequenceContainer::toString(sequenceIndex);
        }
      ;
      std::string toString(const std::string& name) const
          throw (bpp::SequenceNotFoundException)
        {
          return bpp::AlignedSequenceContainer::toString(name);
        }
      ;
      int& operator()(unsigned int sequenceIndex, unsigned int elementIndex)
        {
          return bpp::AlignedSequenceContainer::operator ()(sequenceIndex,
              elementIndex);
        }
      ;
      const int& operator()(unsigned int sequenceIndex,
          unsigned int elementIndex) const
        {
          return bpp::AlignedSequenceContainer::operator ()(sequenceIndex,
              elementIndex);
        }
      ;
      int& operator()(const std::string& sequenceName,
          unsigned int elementIndex)
        {

          return bpp::AlignedSequenceContainer::operator ()(sequenceName,
              elementIndex);
        }
      ;
      const int& operator()(const std::string& sequenceName,
          unsigned int elementIndex) const
        {
          return bpp::AlignedSequenceContainer::operator ()(sequenceName,
              elementIndex);
        }
      ;
      int& valueAt(unsigned int sequenceIndex, unsigned int elementIndex)
          throw (bpp::IndexOutOfBoundsException)
        {
          return bpp::AlignedSequenceContainer::valueAt(sequenceIndex,
              elementIndex);
        }
      ;
      const int& valueAt(unsigned int sequenceIndex,
          unsigned int elementIndex) const
              throw (bpp::IndexOutOfBoundsException)
        {
          return bpp::AlignedSequenceContainer::valueAt(sequenceIndex,
              elementIndex);
        }
      ;
      int& valueAt(const std::string& sequenceName, unsigned int elementIndex)
          throw (bpp::SequenceNotFoundException, bpp::IndexOutOfBoundsException)
        {
          return bpp::AlignedSequenceContainer::valueAt(sequenceName, elementIndex);
        }
      ;
      const int& valueAt(const std::string& sequenceName,
          unsigned int elementIndex) const
              throw (bpp::SequenceNotFoundException,
              bpp::IndexOutOfBoundsException)
        {
          return bpp::AlignedSequenceContainer::valueAt(sequenceName, elementIndex);
        }
      ;
      const bpp::Sequence& getSequence(unsigned int sequenceIndex) const
          throw (bpp::IndexOutOfBoundsException)
        {
          return bpp::AlignedSequenceContainer::getSequence(sequenceIndex);
        }
      ;
      const bpp::Sequence& getSequence(const std::string& name) const
          throw (bpp::SequenceNotFoundException)
        {
          return bpp::AlignedSequenceContainer::getSequence(name);
        }
      ;
      const bpp::Comments& getComments(unsigned int sequenceIndex) const
          throw (bpp::IndexOutOfBoundsException)
        {
          return bpp::AlignedSequenceContainer::getComments(sequenceIndex);
        }
      ;
      const bpp::Comments& getComments(const std::string& name) const
          throw (bpp::SequenceNotFoundException)
        {
          return bpp::AlignedSequenceContainer::getComments(name);
        }
      ;
      const std::vector<int>& getContent(unsigned int sequenceIndex) const
          throw (bpp::IndexOutOfBoundsException)
        {
          return bpp::AlignedSequenceContainer::getContent(sequenceIndex);
        }
      ;
      const std::vector<int>& getContent(const std::string& name) const
          throw (bpp::SequenceNotFoundException)
        {
          return bpp::AlignedSequenceContainer::getContent(name);
        }
      ;
      void setSequencesNames(const std::vector<std::string> & names,
          bool checkNames) throw (bpp::Exception)
        {
          bpp::AlignedSequenceContainer::setSequencesNames(names, checkNames);
        }
      ;
      std::vector<std::string> getSequencesNames() const
        {
          return bpp::AlignedSequenceContainer::getSequencesNames();
        }
      ;
      unsigned int getNumberOfSequences() const
        {
          return bpp::AlignedSequenceContainer::getNumberOfSequences();
        }
      ;
      void setComments(const std::string& name, const bpp::Comments& comments)
          throw (bpp::SequenceNotFoundException)
        {
          bpp::AlignedSequenceContainer::setComments(name, comments);
        }
      ;
      void setComments(unsigned int sequenceIndex,
          const bpp::Comments & comments) throw (bpp::IndexOutOfBoundsException)
        {
          bpp::AlignedSequenceContainer::setComments(sequenceIndex, comments);
        }
      ;
      bpp::Sequence* removeSequence(const std::string& name)
          throw (bpp::SequenceNotFoundException)
        {
          return bpp::AlignedSequenceContainer::removeSequence(name);
        }
      ;
      bpp::Sequence* removeSequence(unsigned int sequenceIndex)
          throw (bpp::IndexOutOfBoundsException)
        {
          return bpp::AlignedSequenceContainer::removeSequence(sequenceIndex);
        }
      ;
      void deleteGeneralComments()
        {
          bpp::AlignedSequenceContainer::deleteGeneralComments();
        }
      ;
      void setGeneralComments(const bpp::Comments& comments)
        {
          bpp::AlignedSequenceContainer::setGeneralComments(comments);
        }
      ;
      const bpp::Comments& getGeneralComments() const
        {
          return bpp::AlignedSequenceContainer::getGeneralComments();
        }
      ;
      const std::string& getName(unsigned int sequenceIndex) const
          throw (bpp::IndexOutOfBoundsException)
        {
          return bpp::AlignedSequenceContainer::getName(sequenceIndex);
        }
      ;
      const bpp::Alphabet* getAlphabet() const
        {
          return bpp::AlignedSequenceContainer::getAlphabet();
        }
      ;
      unsigned int getSequencePosition(const std::string & name) const
          throw (bpp::SequenceNotFoundException)
        {
          return bpp::AlignedSequenceContainer::getSequencePosition(name);
        }
      ;
      bool hasSequence(const std::string& name) const
        {
          return bpp::AlignedSequenceContainer::hasSequence(name);
        }
      ;

      };
  }
#endif /*ALIGNMENT_H_*/
