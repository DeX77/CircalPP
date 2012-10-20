#include "RandomSequence.h"

namespace Circal
  {

    RandomSequence::RandomSequence(const uint &sequence_size, const bpp::Alphabet* alpha) :
      SequenceProxy("Random Sequence", std::vector<int>(), alpha)
      {
        int k=-10;
        for (uint i=0; i< sequence_size; i++)
          {
            //This makes sure we add a valid symbol
            k = rand() % 4;
            this->addElement(k);
          }
      }

    RandomSequence::~RandomSequence()
      {
      }

  }
