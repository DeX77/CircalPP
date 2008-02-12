#include "RandomSequence.h"

namespace Circal
  {

    RandomSequence::RandomSequence(const uint &size, const bpp::Alphabet* alpha) :
      bpp::Sequence("Random Sequence", "", alpha)
      {
        int k=-10;
        for (uint i=0; i< size; i++)
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
