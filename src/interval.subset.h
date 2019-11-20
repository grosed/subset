
#ifndef ___INTERVAL_SUBSET_H___
#define ___INTERVAL_SUBSET_H___

#include <list>

namespace subset
{


  struct interval_subset
  {
    unsigned int s;
    unsigned int e;
    unsigned int b;
    double Xb;
    std::list<unsigned int> affected;
  };


} // namespace subset
 

#endif
