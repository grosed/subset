
#ifndef ___MAX_CUSUM_H___
#define ___MAX_CUSUM_H___

#include <tuple>
#include <vector>

namespace subset
{
  std::tuple<unsigned int,double> max_cusum(const std::vector<double>&,const unsigned int&,const unsigned int&);
} // namespace subset


#endif
