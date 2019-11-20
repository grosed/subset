
#ifndef ___SUBSET_H___
#define  ___SUBSET_H___

#include <tuple>
#include <vector>
#include <list>

namespace subset
{


  std::tuple<std::vector<double>,unsigned int,double,std::list<unsigned int> > subset(const std::vector<std::vector<double> >&,
										      const unsigned int&,
										      const unsigned int&,
										      const std::vector<double>&,
										      const std::vector<double>&);
  

} // namespace subset



#endif
