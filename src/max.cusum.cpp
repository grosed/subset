
#include <cstdlib>
#include <iostream>
#include <vector>
#include <functional>
#include <tuple>
#include <list>
#include <set>
#include <algorithm>
#include <fstream>
#include <limits>
#include <numeric>
#include <cmath>


namespace subset
{

  std::tuple<unsigned int,double> max_cusum(const std::vector<double>& X,const unsigned int& s,const unsigned int& e)
  {
    double Xrhs = 0.0;
    for(unsigned int b = s; b < e; b++)
      {
	Xrhs += X[b];
      }
    double Xlhs = 0.0;
    auto n = e - s + 1;
    auto max = -std::numeric_limits<double>::infinity();
    unsigned int argmax = s;
    for(unsigned int b = s; b < e; b++)
      {
	Xlhs += X[b];
	Xrhs -= X[b];
	auto Xbse = std::fabs(std::sqrt((double)(e-b) / (double)(n*(b-s+1)))*Xlhs - std::sqrt((double)(b-s+1) / (double)(n*(e-b)))*Xrhs);
	if(Xbse > max)
	  {
	    max = Xbse;
	    argmax = b;
	  }
      }
    return std::make_pair(argmax,max);
  }

 
} // namespace subset
