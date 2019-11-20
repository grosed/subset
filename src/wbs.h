
#ifndef ___WBS_H___
#define ___WBS_H___

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
#include "interval.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>


#include <iostream>

namespace subset
{

  struct wbs_state
  {
    typedef Interval interval_type;
    typedef std::vector<interval_type> intervals_type;
    intervals_type Is;  
    // std::list<unsigned int> bs;
    std::list<interval_type> bs;
    unsigned int seed;
    unsigned int M;
    unsigned int T;
    std::function<std::tuple<unsigned int,double>(const unsigned int&, const unsigned int&) > X_tilde;
    double zeta;
  };
  
  template <class state_type>
    state_type create_wbs_intervals(state_type S)
    {
      // initialise random nuber generation
      boost::random::mt19937 gen(S.seed);
      boost::random::uniform_int_distribution<> dist(0,S.T-1);
      S.Is = typename state_type::intervals_type(S.M);
      
      for(int i = 0; i < S.M; i++)
	{
	  typename state_type::interval_type I;
	  I.s = dist(gen);
	  I.e = dist(gen);
	  while(I.e == I.s)
	    {
	      I.e = dist(gen);
	    }
	  if(I.s > I.e)
	    {
	      auto t = I.e;
	      I.e = I.s;
	      I.s = t;
	    }
	  // calculate cusums and break points
	  auto res = S.X_tilde(I.s,I.e);
	  I.b = std::get<0>(res);
	  I.Xb = std::get<1>(res);
	  S.Is[i] = I;
	}
      return S;
    }
  


  // wild binary segmentation
  template <class state_type>
    state_type wbs(state_type S)
    {
      std::function<void(const int& s, const int& e)> wbs_recurse;
      wbs_recurse = [&S,&wbs_recurse](const int& s, const int& e)
	{
	  typename state_type::interval_type argmax; // the interval
	  auto max = -std::numeric_limits<double>::infinity();
	  // the following is the optional part of the WBS algorithm in Fryzlewicz
	  // auto enclosing = S.X_tilde(s,e);
	  //if(std::get<1>(enclosing) > S.zeta && std::get<1>(enclosing) > max)
	  //  {
	  //    max = std::get<1>(enclosing);
	  //    argmax.s = s;
	  //    argmax.e = e;
	  //    argmax.b = std::get<0>(enclosing);
	  //  }
	  for(auto& I : S.Is)
	    {
	      if(I.s >= s && I.s < e && I.e > s && I.e <= e)
		{
		  if(I.Xb > S.zeta && I.Xb > max)
		    {
		      max = I.Xb;
		      argmax = I; // hope we have a sensible copy constructor !!!??
		      // argmax.s = I.s;
		      // argmax.e = I.e;
		      // argmax.b = I.b;
		    }
		}
	    }
	  
	  if(max > -std::numeric_limits<double>::infinity())
	    {
	      // S.bs.push_back(argmax.b);
	      S.bs.push_back(argmax);
	      wbs_recurse(s,argmax.b);
	      wbs_recurse(argmax.b+1,e);
	    }
	};
      int level = 0;
      wbs_recurse(0,S.T-1);
      
      return S;
}
  
  
} // namespace subset




#endif
