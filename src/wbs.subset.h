
#ifndef ___WBS_SUBSET_H___
#define ___WBS_SUBSET_H___

#include <list>
#include <vector>
#include <tuple>
#include <functional>
#include "interval.subset.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

namespace subset
{
  
  struct wbs_subset_state
  {
    typedef interval_subset interval_type;
    typedef std::vector<interval_type> intervals_type;
    intervals_type Is;  
    // std::list<unsigned int> bs;
    std::list<interval_type> bs;
    unsigned int seed;
    unsigned int M;
    unsigned int T;
    std::function<std::tuple<std::vector<double>,unsigned int,double,std::list<unsigned int> >(const unsigned int&, const unsigned int&) > X_tilde;
    double zeta;
  };


  template <class state_type>
    state_type create_wbs_subset_intervals(state_type S)
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
	  // calculate statistics and break points
	  auto res = S.X_tilde(I.s,I.e);
	  I.b = std::get<1>(res) + I.s;
	  I.Xb = std::get<2>(res);
	  I.affected = std::get<3>(res);
	  S.Is[i] = I;
	}
      return S;
    }

  


  
}



# endif
