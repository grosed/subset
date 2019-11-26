

#include <iostream>
#include <tuple>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <list>

namespace subset
{

  std::tuple<std::vector<double>,unsigned int,double,std::list<unsigned int> > subset(const std::vector<std::vector<double> >& Y,
										      const unsigned int& s,
										      const unsigned int& e,
										      const std::vector<double>& betas,
										      const std::vector<double>& alphas,
										      const std::vector<double>& thresholds)
  {
    // assumes error checking done by calling environment
    auto p = Y.size();
    auto n = e - s + 1;
    auto beta = betas[e-s];
    auto alpha = alphas[e-s];
    std::vector<std::vector<double> > S(p);    
    // create the cumulative sums
    std::transform(Y.begin(),Y.end(),S.begin(),[&n,&s,&e](auto& y)
		   {
		     std::vector<double> cusum(n);
		     auto left = y.begin();
		     auto right = y.begin();
		     std::advance(left,s);
		     std::advance(right,e+1);
		     std::partial_sum(left,right, cusum.begin(), std::plus<double>());
		     return(cusum);
		   });
    
    std::vector<std::vector<double> > C(p);
    std::transform(S.begin(),S.end(),C.begin(),[&n](auto& s)
		   {
		     unsigned int l = 1;
		     unsigned int m = n-1;
		     std::vector<double> res(n-1);
		     std::transform(s.begin(),s.end()-1,res.begin(),[&l,&m,&n,&s](auto& v)
				    {
				      return (v*v/(double)(l++)) + ((s[n-1]-v)*(s[n-1]-v)/(double)(m--)) - (s[n-1]*s[n-1]/(double)(n));   
				    });
		     return res;
		   });
    
    std::vector<std::vector<double> > D(p);
    std::transform(C.begin(),C.end(),D.begin(),[&alpha](auto& c)
		   {
		     std::vector<double> d(c.size());
		     std::transform(c.begin(),c.end(),d.begin(),[&alpha](auto& x)
				    {
				      return x > alpha ? x - alpha : 0.0;
				    });
		     return d;
		   });
    
    
    // auto thresh = (double)(p) + std::sqrt(2.0*(double)(p)*beta);
    auto thresh = thresholds[e-s];
    std::vector<double> Ctot(n-1,0.0);
    for(auto& c : C)
      {
	std::transform(c.begin(),c.end(),Ctot.begin(),Ctot.begin(),[](auto& a, auto& b)
		       {
			 return a + b;
		       });
      }
    std::transform(Ctot.begin(),Ctot.end(),Ctot.begin(),[&thresh](auto& x)
		   {
		     return x - thresh;
		   });
    
    
    std::vector<double> Dtot(n-1,0.0);
    for(auto& d : D)
      {
	std::transform(d.begin(),d.end(),Dtot.begin(),Dtot.begin(),[](auto& a, auto& b)
		       {
			 return a + b;
		       });
      }
    
    auto calcValues = [&beta](auto& x)
      {
	std::vector<double> values(x.size());
	std::transform(x.begin(),x.end(),values.begin(),[&beta](auto& val)
      {
	return val - beta;
      });
	return values;
      };
    
    std::vector<double> Dvalues = calcValues(Dtot);
    std::vector<double> Cvalues = calcValues(Ctot);
    
    std::vector<double> savings(Dvalues.size());
    std::transform(Cvalues.begin(),Cvalues.end(),Dvalues.begin(),savings.begin(),[](auto& a, auto& b)
		   {
		     return a > b ? a : b;
		   });
    
    auto itmax = std::max_element(savings.begin(),savings.end());
    auto argmax = std::distance(savings.begin(),itmax);
    if(*itmax < 0)
      {
	return std::make_tuple(Dtot,0,0.0,std::list<unsigned int>());
      }
    else if(Dtot[argmax] < *itmax)
      {
	std::list<unsigned int> affected;
	for(unsigned int i = 0; i < p; i++)
	  {
	    affected.push_back(i);
	  }
	return std::make_tuple(savings,argmax,*itmax,affected);
      }
    else
      {
	std::list<unsigned int> affected;
	for(unsigned int i = 0; i < D.size(); i++)
	  {
	    if(D[i][argmax] > 0.0)
	      {
		affected.push_back(i);
	      }
	  }
	return std::make_tuple(savings,argmax,*itmax,affected);
      }  
    
  }

} // namespce subset


