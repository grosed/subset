

#include <Rcpp.h>

#include <list>
#include "wbs.h"
#include "max.cusum.h"
#include "subset.h"
#include "wbs.subset.h"


//[[Rcpp::export]]
Rcpp::List marshall_subset(const std::vector<std::vector<double> >& Y,
			   const unsigned int& s,
			   const unsigned int& e,
			   const std::vector<double>& betas,
			   const std::vector<double>& alphas)
{ 

  auto result = subset::subset(Y,s,e,betas,alphas);
  return Rcpp::List::create(Rcpp::Named("savings") = std::get<0>(result),
			    Rcpp::Named("cpt") = std::get<1>(result),
			    Rcpp::Named("cost") = std::get<2>(result),
			    Rcpp::Named("affected") = std::get<3>(result));
  
}



//[[Rcpp::export]]
Rcpp::List marshall_wbs_subset(const std::vector<std::vector<double> >& data,
			       const unsigned int& M,
			       const double& zeta,
			       const std::vector<double>& betas,
			       const std::vector<double>& alphas,
			       const unsigned int& seed)
{ 

  using namespace subset;
  using namespace std::placeholders;
  wbs_subset_state S;
  S.M = M;
  S.T = data[0].size();
  S.seed = seed;
  S.X_tilde = std::bind(subset::subset,data,_1,_2,betas,alphas);
  S.zeta = zeta;
  S = create_wbs_subset_intervals(std::move(S));
  S = wbs(std::move(S));
  // extract the breakpoints from the state
  std::list<unsigned int> breakpoints;
  std::list<std::list<unsigned int> > variates;
  for(auto& I : S.bs)
    {
      breakpoints.push_back(I.b);
      variates.push_back(I.affected);
    }

  return Rcpp::List::create(Rcpp::Named("breakpoints") = breakpoints,
			    Rcpp::Named("affected") = variates);
  
}







//[[Rcpp::export]]
std::list<unsigned int> marshall_wbs(const std::vector<double>& data,const unsigned int& M,const double& zeta,const unsigned int& seed)
{
  using namespace subset;
  using namespace std::placeholders;
  wbs_state S;
  S.M = M;
  S.T = data.size();
  S.seed = seed;
  S.X_tilde = std::bind(max_cusum,data,_1,_2);
  S.zeta = zeta;
  S = create_wbs_intervals(std::move(S));  
  S = wbs(std::move(S));
  // extract the breakpoints from the state
  std::list<unsigned int> breakpoints;
  for(auto& I : S.bs)
    {
      breakpoints.push_back(I.b);
    }
  return breakpoints;
}


