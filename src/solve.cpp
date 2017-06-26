#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <RcppEigen.h>
#include <vector>
#include <amgcl/runtime.hpp>

typedef Eigen::Map< Eigen::VectorXd > mV ;
typedef Eigen::SparseMatrix<double> sM;

// [[Rcpp::export]]
Rcpp::NumericVector amgsolve(const std::vector<int>    ptr,const std::vector<int>  col,const std::vector<double> val,const std::vector<double> rhs,std::vector<double> guess,double tol=1e-5,int maxiter=100,const std::string coarsening="smoothed_aggregation",const std::string relax="spai0",const std::string solver="bicgstab")
{
  Rcpp::NumericVector out;
  int    iters;
  double error;

  boost::property_tree::ptree prm;
  prm.put("precond.coarsening.type",coarsening);
  prm.put("precond.relax.type",relax);
  prm.put("solver.type",solver);
  prm.put("solver.tol",tol);
  prm.put("solver.maxiter",maxiter);
  std::size_t n = rhs.size();
  typedef amgcl::backend::builtin<double> Backend;
  amgcl::make_solver<
      amgcl::runtime::amg<Backend>,
      amgcl::runtime::iterative_solver<Backend>
      > solve(boost::tie(n, ptr, col, val), prm);
    

  boost::tie(iters, error) = solve(rhs, guess);
  out = Rcpp::wrap(guess);
  out.attr("iter") = iters;
  out.attr("error") = error;

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector amgsolve_ns(const std::vector<int>    ptr,const std::vector<int>  col,const std::vector<double> val,const std::vector<double> rhs,std::vector<double> guess,Rcpp::NumericMatrix nS, double tol=1e-5,int maxiter=100,const std::string coarsening="smoothed_aggregation",const std::string relax="spai0",const std::string solver="bicgstab",int max_levels=4,int coarse_enough=1)
{
  Rcpp::NumericVector out;
  int    iters;
  double error;

  boost::property_tree::ptree prm;
  prm.put("precond.coarsening.type",coarsening);
  prm.put("precond.relax.type",relax);
  prm.put("solver.type",solver);
  prm.put("solver.tol",tol);
  prm.put("solver.maxiter",maxiter);
  prm.put("precond.coarse_enough",coarse_enough);
  prm.put("precond.max_levels",max_levels);
  //Add null space
  std::vector<double> nsvec = Rcpp::as<std::vector<double> >(nS);
  prm.put("precond.coarsening.nullspace.cols", nS.ncol());
  prm.put("precond.coarsening.nullspace.rows", nS.nrow());
  prm.put("precond.coarsening.nullspace.B", &nsvec[0]);

  std::size_t n = rhs.size();
  typedef amgcl::backend::builtin<double> Backend;
  amgcl::make_solver<
      amgcl::runtime::amg<Backend>,
      amgcl::runtime::iterative_solver<Backend>
      > solve(boost::tie(n, ptr, col, val), prm);
    

  boost::tie(iters, error) = solve(rhs, guess);
  out = Rcpp::wrap(guess);
  out.attr("iter") = iters;
  out.attr("error") = error;

  return out;
}
