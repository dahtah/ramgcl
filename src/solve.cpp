#include <ramgcl.h>

using namespace Rcpp;


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


// Return the preconditioned solver only
// [[Rcpp::export]]
Rcpp::XPtr<Solver > amgsolver(int n, const std::vector<int>    ptr,const std::vector<int>  col,const std::vector<double> val,double tol=1e-5,int maxiter=100,const std::string coarsening="smoothed_aggregation",const std::string relax="spai0",const std::string solver="bicgstab")
{
  int    iters;
  double error;

  boost::property_tree::ptree prm;
  prm.put("precond.coarsening.type",coarsening);
  prm.put("precond.relax.type",relax);
  prm.put("solver.type",solver);
  prm.put("solver.tol",tol);
  prm.put("solver.maxiter",maxiter);
  std::size_t m = (std::size_t) n;
  //  typedef amgcl::backend::builtin<double> Backend;
  Solver *solve_ptr;
  solve_ptr = new Solver(boost::tie(m, ptr, col, val), prm);
    
  solve_ptr->precond();
  //  boost::tie(iters, error) = solve(rhs, guess);
  // out = Rcpp::wrap(guess);
  // out.attr("iter") = iters;
  // out.attr("error") = error;
  Rcpp::XPtr< Solver > p(solve_ptr,true);
  return p;
}



// Solve using a system matrix that's not the same as the one used in the preconditioner
// [[Rcpp::export]]
Rcpp::NumericVector solve_newmat(Rcpp::XPtr< Solver> solve,const std::vector<int>    ptr,const std::vector<int>  col,const std::vector<double> val,const std::vector<double> rhs,std::vector<double> guess)
{
  Rcpp::NumericVector out;
  int    iters;
  double error;
  std::size_t n = rhs.size();
  boost::tie(iters, error) = (*solve)(boost::tie(n, ptr, col, val),rhs, guess);
  out = Rcpp::wrap(guess);
  out.attr("iter") = iters;
  out.attr("error") = error;

  return out;
}



// solver for systems with a null space (incomplete)
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
  //  typedef amgcl::backend::builtin<double> Backend;
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
Rcpp::NumericVector solve_rhs(Rcpp::XPtr< Solver> solve,const std::vector<double> rhs,std::vector<double> guess)
{
  Rcpp::NumericVector out;
  int    iters;
  double error;

  boost::tie(iters, error) = (*solve)(rhs, guess);
  out = Rcpp::wrap(guess);
  out.attr("iter") = iters;
  out.attr("error") = error;
  return out;
}




// Rcpp::NumericMatrix solve_mrhs(Rcpp::XPtr< Solver> solve,Rcpp::NumericMatrix rhs,Rcpp::NumericMatrix guess)
// {
//   int nrhs = rhs.ncol();
//   Rcpp::NumericMatrix out;
//   int    iters;
//   double error;

//   for (int i = 0; i < nrhs; i++)
//     {
//       std::vector< double> r = rhs(_,i);
//       std::vector< double> g =guess(_,i);
//       boost::tie(iters, error) = (*solve)(r,g);
//     }
//   out = Rcpp::wrap(guess);
//   return out;
// }
