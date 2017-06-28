// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/ramgcl.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// amgsolve
Rcpp::NumericVector amgsolve(const std::vector<int> ptr, const std::vector<int> col, const std::vector<double> val, const std::vector<double> rhs, std::vector<double> guess, double tol, int maxiter, const std::string coarsening, const std::string relax, const std::string solver);
RcppExport SEXP ramgcl_amgsolve(SEXP ptrSEXP, SEXP colSEXP, SEXP valSEXP, SEXP rhsSEXP, SEXP guessSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP coarseningSEXP, SEXP relaxSEXP, SEXP solverSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int> >::type ptr(ptrSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type col(colSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type val(valSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type guess(guessSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const std::string >::type coarsening(coarseningSEXP);
    Rcpp::traits::input_parameter< const std::string >::type relax(relaxSEXP);
    Rcpp::traits::input_parameter< const std::string >::type solver(solverSEXP);
    rcpp_result_gen = Rcpp::wrap(amgsolve(ptr, col, val, rhs, guess, tol, maxiter, coarsening, relax, solver));
    return rcpp_result_gen;
END_RCPP
}
// amgsolver
Rcpp::XPtr<Solver > amgsolver(int n, const std::vector<int> ptr, const std::vector<int> col, const std::vector<double> val, double tol, int maxiter, const std::string coarsening, const std::string relax, const std::string solver);
RcppExport SEXP ramgcl_amgsolver(SEXP nSEXP, SEXP ptrSEXP, SEXP colSEXP, SEXP valSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP coarseningSEXP, SEXP relaxSEXP, SEXP solverSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type ptr(ptrSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type col(colSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type val(valSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const std::string >::type coarsening(coarseningSEXP);
    Rcpp::traits::input_parameter< const std::string >::type relax(relaxSEXP);
    Rcpp::traits::input_parameter< const std::string >::type solver(solverSEXP);
    rcpp_result_gen = Rcpp::wrap(amgsolver(n, ptr, col, val, tol, maxiter, coarsening, relax, solver));
    return rcpp_result_gen;
END_RCPP
}
// solve_newmat
Rcpp::NumericVector solve_newmat(Rcpp::XPtr< Solver> solve, const std::vector<int> ptr, const std::vector<int> col, const std::vector<double> val, const std::vector<double> rhs, std::vector<double> guess);
RcppExport SEXP ramgcl_solve_newmat(SEXP solveSEXP, SEXP ptrSEXP, SEXP colSEXP, SEXP valSEXP, SEXP rhsSEXP, SEXP guessSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< Solver> >::type solve(solveSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type ptr(ptrSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type col(colSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type val(valSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type guess(guessSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_newmat(solve, ptr, col, val, rhs, guess));
    return rcpp_result_gen;
END_RCPP
}
// amgsolve_ns
Rcpp::NumericVector amgsolve_ns(const std::vector<int> ptr, const std::vector<int> col, const std::vector<double> val, const std::vector<double> rhs, std::vector<double> guess, Rcpp::NumericMatrix nS, double tol, int maxiter, const std::string coarsening, const std::string relax, const std::string solver, int max_levels, int coarse_enough);
RcppExport SEXP ramgcl_amgsolve_ns(SEXP ptrSEXP, SEXP colSEXP, SEXP valSEXP, SEXP rhsSEXP, SEXP guessSEXP, SEXP nSSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP coarseningSEXP, SEXP relaxSEXP, SEXP solverSEXP, SEXP max_levelsSEXP, SEXP coarse_enoughSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int> >::type ptr(ptrSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type col(colSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type val(valSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type guess(guessSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type nS(nSSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const std::string >::type coarsening(coarseningSEXP);
    Rcpp::traits::input_parameter< const std::string >::type relax(relaxSEXP);
    Rcpp::traits::input_parameter< const std::string >::type solver(solverSEXP);
    Rcpp::traits::input_parameter< int >::type max_levels(max_levelsSEXP);
    Rcpp::traits::input_parameter< int >::type coarse_enough(coarse_enoughSEXP);
    rcpp_result_gen = Rcpp::wrap(amgsolve_ns(ptr, col, val, rhs, guess, nS, tol, maxiter, coarsening, relax, solver, max_levels, coarse_enough));
    return rcpp_result_gen;
END_RCPP
}
// solve_rhs
Rcpp::NumericVector solve_rhs(Rcpp::XPtr< Solver> solve, const std::vector<double> rhs, std::vector<double> guess);
RcppExport SEXP ramgcl_solve_rhs(SEXP solveSEXP, SEXP rhsSEXP, SEXP guessSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< Solver> >::type solve(solveSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type guess(guessSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_rhs(solve, rhs, guess));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"ramgcl_amgsolve", (DL_FUNC) &ramgcl_amgsolve, 10},
    {"ramgcl_amgsolver", (DL_FUNC) &ramgcl_amgsolver, 9},
    {"ramgcl_solve_newmat", (DL_FUNC) &ramgcl_solve_newmat, 6},
    {"ramgcl_amgsolve_ns", (DL_FUNC) &ramgcl_amgsolve_ns, 13},
    {"ramgcl_solve_rhs", (DL_FUNC) &ramgcl_solve_rhs, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ramgcl(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
