#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ramgcl_amgsolve_ns(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ramgcl_amgsolver(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ramgcl_solve_newmat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ramgcl_solve_rhs(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ramgcl_amgsolve_ns",  (DL_FUNC) &ramgcl_amgsolve_ns,  13},
    {"ramgcl_amgsolver",    (DL_FUNC) &ramgcl_amgsolver,     9},
    {"ramgcl_solve_newmat", (DL_FUNC) &ramgcl_solve_newmat,  6},
    {"ramgcl_solve_rhs",    (DL_FUNC) &ramgcl_solve_rhs,     3},
    {NULL, NULL, 0}
};

void R_init_ramgcl(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
