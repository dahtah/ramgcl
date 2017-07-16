#include <RcppEigen.h>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <vector>
#include <amgcl/runtime.hpp>

typedef Eigen::Map< Eigen::VectorXd > mV ;
typedef Eigen::SparseMatrix<double> sM;
typedef amgcl::backend::builtin<double> Backend;
typedef amgcl::make_solver<amgcl::runtime::amg<Backend>,amgcl::runtime::iterative_solver<Backend> > Solver;

