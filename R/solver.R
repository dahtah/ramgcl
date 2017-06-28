asolve <- function(A,b,guess,null,tol=1e-5,maxiter=100,coarsen="smoothed_aggregation",relax="spai0",max_level=4,solver="gmres")
{
    assertChoice(coarsen,c("ruge_stuben","aggregation","smoothed_aggregation","smoothed_aggr_emin"))
    assertChoice(relax,c("gauss_seidel","multicolor_gauss_seidel","ilu0","parallel_ilu0","iluk",
                         "ilut","damped_jacobi","spai0","spai1","chebyshev"))
    if (is(b,"dgeMatrix")) b <- as.vector(b)
    assertVector(b,len=nrow(A))
    if (missing(guess))
    {
        guess <- rep(0,length(b))
    }
    else
    {
        
        if ((length(guess) == 1) && (is.character(guess)))
        {
            #Provide other guesses?
            assertChoice(guess,"b")
            if (guess=="b")
            {
                guess <- b
            }
        }
        else
        {
            assertVector(guess,len=nrow(A))
        }
    }

    #R gives us CSC format, amgcl wants CSR
    At <- Matrix::t(A)

    if (!missing(null))
    {
        if (is.vector(null)) null <- matrix(null,length(null),1)
        assertMatrix(null,min.cols=1,nrows=nrow(A))
        if (!isSymmetric(A))
            {
                stop("Problems with null spaces are only supported on symmetric systems")
            }
        #Project b onto the subspace 
        bs <- svd(null,ncol(null),0)$u
        bproj <- as.vector(b-bs%*%(t(bs)%*%b))
        out <- amgsolve_ns(At@p,At@i,At@x,rhs=bproj,guess=0*bproj,tol=tol,maxiter=maxiter,coarsen=coarsen,relax=relax,nS=null,max_level=max_level,solver=solver)
        attr(out,"b.proj") <- bproj
        out
    }
    else
    {
        amgsolve(At@p,At@i,At@x,rhs=b,guess=guess,tol=tol,maxiter=maxiter,coarsen=coarsen,relax=relax,solver=solver)
    }
}


##' Solve a linear system iteratively using  algebraic multigrid preconditioning
##'
##' Interface to Denis Demidov's amgcl library. There are two main functions, "precond" and "psolve". precond builds the preconditioner (an approximation to the system matrix) and sets up the solver. psolve runs the preconditioner on a specific right-hand side.
##'
##' You can change the system matrix when running psolve (using the optional A argument ). The original preconditioner will be retained, which means you save the cost of setting up a new preconditioner, but the preconditioner is no longer optimally suited to the system. Use this only for small modifications to the system matrix.
##' 
##' @param A sparse, square matrix (from Matrix package)
##' @param b right-hand side (numeric, real)
##' @param guess (optional) initial value for the iteration
##' @param tol tolerance
##' @param maxiter maximum number of iterations
##' @param coarsen coarsening method to use (one of the following: ruge_stuben,aggregation,smoothed_aggregation,smoothed_aggr_emin)
##' @param relax relaxation to use (one of the following: gauss_seidel,ilu0,iluk,ilut,damped_jacobi,spai0,spai1,chebyshev)
##' @param solver which iterative solver to use (cg, idrs, lgmres, gmres, bicgstab,bicgstabl). cg only works on SPD matrices, the rest allow non-symmetric matrices. 
##' @return a vector x such that Ax=b (or as far as the iteration would allow), with attributes "iter" (number of iterations), and "error" (root mean squared residual).
##' @author Simon Barthelme
##' @export
##' @examples
##' A <- sparseMatrix(1:2,1:2,x=rnorm(2))
##' b <- rnorm(2)
##' P <- precond(A)
##' psolve(P,b)
##' #Compare to: 
##' solve(A,b)
##'
##' 
#' @export
precond <- function(A,tol=1e-5,maxiter=100,coarsen="smoothed_aggregation",relax="spai0",max_level=4,solver="gmres")
{
    assertChoice(coarsen,c("ruge_stuben","aggregation","smoothed_aggregation","smoothed_aggr_emin"))
    assertChoice(relax,c("gauss_seidel","ilu0","iluk",
                         "ilut","damped_jacobi","spai0","spai1","chebyshev"))
    assertChoice(solver,c("cg","idrs","lgmres","gmres","bicgstab","bicgstabl"))
    #R gives us CSC format, amgcl wants CSR
    At <- Matrix::t(A)
    xptr <- amgsolver(nrow(A),At@p,At@i,At@x,tol=tol,maxiter=maxiter,coarsen=coarsen,relax=relax,solver=solver)
    P <- list(xptr=xptr,coarsen=coarsen,relax=relax,max_level=max_level,solver=solver,n=nrow(A))
    class(P) <- "amgcl"
    P
}


#' @export
print.amgcl <- function(x,...)
{
    msg <- sprintf("amgcl preconditioner for matrix of size %i \n",x$n)
    cat(msg)
    invisible(x)
}


#' @describeIn precond solve a preconditioned system
#' @export
psolve <- function(P,b,guess,A)
{
    assertClass(P,"amgcl")
    assertAtomic(b,len=P$n)
    if (missing(guess))
    {
        guess <- 0*b
    }
    assertAtomic(guess,len=P$n)
    if (missing(A))
    {
        solve_rhs(P$xptr,b,guess)
    }
    else
    {
        assertClass(A,"dgCMatrix")
        if (nrow(A) != P$n || ncol(A) != P$n)
        {
            stop("A must be a square matrix of the same dimension as the preconditioner")
        }
        At <- Matrix::t(A)
        solve_newmat(P$xptr,At@p,At@i,At@x,b,guess)
    }
}


