##' Solve linear system iteratively using  algebraic multigrid preconditioning
##'
##' Interface to Denis Demidov's amgcl library
##' @param A sparse, square matrix
##' @param b right-hand side
##' @param guess (optional) initial value for the iteration
##' @param null (optional) null space of A, for singular systems
##' @param tol tolerance
##' @param maxiter maximum number of iterations
##' @param coarsen coarsening method to use (one of the following: ruge_stuben,aggregation,smoothed_aggregation,smoothed_aggr_emin)
##' @param relax relaxation to use (one of the following: gauss_seidel,multicolor_gauss_seidel,ilu0,parallel_ilu0,iluk,ilut,damped_jacobi,spai0,spai1,chebyshev)
##' @return a vector x such that Ax=b (or as far as the iteration would allow)
##' @author Simon Barthelme
##' @export
##' @examples
##' A <- sparseMatrix(1:2,1:2,x=rnorm(2))
##' b <- rnorm(2)
##' asolve(A,b)
##' #Compare to: 
##' solve(A,b)
##' #Use on singular matrices:
##' u <- sparseVector(rnorm(3),1:3,3)
##' A <- u%*%t(u)
##' #A has a null space:
##' null <- svd(A,3)$v[,2:3]
##' A%*%null
##' b <- A%*%rnorm(3)
##' asolve(A,b,null=null,max_level=1)
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
            checkChoice(guess,"b")
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
