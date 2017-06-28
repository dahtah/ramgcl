[![Travis-CI Build Status](https://travis-ci.org/dahtah/ramgcl.svg?branch=master)](https://travis-ci.org/dahtah/ramgcl)
# ramgcl

[amgcl](amgcl.readthedocs.io) by Denis Demidov is a C++ library for solving sparse linear systems using iterative methods. The specific technique used is an Algebraic Multigrid preconditioner, along with various classical solvers like Conjugate Gradients. Although it may not work equally well on all matrices, it is blazingly fast on certain problems like finite difference/finite element systems. If you have enough memory, you can solve problems involving millions of variables.

*ramgcl* is an interface to amgcl. If the name sounds to you like a random string of consonants, try pronoucing it "ram-GCL" in your head (that's what I do). 

# Basic usage

To solve an equation, first set up a preconditioner: 

```{r}
library(ramgcl)
#Random sparse, 2x2 matrix
A <- sparseMatrix(1:2,1:2,x=rnorm(2))
#Set up a preconditioner
P <- precond(A)
P
```


Then call psolve:

```{r}
b <- rnorm(2)
psolve(P,b)
##Checking solution:
solve(A,b)
```

See package vignette for more. 

# Direct vs. iterative methods

The direct methods for sparse matrices in R's Matrix package are already pretty good, so when can you expect a speed-up? Here's an example (beware, it's memory heavy):

```{r}
library(igraph)
#One step of the 3D heat equation on a 50x50x50 grid
L <- graph.lattice(rep(50,3)) %>% laplacian_matrix
M <- L+.1*Diagonal(nrow(L))
#Initial condition
u <- rnorm(nrow(L))
C <- Cholesky(M,super=TRUE)
v <- solve(C,u) %>% as.vector
v2 <- precond(M) %>% psolve(u)
all.equal(v,v2)
```

On my machine amgcl takes ~80ms on that problem (compared to a few seconds for direct methods). Bear in mind that direct methods yield much more information, e.g. you get the determinant for free. 


# Current status

Working, with some missing features. ramgcl uses the native back-end, which is parallelised using OpenMP but does not support GPUs. If somebody feels like implementing GPU support via e.g., the ViennaCL back-end, they should definitely get in touch. 
No support yet for specifying a near null-space. 

# Author

Simon Barthelme, Gipsa-lab, CNRS.
