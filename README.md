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

# Author

Simon Barthelme, Gipsa-lab, CNRS.
