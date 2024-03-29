[![CI](https://github.com/xhub/EMP.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/xhub/EMP.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/xhub/EMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/xhub/EMP.jl)


# EMP

**This package is not for public use as it needs to be updated to the latest version of JuMP/MOI**.

This package enables the modeling of *Extending Mathematical Programming* (EMP) concepts within the JuliaOpt ecosystem.
Broadly speaking, EMP enables the modeling of optimization problems with a structure that do not fit the classical minimization problem.
For instance, the following problems can be modeled with EMP
- Variational Inequalities (VI and AVI) or Complementarity Problems (LCP and NCP)
- (Generalized) Nash Equilibrium Problem (NEP and GNEP) and MOPEC (Multiple Optimization Problems with Equilibrium Constraints)
- Optimization Problems with *Optimal Value Function* (OVF) in the problem data. Examples of OVF include coherent risk measures (CVaR), convex regularizers/loss function (Huber, Hinge, l1, ...)

A problem with EMP can be solved by using model transformation to obtain a form amenable to computations by existing solvers.
Currently, using the ReSHOP library is the only option to solve an optimization problem with EMP data structure.

## Design

This package is designed as a thin layer over an existing modeling framework in JuliaOpt. Right now, only JuMP is supported.
The idea is to use this framework to store variables and equations. The additional information to capture the EMP concept is stored
in an EMP master object.

This package relies on ReSHOP.jl to export the model and the EMP information to the ReSHOP library.
The latter is going to perform the necessary model transformations.

*Note: This package is developed independently of the GAMS corporation*
