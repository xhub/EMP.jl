[![Build Status](https://travis-ci.com/xhub/EMP.jl.svg?branch=master)](https://travis-ci.com/xhub/EMP.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/x3hjfgt7wiphvf5i?svg=true)](https://ci.appveyor.com/project/xhub/emp-jl)
[![Coverage Status](https://coveralls.io/repos/github/xhub/EMP.jl/badge.svg?branch=master)](https://coveralls.io/github/xhub/EMP.jl?branch=master)
[![codecov](https://codecov.io/gh/xhub/EMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/xhub/EMP.jl)


# EMP

This package enables the modeling of *Extending Mathematical Progmming* (EMP) concepts within the JuliaOpt ecosystem.
Broadly speaking, EMP enables the modeling of optimization problems with a structure that do not fit the classical minimization problem.
For instance, the following problems can be modeled with EMP
- Variational Inequalities (VI and AVI) or Complementarity Problems (LCP and NCP)
- (Generalized) Nash Equilibrium Problem (NEP and GNEP) and MOPEC (Multiple Optimization Problems with Equilibrium Constraints)
- Optimization Problems with *Optimal Value Function* (OVF) in the problem data. Examples of OVF include coherent risk measures (CVaR), convex regularizers (Huber, Hinge, l1, ...)

A problem with EMP can be solved by using model transformation to obtain a form amenable to computations by existing solvers.
Currently, using the JAMSD library is the only option to solve an optimization problem with EMP data structure.

## Design

This package is designed as a thin layer over an existing modeling framework in JuliaOpt. Right now, only JuMP is supported.
The idea is to use this framework to store variables and equations. The additional information to capture the EMP concept is stored
in an EMP master object.

This package relies on JAMSDWriter.jl to export the model and the EMP information to the JAMSD library.
The latter is going to perform the necessary model transformations.
