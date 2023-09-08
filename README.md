# Differential Equation 

## Overview
This project aims to discretize a differential equation with variable coefficients, \( \frac{d}{dx}(a(x)u'(x)) = f(x) \), under both Dirichlet and Neumann boundary conditions. The exercise is divided into multiple steps, including discretization of the domain, formulation of a discretization scheme, and accuracy testing using Matlab.

## Table of Contents
- [Discretization of Domain](#discretization-of-domain)
- [Discretization Scheme](#discretization-scheme)
- [Compatibility Conditions](#compatibility-conditions)
- [Accuracy Testing](#accuracy-testing)
- [Analysis](#analysis)
- [Conclusions](#conclusions)

## Discretization of Domain
The domain \([0,1]\) is discretized using a mesh with a step size \( h \), resulting in \( n+2 \) nodes. Additional nodes \( x_{j \pm \frac{1}{2}} \) are introduced for each \( j \).

## Discretization Scheme
A discretization scheme is formulated using the nodes, leading to a linear system of equations. The properties of the coefficient matrix are then investigated.

## Compatibility Conditions
For Neumann boundary conditions, the compatibility condition is derived and its implications for solvability are discussed.

## Accuracy Testing
The accuracy of the approximations is tested using Matlab, and the results are discussed.

## Analysis
1. **Dirichlet Boundary Conditions**: The discretization leads to a linear system \( Au = f \), where \( A \) is a tridiagonal, symmetric, and positive-definite matrix. The Matlab operator "\\" is used for efficient solving.
2. **Neumann Boundary Conditions**: The discretization is similar, but the boundary conditions are incorporated into the linear system. One-sided finite difference approximations are used for \( u'(0) \) and \( u'(1) \).

## Conclusions
The proposed discretization method appears to be quite accurate, although the accuracy varies between Dirichlet and Neumann boundary conditions.
