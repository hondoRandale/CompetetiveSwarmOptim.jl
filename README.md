# CompetetiveSwarm.jl
   [![TagBot](https://github.com/hondoRandale/CompetetiveSwarm.jl/actions/workflows/TagBot.yml/badge.svg)](https://github.com/hondoRandale/CompetetiveSwarm.jl/actions/workflows/TagBot.yml)

   [![documentation](https://github.com/hondoRandale/CompetetiveSwarm.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/hondoRandale/CompetetiveSwarm.jl/actions/workflows/documentation.yml)

## Introduction
   This solver is initialized using simpleFWA algorithm.
   in the following iteraion each cluster is updated using an 
   attraction/repulsion mechanism, taking into account distance and
   difference in the objective function value.
___

## calling convention
   each objective function passed to CompetetiveSwarm has to comply with the following
   simple parameter convention f( x; kwargs ) where f is the objective
   function to be minimized. This convention ensures CompetetiveSwarm can be used with
   time-series-problems, classification-problems, regression-problems.
   Univariate as well as multivariate target sets are admissible.

___
## example

```julia
   using CompetetiveSwarm
   using Test
   Easom(x;kwargs) = -cos( x[1] ) * cos( x[2] ) *
                     exp( -( (x[1]-π)^2 + (x[2]-π)^2 ) )
   lower    = [ -10.0f0, -10.0f0 ];
   upper    = [ 10.0f0, 10.0f0 ];
   cFWA( objFunction ) = competetiveSwarm( 16, 16, ();
                                           λ_0           = 7.95f0,
                                           ϵ_A           = 0.5f-2,
                                           C_a           = 1.2f0,
                                           C_r           = 0.8f0,
                                           lower         = lower,
                                           upper         = upper,
                                           objFunction   = objFunction,
                                           maxiter       = 40,
                                           global_forces = false )                             
   solutionCFWA = cFWA( Easom );
   @test isapprox( solutionCFWA.x_b[1], π; atol=0.01 )
   @test isapprox( solutionCFWA.x_b[2], π; atol=0.01 )                             
```

___
## function reference

```@docs
CompetetiveSwarm.competetiveSwarm
```

```@docs
CompetetiveSwarm.competetiveSwarm
```

FWA struct

| Parameter         | Description                           | Type                      |
| :---              | :---                                  | :---                      |
| X                 | each column is the origin of a fw     | Matrix{Float32}           |
| fitnessX          | fitness of each fw                    | Vector{Float32}           |
| A                 | explosion amplitudes                  | Vector{Float32}           |
| S                 | contains all sparks foreach fw        | Vector{ Matrix{Float32} } |
| fitness_sparks    | fitness of each spark                 | Vector{ Vector{Float32} } |
| x_b               | best found solution                   | Vector{Float32}           |
| y_min             | function value at best found solution | Float32                   |
| iter              | number of iterations executed         | Int                       |

