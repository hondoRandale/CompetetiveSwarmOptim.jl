# MIT License

#Copyright (c) 2021 hondoRandale <jules.rasetaharison@tutanota.com> and contributors

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

#######################
## 2D-test functions ##
#######################

## Beale function
limits_Beale_2D(d) = ( repeat( [-4.5f0], d ), repeat( [4.5f0], d ) )
Beale_2D(x;kwargs) = ( 1.5f0 - x[1] + x[1]*x[2] )^2 + ( 2.25f0 - x[1] + x[1]*x[2]^2 )^2 + ( 2.625f0 - x[1] + x[1]*x[2]^3 )^2
min_x_Beale_2D(d)  = [ 3.0f0, 0.5f0 ];

## Himmelblau's function
limits_Himmelblau_2D(d) = ( repeat( [-5.0f0], d ), repeat( [5.0f0], d ) )
Himmelblau_2D(x;kwargs) = ( x[1]^2 + x[2] - 11.0f0 )^2 + ( x[1] + x[2]^2 - 7.0f0  )^2
min_x_Himmelblau_2D(d)  = [ 3.0f0 -2.805118f0 -3.779310f0 3.584428f0;
                            2.0f0 3.131312f0  -3.283186f0 -1.848126f0 ];

## Lévi function N.13
limits_Levi_no13_2D(d) = ( repeat( [-10.0f0], d ), repeat( [10.0f0], d ) )
Levi_no13_2D(x;kwargs) = sin( 3.0f0 * π * x[1] )^2 + ( x[1] - 1.0f0 )^2 * ( 1.0f0 + sin( 1.0f0+3*π*x[2] )^2 ) +
                         ( x[2] - 1.0f0 )^2 * ( 1.0f0 + sin( 2*π*x[2] )^2 )
min_x_Levi_no13_2D(d)  = [ 1.0f0, 1.0f0 ]

## Rosenbrock_function
limits_Rosenbrock_2D(d) = ( repeat( [-10.0f0], d ), repeat( [10.0f0], d ) )
Rosenbrock_2D(x;kwargs) = ( 1.0f0 - x[1] )^2 + 100.0f0 * (x[2] - x[1]^2)^2
min_x_Rosenbrock_2D(d)  = [ 1.0f0, 1.0f0 ]

## Booth_function
limits_Booth_2D(d) = ( repeat( [-10.0f0], d ), repeat( [10.0f0], d ) )
Booth_2D(x;kwargs) = ( x[1] + 2.0f0*x[2] - 7.0f0 )^2 +
                     ( 2.0f0*x[1] + x[2] - 5.0f0 )^2
min_x_Booth_2D(d)  = [ 1.0f0, 3.0f0 ]

## Ackley_function
limits_Ackley_2D(d) = ( repeat( [-5.0f0], d ), repeat( [5.0f0], d ) )
Ackley_2D(x;kwargs) = -20.0f0 * exp( -0.2f0 * sqrt( 0.5f0 * ( x[1]^2 + x[2]^2 ) )  ) -
                      exp( 0.5f0 * ( cos( 2.0f0 * π * x[1] ) + cos( 2.0f0 * π * x[2] ) ) ) +
                      exp( 1.0f0 ) +
                      20.0f0
min_x_Ackley_2D(d)  = [ 0.0f0, 0.0f0 ]

## Matyas_function
limits_Matyas_2D(d) = ( repeat( [-10.0f0], d ), repeat( [10.0f0], d ) )
Matyas_2D(x;kwargs) = 0.26f0 * ( x[1]^2 + x[2]^2  ) - 0.48f0 * x[1] * x[2]
min_x_Matyas_2D(d)  = [ 0.0f0, 0.0f0 ]

## Three_hump_camel_function
limits_Three_hump_camel_2D(d) = ( repeat( [-5.0f0], d ), repeat( [5.0f0], d ) )
Three_hump_camel_2D(x;kwargs) = 2.0f0 * x[1]^2 - 1.05f0 * x[1]^4 +
                                x[1]^6 / 6 + x[1] * x[2] + x[2]^2
min_x_Three_hump_camel_2D(d)  = [ 0.0f0, 0.0f0 ]

## Easom_function
limits_Easom_2D(d) = ( repeat( [-100.0f0], d ), repeat( [100.0f0], d ) )
Easom_2D(x;kwargs) = -cos( x[1] ) * cos( x[2] ) * exp( -( (x[1]-π)^2 + (x[2]-π)^2 ) )
min_x_Easom_2D(d)  = [ Float32(π), Float32(π) ]

## Sphere_function
limits_Sphere_2D(d) = ( repeat( [-100.0f0], d ), repeat( [100.0f0], d ) )
Sphere_2D(x;kwargs) = x[1]^2 + x[2]^2
min_x_Sphere_2D(d)  = [ 0.0f0, 0.0f0 ]

## Schaffer_function_no2
limits_Schaffer_function_no2_2D(d) = ( repeat( [-100.0f0], d ), repeat( [100.0f0], d ) )
Schaffer_function_no2_2D(x;kwargs) = 0.5f0 + ( sin( x[1]^2 - x[2]^2 )^2 - 0.5f0 ) /
                                             ( 1.0f0 + 0.001f0*( x[1]^2 + x[2]^2 ) )^2
min_x_Schaffer_function_no2_2D(d)  = [ 0.0f0, 0.0f0 ]

## Bukin_function_no6
limits_Bukin_function_no6_2D(d) = ( [ -15.0f0, -3.0f0 ], [ -5.0f0, 3.0f0 ] )
Bukin_function_no6_2D(x;kwargs) = 100.0f0 * sqrt( abs( x[2] - 0.01f0 * x[1]^2 ) ) + 0.01f0*abs( x[1] + 10.0f0 )
min_x_Bukin_function_no6_2D(d)  = [ -10.0f0, 1.0f0 ]

lim_functions_2D = [ limits_Beale_2D,
                     limits_Himmelblau_2D, limits_Levi_no13_2D,
                     limits_Ackley_2D, limits_Matyas_2D, limits_Booth_2D,
                     limits_Three_hump_camel_2D, limits_Easom_2D,
                     limits_Sphere_2D, limits_Rosenbrock_2D,
                     limits_Schaffer_function_no2_2D,
                     limits_Bukin_function_no6_2D ]

obj_functions_2D = [ Beale_2D,
                     Himmelblau_2D, Levi_no13_2D,
                     Ackley_2D, Matyas_2D, Booth_2D,
                     Three_hump_camel_2D, Easom_2D,
                     Sphere_2D, Rosenbrock_2D,
                     Schaffer_function_no2_2D,
                     Bukin_function_no6_2D ];

min_functions_2D = [ min_x_Beale_2D,
                     min_x_Himmelblau_2D, min_x_Levi_no13_2D,
                     min_x_Ackley_2D, min_x_Matyas_2D, min_x_Booth_2D,
                     min_x_Three_hump_camel_2D, min_x_Easom_2D,
                     min_x_Sphere_2D, min_x_Rosenbrock_2D,
                     min_x_Schaffer_function_no2_2D,
                     min_x_Bukin_function_no6_2D ]

#######################
## ND-test functions ##
#######################

## Rastrigin function
limits_Rastrigin_ND(d) = ( repeat( [-5.12f0], d ), repeat( [5.12f0], d ) )
Rastrigin_ND(x;kwargs) = 10.0f0 * Float32( length( x ) ) + sum( ( x.^2 ) .- 10.0f0 * cos.( 2*π*x )  )
min_x_Rastrigin_ND(d)  = zeros( Float32, d )

## Sphere function
limits_Sphere_ND(d) = ( repeat( [-5.12f0], d ), repeat( [5.12f0], d ) )
Sphere_ND(x;kwargs) = sum( x.^2 )
min_x_Sphere_ND(d)  = zeros( Float32, d )

## Styblinski–Tang function
limits_Styblinski_Tang_ND(d) = ( repeat( [-5.0f0], d ), repeat( [5.0f0], d ) )
Styblinski_Tang_ND(x;kwargs) = 0.5f0 * sum( (x).^4 .- 16.0f0*(x).^2 .+ 5.0f0.*x )
min_x_Styblinski_Tang_ND(d)  = repeat( [-2.903534f0], d )

## Axis Parallel Hyper-Ellipsoid function
limits_AxisParallelHyperEll_ND(d) = ( repeat( [-1000000.0f0], d ), repeat( [1000000.0f0], d ) )
AxisParallelHyperEll_ND(x;kwargs) = sum( Float32.( 1:1:length( x ) ) .* ( x.^2 ) )
min_x_AxisParallelHyperEll_ND(d)  = zeros( Float32, d )

## Zakharov Function
limits_Zakharov_ND(d) = ( repeat( [-5.0f0], d ), repeat( [10.0f0], d ) )
Zakharov_ND(x;kwargs) = sum( x.^2 ) + sum( 0.5f0 * Float32.( 1:1:length( x ) ) .* x )^2 + sum( 0.5f0 * Float32.( 1:1:length( x ) ) .* x )^4
min_x_Zakharov_ND(d)  = zeros( Float32, d )

## Michalewicz function
limits_Michalewicz_ND(d) = ( repeat( [0.0f0], d ), repeat( [Float32(π)], d ) )
Michalewicz_ND(x;kwargs) = -sum( sin.( x ) .* sin.( ( Float32.( 1:1:length( x ) ) .* ( x.^2 ) ) / π ) )
min_x_Michalewicz_ND(d)  = repeat( [NaN32], d )

## Exponential function
limits_Exponential_ND(d) = ( repeat( [-1.0f0], d ), repeat( [1.0f0], d ) )
Exponential_ND(x;kwargs) = -exp( -0.5f0 * sum( x.^2 ) )
min_x_Exponential_ND(d)  = zeros( Float32, d )

## Griewank function
limits_Griewank_ND(d) = ( repeat( [-600.0f0], d ), repeat( [600.0f0], d ) )
Griewank_ND(x;kwargs) = 1.0f0 + sum( (x).^2 / 4000.0f0 ) - prod( cos.( x ./ sqrt.( Float32.( 1:1:length( x ) ) ) ) )
min_x_Griewank_ND(d)  = zeros( Float32, d )

## Qing function
limits_Qing_ND(d) = ( repeat( [-500.0f0], d ), repeat( [500.0f0], d ) )
Qing_ND(x;kwargs) = sum( ( (x).^2 .- Float32.( 1:1:length( x ) ) ).^2 )
min_x_Qing_ND(d)  = zeros( Float32, d )

## Salomon function
limits_Salomon_ND(d) = ( repeat( [-100.0f0], d ), repeat( [100.0f0], d ) )
Salomon_ND(x;kwargs) = 1.0f0 - cos( 2*π*sqrt( sum( (x).^2 ) ) ) + 0.1f0 * sqrt( sum( (x).^2 ) )
min_x_Salomon_ND(d)  = zeros( Float32, d )

## Schwefel function
limits_Schwefel_ND(d) = ( repeat( [-500.0f0], d ), repeat( [500.0f0], d ) )
Schwefel_ND(x;kwargs) = 418.9829f0 * length( x ) - sum( x .* sin.( sqrt.( abs.( x ) ) ) )
min_x_Schwefel_ND(d)  = repeat( [420.9687f0], d )

## Schwefel 2.20 function
limits_Schwefel_2_2_0_ND(d) = ( repeat( [-100.0f0], d ), repeat( [100.0f0], d ) )
Schwefel_2_2_0_ND(x;kwargs) = sum( abs.( x ) )
min_x_Schwefel_2_2_0_ND(d)  = zeros( Float32, d )

## Schwefel 2.21 function
limits_Schwefel_2_2_1_ND(d) = ( repeat( [-100.0f0], d ), repeat( [100.0f0], d ) )
Schwefel_2_2_1_ND(x;kwargs) = maximum( abs.( x ) )
min_x_Schwefel_2_2_1_ND(d)  = zeros( Float32, d )

## Schwefel 2.22 function
limits_Schwefel_2_2_2_ND(d) = ( repeat( [-100.0f0], d ), repeat( [100.0f0], d ) )
Schwefel_2_2_2_ND(x;kwargs) = sum( abs.( x ) ) + prod( abs.( x ) )
min_x_Schwefel_2_2_2_ND(d)  = zeros( Float32, d )

## Schwefel 2.23 function
limits_Schwefel_2_2_3_ND(d) = ( repeat( [-10.0f0], d ), repeat( [10.0f0], d ) )
Schwefel_2_2_3_ND(x;kwargs) = sum( (x).^10 )
min_x_Schwefel_2_2_3_ND(d)  = zeros( Float32, d )

lim_functions_ND = [ limits_Rastrigin_ND, limits_Sphere_ND,
                     limits_Styblinski_Tang_ND, limits_AxisParallelHyperEll_ND,
                     limits_Zakharov_ND, limits_Michalewicz_ND,
                     limits_Exponential_ND, limits_Griewank_ND,
                     limits_Qing_ND, limits_Salomon_ND,
                     limits_Schwefel_ND, limits_Schwefel_2_2_0_ND,
                     limits_Schwefel_2_2_1_ND, limits_Schwefel_2_2_2_ND,
                     limits_Schwefel_2_2_3_ND ];

obj_functions_ND = [ Rastrigin_ND, Sphere_ND,
                     Styblinski_Tang_ND, AxisParallelHyperEll_ND,
                     Zakharov_ND, Michalewicz_ND,
                     Exponential_ND, Griewank_ND,
                     Qing_ND, Salomon_ND,
                     Schwefel_ND, Schwefel_2_2_0_ND,
                     Schwefel_2_2_1_ND, Schwefel_2_2_2_ND,
                     Schwefel_2_2_3_ND ];

min_functions_ND = [ min_x_Rastrigin_ND, min_x_Sphere_ND,
                     min_x_Styblinski_Tang_ND, min_x_AxisParallelHyperEll_ND,
                     min_x_Zakharov_ND, min_x_Michalewicz_ND,
                     min_x_Exponential_ND, min_x_Griewank_ND,
                     min_x_Qing_ND,  min_x_Salomon_ND,
                     min_x_Schwefel_ND, min_x_Schwefel_2_2_0_ND,
                     min_x_Schwefel_2_2_1_ND, min_x_Schwefel_2_2_2_ND,
                     min_x_Schwefel_2_2_3_ND ];
