module CompetetiveSwarmOptim

  include( "BenchmarkOptimProblems.jl" )
  include( "BenchmarkUtility.jl" )
  using Clustering
  using InvertedIndices
  using LinearAlgebra
  using Random
  using Statistics

  export competetiveSwarm, competetiveSwarmDynBoundary, COMPSWARM,
         eval_benchmark_2D, eval_benchmark_2D_Dyn_Bound,
         obj_functions_2D, retrieveLimits2D,
         retrieveLimitsND, retrieveMin2D, retrieveMinND,
         metaObjCSwarm2D_mod0, metaObjCSwarm2D_mod1, metaObjCSwarm2D_mod2,
         Rastrigin_ND, Sphere_ND,
         Styblinski_Tang_ND, AxisParallelHyperEll_ND,
         Zakharov_ND, Michalewicz_ND,
         Exponential_ND, Griewank_ND,
         Qing_ND, Salomon_ND,
         Schwefel_ND, Schwefel_2_2_0_ND,
         Schwefel_2_2_1_ND, Schwefel_2_2_2_ND,
         Schwefel_2_2_3_ND,
         Beale_2D,
         Himmelblau_2D, Levi_no13_2D,
         Ackley_2D, Matyas_2D, Booth_2D,
         Three_hump_camel_2D, Easom_2D,
         Sphere_2D, Rosenbrock_2D,
         Schaffer_function_no2_2D,
         Bukin_function_no6_2D

  @inline function dynExplosionAmplitude( y_min::Float32,
                                          y_cf::Float32,
                                          A_cf_old::Float32,
                                          C_a::Float32,
                                          C_r::Float32 )
    A_cf = A_cf_old::Float32
    if y_min - y_cf < 0.0f0
      A_cf *= C_a;
    else
      A_cf *= C_r;
    end
    return A_cf
  end

   function explosionAmplitudes!(;A::Vector{Float32},
                                 λ_0::Float32,
                                 ϵ_A::Float32,
                                 y_min::Float32,
                                 cacheA::Vector{Float32},
                                 fitness_fireworks::Vector{Float32} )
    n      = length( fitness_fireworks );
    for ii ∈ 1:1:n
      @inbounds fx         = fitness_fireworks[ii];
      @inbounds cacheA[ii] = ( fx - y_min )::Float32;
    end
    denomA = ( sum( cacheA ) + ϵ_A )::Float32;
    for ii ∈ 1:1:n
      @inbounds A[ii] = ( λ_0 * ( ( cacheA[ii] + ϵ_A ) / denomA ) )::Float32;
    end
  end

  function explosionSparks( ;X_i, A_i, s, lower, upper, rng )
    nSparks = size( s, 2 )::Int;
    d       = size( s, 1 )::Int;
    idz     = bitrand( rng, d );
    for j ∈ 1:1:nSparks
      idz .= bitrand( rng, d );
      for ii ∈ 1:1:d
        if idz[ii]
          @inbounds ΔX_k    = ( A_i * ( rand( rng, Float32, 1 )[1] * 2.0f0 - 1.0f0 ) )::Float32;
          @inbounds val     = ( X_i[ii] + ΔX_k )::Float32;
          if val >= lower[ii] && val <= upper[ii]
            @inbounds s[ii,j] = val::Float32;
          else
            @inbounds s[ii,j] = ( rand( rng, Float32, 1 )[1] * ( upper[ii] - lower[ii] ) + lower[ii] )::Float32;
          end
        end
      end
    end
  end

  @views function evaluateCandidates( ;kwargs::Tuple,
                               f::Function,
                               X,
                               fitness )
    @assert size( X, 2 ) == size( fitness, 1 )
    n = size( X, 2 )::Int;
    @assert n > 0
    for ii ∈ 1:1:n
      @inbounds x           = view( X, :, ii );
      @inbounds fitness[ii] = f( x; kwargs=kwargs );
    end
  end

  function drawCompetitorsPositionsBoundary!(;rng,
                                             X::Matrix{Float32},
                                             rand_mat::Matrix{Float32},
                                             upper::Vector{Float32},
                                             lower::Vector{Float32} )
    @assert length( upper ) == length( lower )
    @assert all( size( X )  == size( rand_mat ) )

    nFireworks = size( X, 2 );
    d          = length( lower );
    rand_mat  .= Float32.( rand( rng, d, nFireworks ) );
    for iter ∈ 1:1:nFireworks
      @inbounds x_sp     = view( X,        :, iter );
      @inbounds rand_vec = view( rand_mat, :, iter );
      for k ∈ 1:1:d
        @inbounds lu      = lower[k]::Float32;
        @inbounds x_sp[k] = ( rand_vec[k] * ( upper[k] - lu ) + lu )::Float32;
      end
    end
  end

  @views function  Δ_internal_forces!( S_i::Matrix{Float32},
                                       Δ_S_i::Matrix{Float32},
                                       fitnessS_i::Vector{Float32},
                                       rng::RandomDevice;
                                       diff_x::SubArray{Float32},
                                       rr::SubArray{Float32},
                                       en_n::Vector{Int} )
    n       = size( S_i, 2 )::Int;
    d       = size( S_i, 1 )::Int;
    max_val = maximum( fitnessS_i )::Float32;
    min_val = minimum( fitnessS_i )::Float32;
    rr     .= rand( rng, Float32, n )::Vector{Float32}; ## bad init size creates problems

    for ii ∈ 1:1:n
      @inbounds x_src    = view( S_i,   :, ii );
      @inbounds δ_x_src  = view( Δ_S_i, :, ii );
      @inbounds notVals  = view( en_n, Not( ii ) );
      for jj ∈ notVals
        @inbounds x_trgt = view( S_i, :, jj );
        diff_x          .= x_trgt .- x_src;
        denom            = norm( x_trgt .- x_src, 2 )::Float32;
        denom2           = ( ( max_val - min_val ) * denom )::Float32;
        if denom > 0.0f0 && abs( denom2 ) != 0.0f0
          @inbounds δ_x_src .+=  ( diff_x / denom ) * tanh( ( fitnessS_i[ii] - fitnessS_i[jj] ) / denom2  ) * rr[ii];
        end
      end
    end
  end

  function  Δ_external_forces!( S_i::Matrix{Float32},
                                X::Matrix{Float32},
                                Δ_S_i::Matrix{Float32},
                                fitnessX::Vector{Float32},
                                fitnessS_i::Vector{Float32},
                                diff_x::SubArray{Float32},
                                rr::SubArray{Float32},
                                rng::RandomDevice )

    n       = size( S_i, 2 )::Int;
    d       = size( S_i, 1 )::Int;
    n2      = size( X, 2 )::Int;
    max_val = maximum( fitnessS_i )::Float32;
    min_val = minimum( fitnessS_i )::Float32;
    rr     .= rand( rng, Float32, n )::Vector{Float32};

    for ii ∈ 1:1:n
      @inbounds x_src    = view( S_i,   :, ii );
      @inbounds δ_x_src  = view( Δ_S_i, :, ii );
      for jj ∈ 1:1:n2
        @inbounds x_trgt = view( X, :, jj );
        diff_x          .= x_trgt .- x_src;
        denom            = norm( x_trgt .- x_src, 2 )::Float32;
        denom2           = ( ( max_val - min_val ) * denom )::Float32;
        if denom > 0.0f0 && abs( denom2 ) != 0.0f0
          @inbounds δ_x_src .+=  ( diff_x / denom ) * tanh( ( fitnessS_i[ii] - fitnessS_i[jj] ) / denom2  ) * rr[ii];
        end
      end
    end

  end

  struct COMPSWARM
    X::Matrix{Float32};
    fitnessX::Vector{Float32};
    A::Vector{Float32};
    S::Vector{ Matrix{Float32} };
    fitnessS::Vector{ Vector{Float32} };
    x_b::Vector{Float32};
    y_min::Float32;
    iter::Int;
  end

  function allSolutions( cswarm::COMPSWARM )
    return cswarm.X, cswarm.fitnessX
  end

  function extractBoundary!( x_b::Vector{Float32};
                             lower::Vector{Float32},
                             upper::Vector{Float32},
                             ρ::Float32 )
    δ_x = ( upper .- lower ) * ρ / ( 2.0f0 )

    @assert all( δ_x .>= 0.0f0 )
    lower .= x_b .- δ_x;
    upper .= x_b .+ δ_x;
    @assert all( lower .<= upper )

  end

  @views function enforceBoundary!(;S_i::Matrix{Float32},
                                    lower::Vector{Float32},
                                    upper::Vector{Float32},
                                    rr::Vector{Float32},
                                    rng::RandomDevice )
    n  = size( S_i, 2 );
    d  = length( lower );
    for ii ∈ 1:1:n
      @inbounds s = view( S_i, :, ii );
      for jj ∈ 1:1:d
        @inbounds if s[jj] < lower[jj] || s[jj] > upper[jj]
          rr   .= rand( rng, Float32, d );
          @inbounds s[jj] = rr[jj] * ( upper[jj] - lower[jj] ) + lower[jj];
        end
      end
    end

  end

  """
           competetiveSwarm( nCompetitors::Int,
                             nSiblings::Int,
                             kwargs...;
                             λ_0::Float32,
                             ϵ_A::Float32,
                             C_a::Float32,
                             C_r::Float32,
                             lower::Vector{Float32},
                             upper::Vector{Float32},
                             objFunction::Function,
                             maxiter::Int,
                             global_forces::Bool=true,
                             ϵ_conv::Float32=0.001f0 )


      novel swarm intelligence algorithm, using fireworks algorithm for initialization.
      The method is based on a simplistic attraction/detraction mechanism:
        - each solution is attracted to solutions with better fitness,
          in analogue fashion each solution is detracted by solution with worse fitness.
        - the attraction/detraction depends on δ fitness as well as the euclidic distance between the two solutions.

      Competetive Swarm algorithm:
        - nCompetitors  number of (competing) clusters used in simulation
        - nSiblings     number of objFunction evaluations per cluster.
        - lower         lower limit solution space
        - upper         upper limit solution space
        - objFunction   objective function to be minimized
        - maxiter       nr of main algo iteration, each iteration triggers
                        nCompetitors * nSiblings objective function evaluations.
        - global_forces use intra-cluster forces in the same fashion as in-cluster forces.
        - ϵ_conv        convergence parameter not functional yet, but will be in future versions
      The Firework initialization is governed by the following parameters:
        - λ_0 maximum sum amplitudes of explosions
        - ϵ_A explosion amplitude smoothing parameter
        - C_a scaling factor up   explosion amplitude
        - C_r scaling factor down explosion amplitude

  """
  function competetiveSwarm( nCompetitors::Int,
                             nSiblings::Int,
                             kwargs...;
                             λ_0::Float32,
                             ϵ_A::Float32,
                             C_a::Float32,
                             C_r::Float32,
                             lower::Vector{Float32},
                             upper::Vector{Float32},
                             objFunction::Function,
                             maxiter::Int,
                             global_forces::Bool=false,
                             ϵ_conv::Float32=0.001f0 )

    @assert nCompetitors > 0
    @assert nSiblings    > 0
    @assert λ_0          > 0.0f0
    @assert ϵ_A          > 0.0f0
    @assert C_a          > 1.0f0
    @assert C_r          < 1.0f0
    @assert C_r          > 0.0f0
    @assert ( length( lower ) == length( upper ) )
    @assert all( lower .< upper )
    @assert maxiter      > 0
    @assert ϵ_conv       > 0.0f0

    d           = length( lower );
    rng         = RandomDevice();
    rr_boundary = Vector{Float32}( undef, d );
    A           = Vector{Float32}( undef, nCompetitors );
    cacheA      = Vector{Float32}( undef, nCompetitors );
    x_b         = Vector{Float32}( undef, d );
    rand_mat    = Matrix{Float32}( undef, d, nCompetitors );
    X           = Matrix{Float32}( undef, d, nCompetitors );
    fitnessX    = Vector{Float32}( undef, nCompetitors );
    S           = Vector{ Matrix{Float32} }( undef, nCompetitors );
    Δ_S_int     = Vector{ Matrix{Float32} }( undef, nCompetitors );
    Δ_S_ext     = Vector{ Matrix{Float32} }( undef, nCompetitors );
    fitnessS    = Vector{ Vector{Float32} }( undef, nCompetitors );
    for ii ∈ 1:1:nCompetitors
      @inbounds fitnessS[ii] = zeros( Float32, nSiblings );
      @inbounds S[ii]        = zeros( Float32, d, nSiblings );
      @inbounds Δ_S_int[ii]  = zeros( Float32, d, nSiblings );
      @inbounds Δ_S_ext[ii]  = zeros( Float32, d, nSiblings );
    end
    ## init randomly
    drawCompetitorsPositionsBoundary!( ;rng, X, rand_mat, upper, lower );
    evaluateCandidates( ;kwargs  = kwargs,
                         f       = objFunction,
                         X       = X,
                         fitness = fitnessX );
    cb    = argmin( fitnessX )::Int;
    @inbounds y_min = fitnessX[cb]::Float32;
    # compute explosion amplitudes
    explosionAmplitudes!( ;A, λ_0, ϵ_A, y_min, cacheA=cacheA, fitness_fireworks=fitnessX );
    @inbounds A[cb] = ( minimum( A[ Not( cb ) ] ) / 2.0f0 )::Float32;
    # use fw algorithm for initialization
    Threads.@threads for jj ∈ 1:1:nCompetitors
      explosionSparks(;X_i   = X[:,jj],
                       A_i   = A[jj],
                       s     = S[jj],
                       lower = lower,
                       upper = upper,
                       rng   = rng );
      evaluateCandidates(;kwargs   = kwargs,
                          f        = objFunction,
                          X        = S[jj],
                          fitness  = fitnessS[jj] );
    end
    rrX   = zeros( Float32, nSiblings, Threads.nthreads() );
    diffX = zeros( Float32,         d, Threads.nthreads() );
    en_n  = collect( 1:1:nCompetitors )::Vector{Int};
    iter  = 1;
    while iter < maxiter
      @inbounds A_cb_old = A[cb];
      explosionAmplitudes!( ;A, λ_0, ϵ_A, y_min, cacheA=cacheA, fitness_fireworks=fitnessX );
      @inbounds A[cb] = A_cb_old;
      @inbounds y_cb  = fitnessX[cb];
      dynExplosionAmplitude( y_min, y_cb, A_cb_old, C_a, C_r );
      Threads.@threads for ii ∈ 1:1:nCompetitors
        rr     = view( rrX,   :, Threads.threadid() );
        diff_x = view( diffX, :, Threads.threadid() );
        @inbounds Δ_internal_forces!( S[ii], Δ_S_int[ii], fitnessS[ii], rng;
                                      diff_x = diff_x,
                                      rr     = rr,
                                      en_n   = en_n );
        if global_forces
          Δ_external_forces!( S[ii], X, Δ_S_ext[ii], fitnessX, fitnessS[ii], diff_x, rr, rng );
        end
      end
      ## update, evaluate candidates S
      Threads.@threads for ii ∈ 1:1:nCompetitors
        for jj ∈ 1:1:nSiblings
          @inbounds S[ii][:,jj] .= X[:,ii] .+ A[ii] * ( Δ_S_int[ii][:,jj] .+ Δ_S_ext[ii][:,jj] );
        end
        enforceBoundary!(;S_i   = S[ii],
                          lower = lower,
                          upper = upper,
                          rr    = rr_boundary,
                          rng   = rng );
        ## evaluate updated candidates
        evaluateCandidates(;kwargs   = kwargs,
                            f        = objFunction,
                            X        = S[ii],
                            fitness  = fitnessS[ii] );
      end
      ## sequential communication phase ##
      for ii ∈ 1:1:nCompetitors
        @inbounds val = minimum( fitnessS[ii] );
        @inbounds idx = argmin( fitnessS[ii] );
        if val < y_min
          y_min = val;
          @inbounds x_b  .= S[ii][:,idx]
        end
        if val < fitnessX[ii]
          @inbounds fitnessX[ii] = val;
          @inbounds X[:,ii]     .= S[ii][:,idx]
        end
      end
      cb = argmin( fitnessX );
      iter += 1;
    end
    return COMPSWARM( X, fitnessX, A , S, fitnessS, x_b, y_min, iter )
  end

  """
           competetiveSwarmDynBoundary( nCompetitors::Int,
                                        nSiblings::Int,
                                        kwargs...;
                                        λ_0::Float32,
                                        ϵ_A::Float32,
                                        C_a::Float32,
                                        C_r::Float32,
                                        lower::Vector{Float32},
                                        upper::Vector{Float32},
                                        objFunction::Function,
                                        maxiter::Int,
                                        global_forces::Bool=false,
                                        ϵ_conv::Float32=0.001f0,
                                        nRefinements::Int=4 )

    Same as competetiveSwarm, but with dynamic iterative refinement of the search space,
    yielding improved results. nRefinements denotes the number of search space refinements applied.

  """
  function competetiveSwarmDynBoundary( nCompetitors::Int,
                                        nSiblings::Int,
                                        kwargs...;
                                        λ_0::Float32,
                                        ϵ_A::Float32,
                                        C_a::Float32,
                                        C_r::Float32,
                                        lower::Vector{Float32},
                                        upper::Vector{Float32},
                                        obj_f::Function,
                                        maxiter::Int,
                                        global_forces::Bool= true,
                                        ϵ_conv::Float32    = 0.01f0,
                                        nRefinements::Int  = 4,
                                        ρ::Float32         = 0.25f0 )
     lo = deepcopy( lower )
     up = deepcopy( upper )

     resSwarm = [competetiveSwarm( nCompetitors, nSiblings, kwargs...;
                                  λ_0, ϵ_A, C_a, C_r, lower=lo, upper=up,
                                  objFunction=obj_f, maxiter, global_forces,
                                  ϵ_conv )];
    extractBoundary!( resSwarm[1].x_b;
                      lower = lo,
                      upper = up,
                      ρ     = ρ )
    iter = 1;
    while iter < nRefinements
      resSwarm .= [competetiveSwarm( nCompetitors, nSiblings, kwargs...;
                                   λ_0, ϵ_A, C_a, C_r, lower=lo, upper=up,
                                   objFunction=obj_f, maxiter, global_forces,
                                   ϵ_conv )];
      extractBoundary!( resSwarm[1].x_b;
                        lower = lo,
                        upper = up,
                        ρ     = ρ )
      iter += 1;
    end
    return resSwarm[1]
  end

  function eval_benchmark_2D_Dyn_Bound( obj_f::Function, kwargs...;
                                        nCompetitors::Int   = 32,
                                        nSiblings::Int      = 32,
                                        λ_0::Float32        = 4.0f0,
                                        ϵ_A::Float32        = 0.01f0,
                                        C_a::Float32        = 1.1f0,
                                        C_r::Float32        = 0.9f0,
                                        maxiter::Int        = 20,
                                        global_forces::Bool = false,
                                        ϵ_conv::Float32     = 1f-6,
                                        nRefinements::Int   = 4  )

   @assert nCompetitors > 0
   @assert nSiblings    > 0
   @assert λ_0          > 0.0f0
   @assert ϵ_A          > 0.0f0
   @assert C_a          > 1.0f0
   @assert C_r          > 0.0f0
   @assert C_r          < 1.0f0
   @assert maxiter      > 0

   lower, upper          = retrieveLimits2D( obj_f );
   cSwarm( objFunction ) = competetiveSwarmDynBoundary( nCompetitors, nSiblings, kwargs;
                                                        λ_0          = λ_0,
                                                        ϵ_A          = ϵ_A,
                                                        C_a          = C_a,
                                                        C_r          = C_r,
                                                        lower        = lower,
                                                        upper        = upper,
                                                        obj_f        = obj_f,
                                                        ϵ_conv       = ϵ_conv,
                                                        maxiter      = maxiter,
                                                        nRefinements = nRefinements )
    return cSwarm( obj_f )
  end

  function eval_benchmark_2D( obj_f::Function, kwargs...;
                              nCompetitors::Int = 32,
                              nSiblings::Int    = 32,
                              λ_0::Float32      = 4.0f0,
                              ϵ_A::Float32      = 0.01f0,
                              C_a::Float32      = 1.1f0,
                              C_r::Float32      = 0.9f0,
                              maxiter::Int      = 20,
                              ϵ_conv::Float32   = 1f-6,
                              global_forces     = false )

   @assert nCompetitors > 0
   @assert nSiblings    > 0
   @assert λ_0          > 0.0f0
   @assert ϵ_A          > 0.0f0
   @assert C_a          > 1.0f0
   @assert C_r          > 0.0f0
   @assert C_r          < 1.0f0
   @assert maxiter      > 0

   lower, upper          = retrieveLimits2D( obj_f );
   cSwarm( objFunction ) = competetiveSwarm( nCompetitors, nSiblings, kwargs;
                                             λ_0          = λ_0,
                                             ϵ_A          = ϵ_A,
                                             C_a          = C_a,
                                             C_r          = C_r,
                                             lower        = lower,
                                             upper        = upper,
                                             objFunction  = obj_f,
                                             ϵ_conv       = ϵ_conv,
                                             maxiter      = maxiter )
    println( "running competetiveSwarm" )
    println()
    println()
    return cSwarm( obj_f )
  end

  function eval_benchmark_ND( obj_f::Function, kwargs...;
                              d::Int            = 8,
                              nCompetitors::Int = 16,
                              nSiblings::Int    = 16,
                              λ_0::Float32      = 4.0f0,
                              ϵ_A::Float32      = 0.01f0,
                              C_a::Float32      = 1.1f0,
                              C_r::Float32      = 0.9f0,
                              maxiter::Int      = 64 )

   @assert nCompetitors > 0
   @assert nSiblings    > 0
   @assert λ_0          > 0.0f0
   @assert ϵ_A          > 0.0f0
   @assert C_a          > 1.0f0
   @assert C_r          > 0.0f0
   @assert C_r          < 1.0f0
   @assert maxiter      > 0

   lower, upper          = retrieveLimitsND( obj_f, d );
   cSwarm( objFunction ) = competetiveSwarm( nCompetitors, nSiblings, kwargs;
                                             λ_0         = λ_0,
                                             ϵ_A         = ϵ_A,
                                             C_a         = C_a,
                                             C_r         = C_r,
                                             lower       = lower,
                                             upper       = upper,
                                             objFunction = objFunction,
                                             maxiter     = maxiter )
    return cSwarm( obj_f )
  end

  function metaObjCSwarm2D_mod0( x, args... )
    nCompetitors = round( Int, x[1] )::Int;
    nSiblings    = round( Int, x[2] )::Int;
    λ_0          = x[3]::Float32;
    ϵ_A          = x[4]::Float32;
    C_a          = x[5]::Float32;
    C_r          = x[6]::Float32;
    maxiter      = round( Int, args[2] )::Int;
    fObj         = args[1]::Function;
    nSamples     = args[3]::Int;
    perf         = zeros( Float32, nSamples );
    for ii in 1:1:nSamples
      perf[ii] = eval_benchmark_2D( fObj ).y_min;
    end
    μ_p = median( perf )::Float32;
    return μ_p + ( 10.0f0 * max( 0.0f0,  ( x[1] * x[2] ) - args[4] ) ) ^ 4
  end

  function metaObjCSwarm2D_mod1( x, args... )
    nCompetitors = round( Int, x[1] )::Int;
    nSiblings    = round( Int, x[2] )::Int;
    λ_0          = x[3]::Float32;
    ϵ_A          = x[4]::Float32;
    C_a          = x[5]::Float32;
    C_r          = x[6]::Float32;

    μ_p = mean( [ eval_benchmark_2D( fObj ).y_min for fObj ∈ obj_functions_2D ] );

    return μ_p + ( 10.0f0 * max( 0.0f0,  ( x[1] * x[2] ) - args[1] ) ) ^ 4
  end

  function metaObjCSwarm2D_mod2( x, args... )
    nCompetitors = round( Int, x[1] )::Int;
    nSiblings    = round( Int, x[2] )::Int;
    λ_0          = x[3]::Float32;
    ϵ_A          = x[4]::Float32;
    C_a          = x[5]::Float32;
    C_r          = x[6]::Float32;

    μ_p = mean( [ eval_benchmark_2D( fObj, global_forces=true ).y_min for fObj ∈ obj_functions_2D ] );

    return μ_p + ( 10.0f0 * max( 0.0f0,  ( x[1] * x[2] ) - args[1] ) ) ^ 4
  end
end
