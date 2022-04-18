using CompetetiveSwarmOptim
using Test

@testset "CompetetiveSwarmOptim.jl" begin
   # 2D-Tests 
   result_2D = sum( [ eval_benchmark_2D( fObj ).y_min for fObj ∈ obj_functions_2D ] );
   min_2D    = sum( [ retrieveMin2D( fObj ) for fObj ∈ obj_functions_2D ] );

   @test result_2D <= min_2D + 2.0f0

   result_2D_dyn = sum( [ eval_benchmark_2D_Dyn_Bound( fObj ).y_min for fObj ∈ obj_functions_2D ] );
   @test result_2D_dyn <= min_2D + 2.0f0

end
