function retrieveMin2D( f::Function )
  idx = findall(x -> x==f, obj_functions_2D )
  @assert length( idx ) > 0
  return f( min_functions_2D[idx][1](2), kwargs=() )
end

function retrieveMinND( f::Function, d::Int )
  idx = findall(x -> x==f, obj_functions_ND )
  @assert length( idx ) > 0
  return f( min_functions_ND[idx][1](d), kwargs=() )
end

function retrieveLimits2D( f::Function )
  idx = findall(x -> x==f, obj_functions_2D )
  print( "retrieving limits for objective function: " )
  show( f )
  println()
  return lim_functions_2D[idx][1]( 2 )
end

function retrieveLimitsND( f::Function, d::Int )
  idx = findall(x -> x==f, obj_functions_ND )
  @assert length( idx ) > 0
  return lim_functions_ND[idx][1]( d )
end
