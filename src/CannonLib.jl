"A Julia library for Minecraft simulations useful for making TNT cannons"
module CannonLib

export cannonangles, best_alignment

"""
    cannonangles(tnt::Int)

Number of possible angles given a maximum TNT count for 1 side of the cannon.
"""
function cannonangles(tnt::Int)
  count::Int = 0
  for x::Int ∈ 2:(tnt-1)
    for z::Int ∈ x:tnt
      gcd(x, z) == 1 && (count += 1) end end
  count*8 + tnt*8
end

"""
    best_alignment(options::Vector{Float64}, eyeHeight::Float64, heights::Vector{Float64}, count::Int64)

    Scans over a list of heights and returns the closest [count] TNT positions, where the TNT is standing on a block of size [options].
"""
function best_alignment(options::Vector{Float64}, eyeHeight::Float32, heights::Vector{Float64}, count::Int64)
  diffs = []
  indexes = []
  e = Float64(0.98f0*0.0625e0)
  heights_ = heights .+ eyeHeight .- e
  for y in heights_
    d = y - floor(y) .- options
    i::Int64 = partialsortperm(d, 1; by = abs)
    push!(diffs, d[i])
    push!(indexes, i)
  end

  out = []
  j::Int64 = 0
  for i::Int64 in sortperm(diffs; by = abs, alg=QuickSort)
    j += 1
    if j > count return out end
    push!(out, [i, heights[i], options[indexes[i]], diffs[i]]) end
  out
end

"""
    best_alignment(blocktype::String, entitytype::String, heights::Vector{Float64}, count::Int64)

    Scans over a list of heights and returns the closest [count] TNT positions, where the TNT is standing on a block in the category [blocktype].
"""
function best_alignment(blocktype::String, entitytype::String, heights::Vector{Float64}, count::Int64)
  eyeheight = 0.25f0*0.85f0
  if entitytype != "ender_pearl" && entitytype != "snowball" && entitytype != "item" && entitytype != "fishing_bobber"
    if entitytype == "tnt" eyeheight = 0f0
    elseif entitytype == "arrow" eyeheight = 0.13f0
    elseif entitytype == "player" eyeheight = 0.4f0
    elseif entitytype == "falling_block" eyeheight = 0.98f0*0.85f0
      else throw(ArgumentError("Invalid entity type")) end
  end

  if blocktype == "movable" return best_alignment([0,1,3,8,9,9.5,10,14,15,16]./16e0, eyeheight, heights, count)
  else return [] end
end

end
