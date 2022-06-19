"A Julia library for Minecraft simulations useful for making TNT cannons"
module CannonLib

include("Entity.jl")
import .Entity.tickprojectilelist
export cannonangles, best_alignment, recursive_bounces

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
function best_alignment(options::Vector{Float64}, eyeHeight::Float32, heights::Vector{Float64})
  deltas = []
  indices = []
  e = Float64(0.98f0*0.0625e0)
  heights_ = heights .+ eyeHeight .- e
  for y in heights_
    d = y - floor(y) .- options
    i::Int64 = partialsortperm(d, 1; by = abs)
    push!(deltas, d[i])
    push!(indices, i)
  end

  p = sortperm(deltas; by = abs, alg=QuickSort)
  range = 1:length(p)
  (indices=p, blockheights=options[indices[range]], deltas=deltas[range])
end

"""
    best_alignment(blocktype::String, entitytype::String, heights::Vector{Float64}, count::Int64)

    Scans over a list of heights and returns the closest [count] TNT positions, where the TNT is standing on a block in the category [blocktype].
"""
function best_alignment(blocktype::String, entitytype::String, heights::Vector{Float64})
  eyeheight = 0.25f0*0.85f0
  if entitytype != "ender_pearl" && entitytype != "snowball" && entitytype != "item" && entitytype != "fishing_bobber"
    if entitytype == "tnt" eyeheight = 0f0
    elseif entitytype == "arrow" eyeheight = 0.13f0
    elseif entitytype == "player" eyeheight = 0.4f0
    elseif entitytype == "falling_block" eyeheight = 0.98f0*0.85f0
      else throw(ArgumentError("Invalid entity type")) end
  end

  if blocktype == "movable" return best_alignment([0,1,3,8,9,9.5,10,14,15,16]./16e0, eyeheight, heights)
  else return [] end
end

function recursive_bounces(pos::Vector{Float64}, vel::Vector{Float64}, tickranges::Vector{UnitRange{Int64}})
  ticks = tickprojectilelist(pos, vel, tickranges[1][length(tickranges[1])])
  positions = ticks[:pos][tickranges[1]]
  velout = ticks[:vel][tickranges[1]]
  if length(tickranges) > 1
    posout = Vector{Float64}[]
    velout = Vector{Float64}[]
    addressout = Vector{Float64}[]
    i = 1
    for pos in positions
      bounces = recursive_bounces(pos, vel, tickranges[2:length(tickranges)])
      for pos in bounces[:pos]
        push!(posout, pos)
      end
      for vel in bounces[:vel]
        push!(velout, vel)
      end
      for addr in bounces[:addr]
        push!(addressout, [i, addr...])
      end
      i += 1
    end
    return (pos=posout, vel=velout, addr=addressout)
  end
  return (pos=positions, vel=velout, addr=tickranges[1])
end

end
