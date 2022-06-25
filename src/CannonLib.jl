"A Julia library for Minecraft simulations useful for making TNT cannons"
module CannonLib

include("Entity.jl")
import .Entity.tickprojectileyposlist
import .Entity.projectile_pos_y_table
export cannonangles, closestin, best_alignment, recursive_bounces, recursive_bounces_lookup

function projectile_pos_y_table_f(range, d)
  out::Vector{Float64} = []
  for i ∈ range
    push!(out, projectile_pos_y_table[i] + d)
  end
  out
end

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

function closestin(y::Float64, options::Vector{Float64})
    d = y - floor(y) .- options
    i::Int64 = partialsortperm(d, 1; by = abs)
    (i, d[i])
end

"""
    best_alignment(options::Vector{Float64}, eyeHeight::Float64, heights::Vector{Float64}, count::Int64)

    Scans over a list of heights and returns the closest [count] TNT positions, where the TNT is standing on a block of size [options].
"""
function best_alignment(options::Vector{Float64}, eyeHeight::Float32, heights::Vector{Float64})
  deltas::Vector{Float64} = []
  indices::Vector{Int64} = []
  for y ∈ heights .+ eyeHeight .- Float64(0.98f0*0.0625e0)
    i, d = closestin(y, options)
    push!(deltas, d)
    push!(indices, i)
  end

  p = sortperm(deltas; by = abs, alg=QuickSort)
  range = 1:length(p)
  (indices=p, blockheights=options[indices[range]], deltas=deltas[range])
end

"""
    best_alignment(blocktype::String, entitytype::String, heights::Vector{Float64}, count::Int64)

    Scans over a list of heights and returns the closest [count] TNT positions, where the TNT is standing on a block ∈ the category [blocktype].
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

  if blocktype == "all" return best_alignment([0,1,1.5,2,3,4,5,6,7,8,9,9.5,10,11,12,13,14,15,16]./16e0, eyeheight, heights)
  if blocktype == "movable_0ticks" return best_alignment([0,1,2,3,4,6,8,9,9.5,10,12,14,15,16]./16e0, eyeheight, heights)
  if blocktype == "movable" return best_alignment([0,1,3,8,9,9.5,10,14,15,16]./16e0, eyeheight, heights)
  else return [] end
end

function recursive_bounces(options::Vector{Float64}, pos::Float64, vel::Float64, tickranges::Vector{UnitRange{UInt8}}, limit::Float64, eyeheight::Float32; explosionheight::Float64 = Float64(0.98f0*0.0625e0), addr::Vector{UInt8} = UInt8[])
  positions::Vector{Float64} = tickprojectileyposlist(pos, vel, tickranges[1])
  posout::Vector{Float64} = []
  addressout::Vector{Vector{UInt8}} = []
  indexout::Vector{UInt8} = []
  deltaout::Vector{Float64} = []
  firstaddress::UInt8 = tickranges[1][1] - 0x1
  l::UInt8 = length(tickranges)
  if l > 0x1
    for i::UInt8 ∈ eachindex(positions)
      pos = positions[i]
      range = tickranges[0x2:l]
      range[0x1] = range[0x1][firstaddress + i:length(range[0x1])]
      bounces = recursive_bounces(options, pos, 1e0, range, limit, eyeheight; explosionheight=explosionheight, addr=vcat(addr, firstaddress + i))
      if length(bounces[:addr]) > 0x0
        addressout = vcat(addressout, bounces[:addr])
        posout = vcat(posout, bounces[:pos])
        indexout = vcat(indexout, bounces[:index])
        deltaout = vcat(deltaout, bounces[:delta])
      end
    end
    return (addr=addressout, pos=posout, index=indexout, delta=deltaout)
  end

  for j in 0x0:length(addr) + 0x1
    for i in eachindex(positions)
      index, d = closestin(positions[i] + eyeheight - explosionheight, options)
      if abs(d) < limit
        push!(addressout, vcat(addr, firstaddress + i, j))
        push!(posout, positions[i])
        push!(indexout, index)
        push!(deltaout, d)
      end
      positions[i] += 0.51e0
    end
  end
  return (addr=addressout, pos=posout, index=indexout, delta=deltaout)
end

end
