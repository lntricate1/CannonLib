include("CannonLib.jl")
using .CannonLib
using Profile

function dothething(i)
  @time positions::Vector{Float64}, addresses::Vector{Vector{Int64}} = recursive_bounces([0,1,3,8,9,9.5,10,14,15,16]./16e0, 0e0, 1e0, [0x1:0x64 for n in 1:i], 1e-5)
  #positions = getindex.(positions, 2)

  println("Top 100 Pearls[$i]")
  @time indices, heights, deltas = best_alignment("movable", "ender_pearl", positions)
  for i ∈ indices[1:100]
    println("$(addresses[i]); pos: $(positions[i]), block: $(heights[i]), delta: $(deltas[i])")
  end
  println()

  println("Top 100 Players[$i]")
  @time indices, heights, deltas = best_alignment("movable", "player", positions)
  for i ∈ indices[1:100]
    println("$(addresses[i]); pos: $(positions[i]), block: $(heights[i]), delta: $(deltas[i])")
  end
  println()
end

for i in 6:10
@profile @time dothething(i)
end
