include("CannonLib.jl")
using .CannonLib
using Profile

function dothething(i)
  positions, velocities, addresses = recursive_bounces(0e0, 1e0, [0x1:0x64 for n in 1:i])
  #positions = getindex.(positions, 2)

  indices, heights, deltas = best_alignment("movable", "ender_pearl", positions)
  println("Top 100 Pearls[$i]")
  for i ∈ indices[1:100]
    println("$(addresses[i]); pos: $(positions[i]), vel: $(velocities[i]) block: $(heights[i]), delta: $(deltas[i])")
  end
  println()

  println("Top 100 Players[$i]")
  indices, heights, deltas = best_alignment("movable", "player", positions)
  for i ∈ indices[1:100]
    println("$(addresses[i]); pos: $(positions[i]), vel: $(velocities[i]) block: $(heights[i]), delta: $(deltas[i])")
  end
  println()
end

for i in 1:3
  @profile @time dothething(i)
end
