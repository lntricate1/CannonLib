"A Julia library for simulating Minecraft entities"
module Entity

export tick_projectile, tick_projectile_list, tick_projectile_y, tick_projectile_y_honey
export tick_tnt, tick_tnt_list, tick_tnt_y

"""
    tick_projectile(pos, vel, ticks)

Accurate projectile position after `ticks` ticks.

Returns a NamedTuple of (pos::Vector{Float64}, vel::Vector{Float64})
"""
function tick_projectile(pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int=1)::NamedTuple{(:pos, :vel), Tuple{Vector{Float64}, Vector{Float64}}}
  g = Float32[0, -0.03, 0]
  for _ ∈ 1:ticks
    pos += vel
    vel *= 0.99f0
    vel += g
  end
  (pos=pos, vel=vel)
end

"""
    tick_projectile_list(pos, vel, ticks)

List of accurate positions of a projectile until `ticks`.

Returns a NamedTuple of Vectors in the form (pos::Vector{Vector{Float64}}, vel::Vector{Vector{Float64}})
"""
function tick_projectile_list(pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int64)::NamedTuple{(:pos, :vel), Tuple{Vector{Vector{Float64}}, Vector{Vector{Float64}}}}
  posout = Vector{Float64}[]
  velout = Vector{Float64}[]
  g = Float32[0, -0.03, 0]
  for _ ∈ 1:ticks
    pos += vel
    vel *= 0.99f0
    vel += g
    push!(posout, pos);
    push!(velout, vel);
  end
  (pos=posout, vel=velout)
end

"""
    tick_projectile_y(pos, vel, ticks)

Accurate projectile Y position after `ticks` ticks.

Returns a NamedTuple of (pos::Float64, vel::Float64)
"""
function tick_projectile_y(pos::Float64, vel::Float64, ticks::Int=1)::Float64
  for _ ∈ 1:ticks
    pos += vel
    vel *= 0.99f0
    vel -= 0.03f0
  end
  (pos=pos, vel=vel)
end

"""
    tick_projectile_y_honey(pos, vel, ticks)

Accurate projectile Y position after `ticks` ticks in honey.

Returns a NamedTuple of (pos::Float64, vel::Float64)
"""
function tick_projectile_y_honey(pos::Float64, vel::Float64, ticks::Int=1)::Float64
  for _ ∈ 1:ticks
    pos += vel
    vel *= 0.99f0
    vel -= 0.03f0
    if vel < -0.08e0 vel = -0.05e0 end
  end
  (pos=pos, vel=vel)
end

"""
    tick_tnt(pos, vel, ticks)

Accurate tnt position after `ticks` ticks.

Returns a NamedTuple of (pos::Vector{Float64}, vel::Vector{Float64})
"""
function tick_tnt(pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int=1)::NamedTuple{(:pos, :vel), Tuple{Vector{Float64}, Vector{Float64}}}
  g = Float64[0, -0.04, 0]
  for _ ∈ 1:ticks
    vel += g
    pos += vel
    vel *= 0.98e0
  end
  (pos=pos, vel=vel)
end

"""
    tick_tnt_list(pos, vel, ticks)

List of accurate positions of a tnt until `ticks`.

Returns a NamedTuple of Vectors in the form (pos::Vector{Vector{Float64}}, vel::Vector{Vector{Float64}})
"""
function tick_tnt_list(pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int64)::NamedTuple{(:pos, :vel), Tuple{Vector{Vector{Float64}}, Vector{Vector{Float64}}}}
  posout = Vector{Float64}[]
  velout = Vector{Float64}[]
  g = Float64[0, -0.04, 0]
  for _ ∈ 1:ticks
    vel += g
    pos += vel
    vel *= 0.98e0
    push!(posout, pos);
    push!(velout, vel);
  end
  (pos=posout, vel=velout)
end

"""
    tick_tnt_y(pos, vel, ticks)

Accurate tnt Y position after `ticks` ticks.

Returns a NamedTuple of (pos::Float64, vel::Float64)
"""
function tick_tnt_y(pos::Float64, vel::Float64, ticks::Int=1)::NamedTuple{(:pos, :vel), Tuple{Float64, Float64}}
  for _ ∈ 1:ticks
    vel -= 0.04e0
    pos += vel
    vel *= 0.98e0
  end
  (pos=pos, vel=vel)
end

end
