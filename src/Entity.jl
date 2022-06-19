"A Julia library for simulating Minecraft entities"
module Entity

export simprojectile, tickprojectile, tickprojectilelist, tickprojectile!,
       simtnt, ticktnt, ticktnt!,
       projectilecalcvel

"""
    projectilecalcvel(dpos::Vector, ticks)

Initial velocity of a projectile given a change in position and number of ticks.
"""
function projectilecalcvel(dpos::Vector, ticks)
  ticks < 0 && throw(DomainError(ticks, "Ticks must be positive"))
  D = 0.99
  g = [0, -0.03, 0]
  s = (D^ticks - 1)/(D - 1)
  (dpos + g*(s-ticks)/(D-1))/s
end

"""
    projectilecalcpos(vel::Vector, ticks)

Change in position of a projectile given an initial velocity and number of ticks (which can be negative).
"""
function projectilecalcpos(vel::Vector, ticks)
  D = 0.99
  g = [0, -0.03, 0]
  s = (D^ticks - 1)/(D - 1)
  vel * s + g * 1/(D - 1) * (s - ticks)
end

"""
    simprojectile(pos, vel, ticks)

State of a projectile in a different time, given an initial state. Ticks can be negative to look into the past.

Returns a tuple of (pos::Vector, vel::Vector)

WARNING: This is meant for optimal performance at high values of abs(ticks). It WILL have floating point rounding errors of increasing severity at bigger time frames, because Minecraft has rounding errors that this does not. For a simulation 100% accurate to the game, use [`tickprojectile`](@ref).
"""
function simprojectile(pos, vel, ticks)
  D = 0.99
  g = [0, -0.03, 0]
  s = (D^ticks - 1)/(D - 1)
  pos1 = pos + vel * s +
    g * 1/(D - 1) * (s - ticks)
  vel1 = vel * D^ticks + g * s
  (pos=pos1, vel=vel1)
end

"""
    tickprojectile(pos, vel, ticks)

Accurate state of a projectile in a different time, given an initial state. Ticks can be negative to look into the past.

Returns a tuple of (pos::Vector{Float64}, vel::Vector{Float64})
"""
function tickprojectile(pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int=1)::NamedTuple{(:pos, :vel), Tuple{Vector{Float64}, Vector{Float64}}}
  pos1 = pos
  vel1 = vel
  D = 0.99f0
  g = Float32[0, -0.03, 0]
  if ticks > 0
    for i ∈ 1:ticks
      pos1 += vel1
      vel1 *= D
      vel1 += g
    end
  else
    for i ∈ 1:(-ticks)
      vel1 -= g
      vel1 /= D
      pos1 -= vel1
    end
  end
  (pos=pos1, vel=vel1)
end

"""
    tickprojectilelist(pos, vel, ticks)

Accurate state of a projectile in a different time, given an initial state.

Returns a NamedTuple of Vectors in the form (pos::Vector{Vector{Float64}}, vel::Vector{Vector{Float64}})
"""
function tickprojectilelist(pos::Vector{Float64}, vel::Vector{Float64}, ticks::UnitRange)::NamedTuple{(:pos, :vel), Tuple{Vector{Vector{Float64}}, Vector{Vector{Float64}}}}
  pos1::Vector{Float64} = pos
  vel1::Vector{Float64} = vel
  posout::Vector{Float64} = []
  velout::Vector{Float64} = []
  for Nothing ∈ 1:(ticks[1]+length(ticks)-1)
    pos1 += vel1
    vel1 *= 0.99f0
    vel1 += [0f0, -0.03f0, 0f0]
    push!(posout, pos1);
    push!(velout, vel1);
  end
  (pos=posout[ticks], vel=velout[ticks])
end

"""
    tickprojectilelist(pos, vel, ticks)

Accurate state of a projectile in a different time, given an initial state.

Returns a NamedTuple of Vectors in the form (pos::Vector{Vector{Float64}}, vel::Vector{Vector{Float64}})
"""
function tickprojectilelist(pos::Float64, vel::Float64, ticks::UnitRange)::NamedTuple{(:pos, :vel), Tuple{Vector{Float64}, Vector{Float64}}}
  pos1::Float64 = pos
  vel1::Float64 = vel
  posout::Vector{Float64} = []
  velout::Vector{Float64} = []
  for Nothing ∈ 1:(ticks[1]+length(ticks)-1)
    pos1 += vel1
    vel1 = vel1 * 0.99f0 -0.03f0
    push!(posout, pos1);
    push!(velout, vel1);
  end
  (pos=posout[ticks], vel=velout[ticks])
end

"""
    tickprojectile!(pos, vel, ticks)

Accurate state of a projectile in a different time, given an initial state. Ticks can be negative to look into the past.

Returns a tuple of (pos::Vector{Float64}, vel::Vector{Float64})
"""
function tickprojectile!(pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int=1)::NamedTuple{(:pos, :vel), Tuple{Vector{Float64}, Vector{Float64}}}
  D = 0.99f0
  g = Float32[0, -0.03, 0]
  if ticks > 0
    for i ∈ 1:ticks
      pos += vel
      vel *= D
      vel += g
    end
  else
    for i ∈ 1:(-ticks)
      vel -= g
      vel /= D
      pos -= vel
    end
  end
  (pos=pos, vel=vel)
end

"""
    simtnt(pos, vel, ticks)

State of a tnt entity in a different time, given an initial state. Ticks can be negative to look into the past.

Returns a tuple of (pos::Vector, vel::Vector)

WARNING: This is meant for optimal performance at high values of abs(ticks). It WILL have floating point rounding errors of increasing severity at bigger time frames, because Minecraft has rounding errors that this does not. For a simulation 100% accurate to the game, use [`ticktnt`](@ref).
"""
function simtnt(pos, vel, ticks)
  D = 0.98
  g = [0, -0.04, 0]
  s = (D^ticks - 1)/(D - 1)
  pos1 = pos + vel * s +
    g * (D/(D - 1) * (s - ticks) + ticks)
  vel1 = vel * D^ticks + g * D*s
  (pos=pos1, vel=vel1)
end

"""
    ticktnt(pos, vel, ticks)

Accurate state of a tnt entity in a different time, given an initial state. Ticks can be negative to look into the past.

Returns a tuple of (pos::Vector{Float64}, vel::Vector{Float64})
"""
function ticktnt(pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int=1)::Tuple{Vector{Float64}, Vector{Float64}}
  pos1 = pos
  vel1 = vel
  D = 0.98e0
  g = Float64[0, -0.04, 0]
  if ticks > 0
    for i ∈ 1:ticks
      vel1 += g
      pos1 += vel1
      vel1 *= D
    end
  else
    for i ∈ 1:(-ticks)
      vel1 -= g
      pos1 -= vel1
      vel1 /= D
    end
  end
  (pos=pos1, vel=vel1)
end

"""
    ticktnt!(pos, vel, ticks)

Accurate state of a tnt entity in a different time, given an initial state. Ticks can be negative to look into the past.

Returns a tuple of (pos::Vector{Float64}, vel::Vector{Float64})
"""
function ticktnt!(pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int=1)::Tuple{Vector{Float64}, Vector{Float64}}
  D = 0.98e0
  g = Float64[0, -0.04, 0]
  if ticks > 0
    for i ∈ 1:ticks
      vel += g
      pos += vel
      vel *= D
    end
  else
    for i ∈ 1:(-ticks)
      vel -= g
      pos -= vel
      vel /= D
    end
  end
  (pos=pos, vel=vel)
end


end
