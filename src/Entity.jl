export ThrownEntity, PersistentProjectileEntity, TNTEntity
export tick, tick_list, tick_y, tick_y_list, tick_y_honey

abstract type Entity end

struct ThrownEntity <: Entity end

struct PersistentProjectileEntity <: Entity end

struct TNTEntity <: Entity end

"""
    tick(entitytype, pos, vel, ticks; NoGravity=false)

Accurate position of an entity of type `entitytype` after `ticks` ticks.

Returns a NamedTuple of (pos::Vector{Float64}, vel::Vector{Float64})
"""
function tick(::Type{ThrownEntity}, pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int=1; NoGravity::Bool=false)
  g = Float32[0, -0.03, 0]
  for _ ∈ 1:ticks
    if sum(vel.^2) > 1e-7
      pos += vel
    end
    vel *= 0.99f0
    if !NoGravity vel += g end
  end
  (pos=pos, vel=vel)
end

"""
    tick_list(entitytype, pos, vel, ticks; NoGravity=false)

List of accurate positions of an entity of type `entitytype` until `ticks`.

Returns a NamedTuple of Vectors in the form (pos::Vector{Vector{Float64}}, vel::Vector{Vector{Float64}})
"""
function tick_list(::Type{ThrownEntity}, pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int64; NoGravity::Bool=false)
  posout = Vector{Float64}[]
  velout = Vector{Float64}[]
  g = Float32[0, -0.03, 0]
  for _ ∈ 1:ticks
    if sum(vel.^2) > 1e-7
      pos += vel
    end
    vel *= 0.99f0
    if !NoGravity vel += g end
    push!(posout, pos);
    push!(velout, vel);
  end
  (pos=posout, vel=velout)
end

"""
    tick_y(entitytype, pos, vel, ticks)

Accurate Y position of an entity of type `entitytype` after `ticks` ticks.

Returns a NamedTuple of (pos::Float64, vel::Float64)
"""
function tick_y(::Type{ThrownEntity}, pos::Float64, vel::Float64, ticks::Int)
  for _ ∈ 1:ticks
    pos += vel
    vel *= Float64(0.99f0)
    vel -= Float64(0.03f0)
  end
  (pos=pos, vel=vel)
end

"""
    tick_y(entitytype, pos, vel)

Accurate Y position of an entity of type `entitytype` in the next tick.

Returns a NamedTuple of (pos::Float64, vel::Float64)
"""
function tick_y(::Type{ThrownEntity}, pos::Float64, vel::Float64)
  pos += vel
  vel *= 0.99f0
  vel -= 0.03f0
  (pos=pos, vel=vel)
end

"""
    tick_y_list(entitytype, pos, vel, ticks)

List of accurate Y positions of an entity of type `entitytype` until `ticks`.

Returns a NamedTuple of (pos::Vector{Float64}, vel::Vector{Float64})
"""
function tick_y_list(::Type{ThrownEntity}, pos::Float64, vel::Float64, ticks::Int64)
  posout = Float64[]
  velout = Float64[]
  for _ ∈ 1:ticks
    pos += vel
    vel *= 0.99f0
    vel -= 0.03f0
    push!(posout, pos);
    push!(velout, vel);
  end
  (pos=posout, vel=velout)
end

"""
    tick_y_honey(entitytype, pos, vel, ticks)

Accurate Y position of an entity of type `entitytype` after `ticks` ticks in honey.

Returns a NamedTuple of (pos::Float64, vel::Float64)
"""
function tick_y_honey(::Type{ThrownEntity}, pos::Float64, vel::Float64, ticks::Int=1)
  for _ ∈ 1:ticks
    pos += vel
    vel *= 0.99f0
    vel -= 0.03f0
    if vel < -0.08e0 vel = -0.05e0 end
  end
  (pos=pos, vel=vel)
end

function tick(::Type{PersistentProjectileEntity}, pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int64=1; NoGravity::Bool=false)
  g = Float32[0, -0.05, 0]
  for _ ∈ 1:ticks
    if sum(vel.^2) > 1e-7
      pos += vel
    end
    vel *= 0.99f0
    if !NoGravity vel += g end
  end
  (pos=pos, vel=vel)
end

function tick_list(::Type{PersistentProjectileEntity}, pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int64=1; NoGravity::Bool=false)
  posout = Vector{Float64}[]
  velout = Vector{Float64}[]
  g = Float32[0, -0.05, 0]
  for _ ∈ 1:ticks
    if sum(vel.^2) > 1e-7
      pos += vel
    end
    vel *= 0.99f0
    if !NoGravity vel += g end
    push!(posout, pos);
    push!(velout, vel);
  end
  (pos=posout, vel=velout)
end

function tick_y(::Type{PersistentProjectileEntity}, pos::Float64, vel::Float64, ticks::Int64)
  for _ ∈ 1:ticks
    pos += vel
    vel *= 0.99f0
    vel -= 0.05f0
  end
  (pos=pos, vel=vel)
end

function tick_y_list(::Type{PersistentProjectileEntity}, pos::Float64, vel::Float64, ticks::Int64)
  posout = Float64[]
  velout = Float64[]
  for _ ∈ 1:ticks
    pos += vel
    vel *= 0.99f0
    vel -= 0.05f0
    push!(posout, pos);
    push!(velout, vel);
  end
  (pos=posout, vel=velout)
end

function tick_y_honey(::Type{PersistentProjectileEntity}, pos::Float64, vel::Float64, ticks::Int64)
  for _ ∈ 1:ticks
    pos += vel
    vel *= 0.99f0
    vel -= 0.05f0
    if vel < -0.08e0 vel = -0.05e0 end
  end
  (pos=pos, vel=vel)
end

function tick(::Type{TNTEntity}, pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int64=1; NoGravity::Bool=false)
  g = Float64[0, -0.04, 0]
  for _ ∈ 1:ticks
    if !NoGravity vel += g end
    if sum(vel.^2) > 1e-7
      pos += vel
    end
    vel *= 0.98e0
  end
  (pos=pos, vel=vel)
end

function tick_list(::Type{TNTEntity}, pos::Vector{Float64}, vel::Vector{Float64}, ticks::Int64=1; NoGravity::Bool=false)
  posout = Vector{Float64}[]
  velout = Vector{Float64}[]
  g = Float64[0, -0.04, 0]
  for _ ∈ 1:ticks
    if !NoGravity vel += g end
    if sum(vel.^2) > 1e-7
      pos += vel
    end
    vel *= 0.98e0
    push!(posout, pos);
    push!(velout, vel);
  end
  (pos=posout, vel=velout)
end

function tick_y(::Type{TNTEntity}, pos::Float64, vel::Float64, ticks::Int64)
  for _ ∈ 1:ticks
    vel -= 0.04e0
    pos += vel
    vel *= 0.98e0
  end
  (pos=pos, vel=vel)
end

function tick_y_list(::Type{TNTEntity}, pos::Float64, vel::Float64, ticks::Int64)
  posout = Float64[]
  velout = Float64[]
  for _ ∈ 1:ticks
    vel -= 0.04e0
    pos += vel
    vel *= 0.98e0
    push!(posout, pos);
    push!(velout, vel);
  end
  (pos=posout, vel=velout)
end

function tick_y_honey(::Type{TNTEntity}, pos::Float64, vel::Float64, ticks::Int=1)
  for _ ∈ 1:ticks
    vel -= 0.04e0
    pos += vel
    vel *= 0.98e0
    if vel < -0.08e0 vel = -0.05e0 end
  end
  (pos=pos, vel=vel)
end
