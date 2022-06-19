"A Julia library for simulating Minecraft explosions"
module Explosion

import LinearAlgebra:cross,dot,rank,I
import RowEchelon:rref

export explosion, explosioncount, ec1

"""
    explosion(explosionpos::Vector{Float64}, entitypos::Vector{Float64}, eyeheight::Float64, power::Float32, exposure::Rational)::Vector{Float64}

Velocity vector applied on entity at `entitypos` with eye height `eyeheight` by an explosion at `explosionpos` with power `power` and exposure `exposure`.

# Arguments
- `explosionpos::Vector{Float64}`: Position of the explosion.
- `entitypos::Vector{Float64}`: Position of the affected entity.
- `eyeheight::Float32: Eye height of the affected entity.
- `power::Float32`: Explosion power.
- `exposure::Rational`: Exposure in the form (hit rays//total rays).
"""
function explosion(explosionpos::Vector{Float64}, entitypos::Vector{Float64}, eyeheight::Float64, power::Float32, exposure::Rational)::Vector{Float64}
  power < 0 && throw(DomainError(power, "explosion power must be nonnegative"))
  exposure < 0//1 && throw(DomainError(exposure, "exposure cannot be negative"))
  exposure > 1//1 && throw(DomainError(exposure, "exposure cannot be greater than 1"))

  dis = Float64(Float32(sqrt(sum((explosionpos - entitypos).^2)))/(2f0 * power))
  dis > 1e0 && return [0,0,0]

  dir = (entitypos + [0, eyeheight, 0] - explosionpos)
  dirmag = Float64(Float32(sqrt(sum(dir.^2))))
  dirmag == 0e0 && return [0,0,0]

  dir/dirmag * (1e0 - dis) * Float64(Float32(exposure))
end

"""
    explosion(explosionpos::Vector{Float64}, entitypos::Vector{Float64}, eyeheight::Float64, power::Float32, exposure::Rational, count::Int)::Vector{Float64}

Velocity vector applied on entity at `entitypos` with eye height `eyeheight` by an explosion at `explosionpos` with power `power` and exposure `exposure`.

# Arguments
- `explosionpos::Vector{Float64}`: Position of the explosion.
- `entitypos::Vector{Float64}`: Position of the affected entity.
- `eyeheight::Float32: Eye height of the affected entity.
- `power::Float32`: Explosion power.
- `exposure::Rational`: Exposure in the form (hit rays//total rays).
- `count::Int`: Number of explosions.
"""
function explosion(explosionpos::Vector{Float64}, entitypos::Vector{Float64}, eyeheight::Float64, power::Float32, exposure::Rational, count::Int)::Vector{Float64}
  count < 0 && throw(DomainError(count, "explosion count must be nonnegative"))
  vel1 = explosion(explosionpos, entitypos, eyeheight, power, exposure)
  vel = Float64[0,0,0]
  for i ∈ count vel += vel1 end
  vel
end

"""
    explosion(explosionpos::Vector{Float64}, entitypos::Vector{Float64}; type::String="tnt", type1::String="ender_pearl", exposure::Rational=8//8, count::Int = 1)::Vector{Float64}

Velocity vector applied on entity `type1` at `entitypos` by an explosion at `explosionpos` caused by `type`, with an exposure of `exposure`.

# Arguments
- `explosionpos::Vector{Float64}`: Position of the explosion.
- `entitypos::Vector{Float64}`: Position of the affected entity.
- `type::String`: The source of the explosion.
- `type1::String`: The type of entity affected.
- `exposure::Rational`: Exposure in the form (hit rays//total rays).
- `count::Int`: Number of explosions.
"""
function explosion(explosionpos::Vector{Float64}, entitypos::Vector{Float64}; type::String="tnt", type1::String="ender_pearl", exposure::Rational=8//8, count::Int = 1)::Vector{Float64}
  power = 4f0
  dy = Float64(0.98f0*0.0625e0)
  # THE explosives
  # tnt is default
  if type != "tnt"
    if type == "fireball" || type == "wither_skull"
      power = 1f0
      dy = 0e0
    elseif type == "creeper"
      power = 3f0
      dy = 0e0
    elseif type == "bed" || type == "respawn_anchor"
      power = 5f0
      dy = 0f0
    elseif type == "end_crystal" || type == "charged_creeper"
      power = 6f0
      dy = 0e0
    elseif type == "wither"
      power = 7f0
      dy = 0e0
    else throw(ArgumentError("Invalid explosion source"))
    end
  end

  eyeheight = 0.25f0*0.85f0
  # THE targets
  # ender_pearl, snowball, item and fishing_bobber are the same and default
  if type1 != "ender_pearl" && type1 != "snowball" && type1 != "item" && type1 != "fishing_bobber"
    if type1 == "tnt" eyeheight = 0f0
    elseif type1 == "arrow" eyeheight = 0.13f0
    elseif type1 == "player" eyeheight = 0.4f0
    elseif type1 == "falling_block" eyeheight = 0.98f0*0.85f0
      else throw(ArgumentError("Invalid entity type")) end
  end

  count == 1 ? explosion(explosionpos + [0, dy, 0], entitypos, Float64(eyeheight), power, exposure) :
    explosion(explosionpos + [0, dy, 0], entitypos, Float64(eyeheight), power, exposure, count)
end

"""
    explosioncount(explosionpos::Vector{Float64}, entitypos::Vector{Float64}, targetvel::Vector{Float64}, eyeheight::Float64, power::Float32, exposure::Rational)::Float64

Explosion count required to approach at least velocity `targetvel`.

# Arguments
- `explosionpos::Vector{Float64}`: Position of the explosion.
- `entitypos::Vector{Float64}`: Position of the affected entity.
- `targetvel::Vector{Float64}`: Target entity velocity.
- `eyeheight::Float32: Eye height of the affected entity.
- `power::Float32`: Explosion power.
- `exposure::Rational`: Exposure in the form (hit rays//total rays).
"""
function explosioncount(entitypos::Vector{Float64}, targetvel::Vector{Float64}, eyeheight::Float64, explosionpos::Vector{Float64}, power::Float32, exposure::Rational)::Float64
  dot(explosion(explosionpos, entitypos, eyeheight, power, exposure), targetvel)
end

"""
explosioncount(explosionpos::Vector{Float64}, entitypos::Vector{Float64}, targetvel::Vector{Float64}, explosions::Vector{Tuple{Vector{Float64}, Float32, Rational}})::Float64

Explosion count required to approach velocity `targetvel`.

# Arguments
- `entitypos::Vector{Float64}`: Position of the affected entity.
- `targetvel::Vector{Float64}`: Target entity velocity.
- `eyeheight::Float32: Eye height of the affected entity.
- `explosions::Vector{Tuple{Vector{Float64}, Float32, Rational}}`: An array of tuples `(explosionpos::Vector{Float64}, power::Float32, exposure::Rational)`.
"""
function explosioncount(entitypos::Vector{Float64}, targetvel::Vector{Float64}, eyeheight::Float64, explosions::Vector{Tuple{Vector{Float64}, Float32, Rational{Int}}})
  length(explosions) == 1 && return dot(explosion(entitypos, targetvel, eyeheight, explosions[1][1], explosions[1][2], explosions[1][3]), targetvel)
  vels = [Explosion.explosion(a[1], entitypos, eyeheight, a[2], a[3]) for a ∈ explosions]
  planes = []
  # Find all planes formed by 2 velocity vectors. If more than 2 vectors share a plane, keep only those with highest dot product with targetvel.
  for i ∈ 1:length(vels)
    v = vels[i]
    d = dot(targetvel, v)
    for j ∈ 1:length(vels)
      i == j && continue
      v1 = vels[j]
      c = cross(v, v1)
      c == Float64[0,0,0] && continue
      c /= sqrt(sum(c.^2))
      d1 = dot(targetvel, v1)
      plane = findfirst((x) -> c == x[1][1] || c == -x[1][1], planes) # Check if plane exists
      if plane === nothing
        c1 = abs(dot(targetvel, c)) # Add plane: ((plane cross product, dot with targetvel), (dot, vel, pos), (dot, vel, pos))
        push!(planes, d > d1 ? ((c, c1), (d, v, explosions[i][1]), (d1, v1, explosions[j][1])) : ((c, c1), (d1, v1, explosions[j][1]), (d, v, explosions[i][1])))
      else ## Replace only if more efficient
        if planes[plane][2][1] < d1
          planes[plane] = cross(planes[plane][2][2], v1) != [0e0,0e0,0e0] ?
            (planes[plane][1], (d1, v1, explosions[j][1]), planes[plane][2]) :
            (planes[plane][1], (d1, v1, explosions[j][1]), planes[plane][3])
        continue end
        if planes[plane][3][1] < d1 && cross(planes[plane][2][2], v1) != [0e0,0e0,0e0]
          planes[plane] = (planes[plane][1], planes[plane][2], (d1, v1, explosions[j][1]))
        continue end
      end
    end
  end

  # If only 1 plane exists, aim for the point in the plane closest to the target
  if length(planes) == 1
    plane = planes[1]
    v = (plane[2][2], plane[3][2])
    P = targetvel - plane[1][1]*plane[1][2]
    CX = cross(plane[1][1], v[1])
    T2 = dot(P, CX)/dot(v[2], CX)
    T1 = sqrt(sum((P - T2*v[2]).^2) / sum(v[1].^2))

    T1 == 0 && return (plane[3][3], T2)
    T2 == 0 && return (plane[2][3], T1)
    return [(plane[2][3], T1), (plane[3][3], T2)]
  end

  # If more than 1 plane exists, use the second most efficient plane to reach the most efficient plane
  sort!(planes, by = x -> x[1][2])
  #println()
  #for p in planes println(p) end
  #println()
  tnt3 = abs(dot(planes[2][2][2], planes[1][1][1])) > 0 ? planes[2][2] : planes[2][3] # Pick the vel that does not share plane 1
  v = tnt3[2]
  plane = planes[1]

  T3 = dot(targetvel, plane[1][1])/dot(v, plane[1][1])

  P = targetvel - v * T3
  #println("v: $v, P: $P, T3: $T3")

  v = (plane[2][2], plane[3][2])
  #println("v: $v")
  CX = cross(plane[1][1], v[1])
  T2 = dot(P, CX)/dot(v[2], CX)
  T1 = sqrt(sum((P - T2*v[2]).^2) / sum(v[1].^2))# * (plane[2][1] < 0 ? -1 : 1)

  tnt = [(plane[2][3], T1), (plane[3][3], T2), (tnt3[3], T3)]
  T2 == 0 && deleteat!(tnt, 2)
  T1 == 0 && popfirst!(tnt)
  T3 == 0 && pop!(tnt)
  tnt
end

"""
    explosioncount(explosionpos::Vector{Float64}, entitypos::Vector{Float64}, targetvel::Vector{Float64}, explosions::Tuple{Vector{Float64}, Float32, Rational}...)::Float64

Explosion count required to approach velocity `targetvel`.

# Arguments
- `entitypos::Vector{Float64}`: Position of the affected entity.
- `targetvel::Vector{Float64}`: Target entity velocity.
- `eyeheight::Float32: Eye height of the affected entity.
- `explosions::Tuple{Vector{Float64}, Float32, Rational}`: Any number of tuples `(explosionpos::Vector{Float64}, power::Float32, exposure::Rational)`.
"""
function explosioncount(entitypos::Vector{Float64}, targetvel::Vector{Float64}, eyeheight::Float64, explosions::Tuple{Vector{Float64}, Float32, Rational{Int}}...)
  explosioncount(entitypos, targetvel, eyeheight, explosions)
end

function ns(a)
  s = size(a)[2]
  r = rank(a)
  r == s && return 0
  b = [-rref(a)[1:r, :][:, r+1:s]; I(s-r)]
end

function ec1(entitypos::Vector{Float64}, targetvel::Vector{Float64}, eyeheight::Float64, explosions::Vector{Tuple{Vector{Float64}, Float32, Rational{Int}}})
  vels::Vector{Tuple{Vector{Float64}, Float64}} = []
  for e ∈ explosions
    v = Explosion.explosion(e[1], entitypos, eyeheight, e[2], e[3])
    d = dot(v, targetvel)
    d > 0 && push!(vels, (v, d))
  end
  sort!(vels, by = x -> -x[2])
  print(vels)
end

end
