"A Julia library for Minecraft simulations useful for making TNT cannons"
module CannonLib

using Memoize

include("Entity.jl")
export blockdict, eyeheightdict, cannondict;
export cannon_angles, cannon_destinations, cannon_encode, get_slime_bounces, get_slime_bounces_with_honey, sim_bounces

const blockdict = Dict([
"all" => (0,1,1.5,2,3,4,5,6,7,8,9,9.5,10,11,12,13,14,15,16)./16 .+ 0.0625e0*0.98f0,
"movable" => (0,1,3,8,9,9.5,10,14,15,16)./16 .+ 0.0625e0*0.98f0,
"movable_0ticks" => (0,1,2,3,4,6,8,9,9.5,10,12,14,15,16)./16 .+ 0.0625e0*0.98f0,
])

const eyeheightdict = Dict([
"ender_pearl" => 0.2125f0,
"snowball" => 0.2125f0,
"player" => 0.4f0,
])

const cannondict = Dict([
"day6" => (efficiencyaxis=0.5948991084261209e0, basketparity=false, bits=[1764 588 588 252 168 84 42 21 14 7 4 2 1; ones(Int64, 13)'], sides=("red", "blue")),
"ghost" => (efficiencyaxis=0.0, basketparity=false, bits=[5219 448 28 7 1;4 4 4 2 3], sides=("red", "blue")),
])

"""
    cannon_angles(tnt::Int)

Number of possible angles given a maximum TNT count for 1 side of the cannon. This is basically just a fast sum of totients from 1 to `tnt`; see https://oeis.org/A002088.

For a differential basket, just use `cannon_angles(totaltnt) ÷ 2`
"""
function cannon_angles(tnt::Integer)
  if tnt < 1 throw(ArgumentError("ϕ(0) is not defined")) end
  if tnt == 1 return 8 end
  return 8 * (_cannon_angles(tnt) + 1)
end

# qwr (https://math.stackexchange.com/users/122489/qwr), How to calculate these totient summation sums efficiently?, URL (version: 2019-04-30): https://math.stackexchange.com/q/1740370
# Andy (https://math.stackexchange.com/users/91609/andy), How to calculate these totient summation sums efficiently?, URL (version: 2018-05-29): https://math.stackexchange.com/q/475219
# function _cannon_angles!(n::Integer, vec::Vector{Integer}, dict::Dict{Integer, Integer})
#   if n < 3 return 1 end
#   if n == 3 return 3 end
#   sum = muladd(n, n, 0-n) >> 1 # No idea why 0-n is faster than -n
#   a = 2
#   b = n >> 1
#   while true
#     sum -= b < length(vec) ? vec[b] :
#       haskey(dict, b) ? dict[b] :
#       dict[b] = _cannon_angles!(b, vec, dict)
#     if a == b return sum end
#
#     preva = a
#     sum -= (b - (b = n ÷ (a += 1))) * (
#       preva <= length(vec) ? vec[preva] :
#         (x = _cannon_angles!(preva, vec, dict); push!(vec, x); x)
#     )
#     if preva == b return sum end
#   end
# end

@memoize function _cannon_angles(n::Integer)
  if n < 3 return 1 end
  if n == 3 return 3 end
  sum = muladd(n, n, 0-n) >> 1 # No idea why 0-n is faster than -n
  a = 2
  b = n >> 1
  while true
    sum -= _cannon_angles(b)
    if a == b return sum end

    preva = a
    sum -= (b - (b = n ÷ (a += 1))) * _cannon_angles(preva)
    if preva == b return sum end
  end
end

"""
    cannon_destinations(efficiencyaxis, startx, startz, endx, endz, ticks)

Given an efficiency, start position, and target position, returns the tnt count per basket corner.

Returns (tnt_northwest::Vector{Int64}, tnt_southwest::Vector{Int64}, final_pos_x::Vector{Float64}, final_pos_z::Vector{Float64}, distance::Vector{Float64}).

# Arguments
- `efficiencyaxis::Float64`: The velocity given to the projectile in each horizontal axis by 1 tnt. Can be obtained by looking at the output of `/log explosions full` to see applied velocity.
- `startx::Float64`: The projectile starting x position.
- `startz::Float64`: The projectile starting z position.
- `endx::Float64`: The projectile ending x position.
- `endz::Float64`: The projectile ending z position.
- `ticks::Int64`: The maximum number of ticks for the projectile to fly.
"""
function cannon_destinations(efficiencyaxis::Float64, startx::Float64, startz::Float64, endx::Float64, endz::Float64, ticks::Int)
  tnt_northwest, tnt_southwest = Vector{Int}(undef, ticks), Vector{Int}(undef, ticks)
  pos_x, pos_z = Vector{Float64}(undef, ticks), Vector{Float64}(undef, ticks)
  dist = Vector{Float64}(undef, ticks)

  da = endx - startx + endz - startz
  db = endx - startx - (endz - startz)
  da, db = (da, db) ./ 2efficiencyaxis
  for t ∈ 1:ticks
    deltap_over_v = 100(1 - 0.99^t)
    tntnw, tntsw = round.(Int, (da, db) ./ deltap_over_v)
    finalx = muladd(tntnw + tntsw, deltap_over_v * efficiencyaxis, startx)
    finalz = muladd(tntnw - tntsw, deltap_over_v * efficiencyaxis, startz)

    tnt_northwest[t], tnt_southwest[t] = tntnw, tntsw
    pos_x[t], pos_z[t] = finalx, finalz
    dist[t] = sqrt((finalx - endx)^2 + (finalz - endz)^2)
  end
  println("ticks, tnt_northwest, tnt_southwest, pos_x, pos_z, distance:")
  display(Any[[1:ticks...] tnt_northwest tnt_southwest pos_x pos_z dist])

  return (tnt_northwest=tnt_northwest, tnt_southwest=tnt_southwest, final_pos_x=pos_x, final_pos_z=pos_z, distance=dist)
end

"""
    cannon_destinations(cannon, startx, startz, endx, endz, ticks)

Given a cannon, start position, and target position, returns the tnt count per basket corner.

Returns (tnt_northwest::Vector{Int64}, tnt_southwest::Vector{Int64}, final_pos_x::Vector{Float64}, final_pos_z::Vector{Float64}, distance::Vector{Float64}).

# Arguments
- `cannon::String`: The name of the cannon.
- `startx::Float64`: The projectile starting x position.
- `startz::Float64`: The projectile starting z position.
- `endx::Float64`: The projectile ending x position.
- `endz::Float64`: The projectile ending z position.
- `ticks::Int64`: The maximum number of ticks for the projectile to fly.
"""
cannon_destinations(cannon::String, startx::Float64, startz::Float64, endx::Float64, endz::Float64, ticks::Int) = cannon_destinations(cannondict[cannon][:efficiencyaxis], startx, startz, endx, endz, ticks)

function cannon_encode(tnt_northwest::Int, tnt_southwest::Int, parity::Bool, bits::Matrix{Int})
  nbits = size(bits)[2]
  dirindex = 2sign(tnt_northwest) + sign(tnt_southwest) # Get direction from signs of tnt
  tnt_northwest, tnt_southwest = abs.((tnt_northwest, tnt_southwest))
  swap = (dirindex & 0x2 == 0x2) ⊻ parity # Determine if cannon sides should swap
  dir = ("West", "North", "South", "East")[dirindex+1]

  bits_west, bits_east = Vector{Int}(undef, nbits), Vector{Int}(undef, nbits)
  for i ∈ 1:nbits
    # Calculate bits i up to the max bits[2,:]
    bit_west, bit_east = min.((tnt_northwest, tnt_southwest) .÷ bits[1,i], bits[2,i])
    tnt_northwest -= bit_west * bits[1,i]
    tnt_southwest -= bit_east * bits[1,i]

    bits_west[i], bits_east[i] = bit_west, bit_east
  end
  return swap ?
    (dir=dir, bits_west=bits_west, bits_east=bits_east) :
    (dir=dir, bits_west=bits_east, bits_east=bits_west)
end

cannon_encode(tnt_northwest::Int, tnt_southwest::Int, parity::Bool, bits::Vector{Int}) = cannon_encode(tnt_northwest, tnt_southwest, parity, [bits; ones(length(bits))'])

"""
    cannon_encode(cannon, tnt_northwest, tnt_southwest, rotation, mirror)

Given a cannon, tnt counts, and rotation/mirror, prints the bits to put in the ROM.

See also [`cannon_destinations`](@ref).

# Arguments
- `cannon::String`: The name of the cannon.
- `tnt_northwest::Int64`: The number of tnt in the northwest corner of the basket, or southeast if negative.
- `tnt_northeast::Int64`: The number of tnt in the northeast corner of the basket, or southwest if negative.
- `rotation::String`: The schematic rotation. Can be `"CW_90"`, `"CW_180"`, `"CCW_90"`, or `"NONE"`.
- `mirror::String`: The schematic rotation. Can be `"FRONT_BACK"`, `"LEFT_RIGHT"`, or `"NONE"`.
"""
function cannon_encode(cannon::String, tnt_northwest::Int, tnt_southwest::Int, rotation::String, mirror::String)
  parity = mirror == "FRONT_BACK" ? 0x3 : # Generate transformation from mirror
    mirror == "LEFT_RIGHT" ? 0x1 : 0x0
  case = rotation == "CW_90" ? 0x3 : # Generate transformation from rotation
    rotation == "CW_180" ? 0x2 :
    rotation == "CCW_90" ? 0x1 : 0x0

  c = cannondict[cannon]
  parity = c[:basketparity] ⊻ parity # Apply mirror transformation (litematica applies mirror before rotation)
  parity ⊻= Bool(parity & 1) && Bool(case & 1) ? case ⊻ 0x2 : case # Apply rotation transformation

  dir, bits_west, bits_east = cannon_encode(tnt_northwest, tnt_southwest, Bool(parity & 0x1), c[:bits])
  e, w = [c[:bits][1,:] bits_east], [c[:bits][1,:] bits_west]

  println("Direction: $dir")
  println("$(c[:sides][1]) side:")
  display(w')
  println("$(c[:sides][2]) side:")
  display(e')

  return parity & 0x2 == 0x2 ?
    (dir=dir, sides=c[:sides], side1=e, side2=w) :
    (dir=dir, sides=c[:sides], side1=w, side2=e)
end

"""
    cannon_encode(cannon, tnt_northwest, tnt_southwest)

Given a cannon, tnt counts, and rotation/mirror, prints the bits to put in the ROM.

# Arguments
- `cannon::String`: The name of the cannon.
- `tnt_northwest::Int64`: The number of tnt in the northwest corner of the basket, or southeast if negative.
- `tnt_northeast::Int64`: The number of tnt in the northeast corner of the basket, or southwest if negative.
"""
cannon_encode(cannon::String, tnt_northwest::Int64, tnt_southwest::Int64) = cannon_encode(cannon::String, tnt_northwest::Int64, tnt_southwest::Int64, "NONE", "NONE")

# function get_slime_bounces_with_honey(pos_::Float64, honeyticks::Int64, ticks::Tuple, threshold::Float64; entitytype::Type=ThrownEntity, explosionheight=Float64(0.98f0*0.0625f0), minpos=-256.0, maxpos=256.0, maxticks=999, prehoneyticks=0)
#   prevI = Tuple(0x0 for i ∈ ticks)
#   outdelta = Float64[]
#   outpos = Float64[]
#   outblock = Float64[]
#   outbounces = Tuple[]
#   outhoney = Int64[]
#   pos_ = tick_y(ThrownEntity, pos_, 0.0, prehoneyticks)[:pos]
#   starting_positions = tick_y_honey.(ThrownEntity, pos_, 0.0, 1:honeyticks)[:pos]
#   for i ∈ eachindex(starting_positions)
#     pos, bounces, delta, block = _get_slime_bounces((starting_positions[i], Float64.(prevI)...), prevI, CartesianIndices(ticks), threshold, entitytype, explosionheight, minpos, maxpos, maxticks)
#     outpos = vcat(outpos, pos)
#     for _ ∈ 1:length(pos) push!(outhoney, i) end
#     outbounces = vcat(outbounces, bounces)
#     outdelta = vcat(outdelta, delta)
#     outblock = vcat(outblock, block)
#   end
#   outpos, outhoney, outbounces, outdelta, outblock
# end

function _slime_bounce_pos(pos::Float64)
  pos -= floor(pos)
  return 0.75 <= pos ? 0.51 :
    0.5 <= pos ? 0. :
    0.25 <= pos ? 1.51 :
    1.
end

function _slime_bounce_short_pulse(pos::Float64)
  if 0.5 <= pos - floor(pos) < 0.75
    return true
  end
  return false
end

function _nextpostuple(entitytype::Type, index::Tuple, previndex::Tuple, positions::NTuple{N, Float64})::NTuple{N, Float64} where N
  ntuple(i -> i != 1 && index[N-i+1] != previndex[N-i+1] ? tick_y(entitytype, positions[i-1] + _slime_bounce_pos(positions[i-1]), 1e0, index[N-i+1])[:pos] : positions[i], Val(N))
end

@inline function _nearestblockpos(pos::Float64)
  block = round(16pos)/16
  return block, pos - block
end

@inline function _nearestblockpos(::Type{ThrownEntity}, pos::Float64, explosionheight::Float64)::Tuple{Float64, Float64}
  return _nearestblockpos(pos + 0.2125f0 - explosionheight)
end

@inline function _nearestblockpos(::Type{PersistentProjectileEntity}, pos::Float64, explosionheight::Float64)::Tuple{Float64, Float64}
  return _nearestblockpos(pos + 0.13f0 - explosionheight)
end

@inline function _nearestblockpos(::Type{TNTEntity}, pos::Float64, explosionheight::Float64)::Tuple{Float64, Float64}
  return _nearestblockpos(pos - explosionheight)
end

export get_slime_bounces_new

function get_slime_bounces_new(pos::Float64, ticks::Tuple, threshold::Float64; entitytype::Type=ThrownEntity, explosionheight=Float64(0.98f0*0.0625f0), minpos=pos, maxpos=256.0, maxticks=sum(ticks))
  prevI = ticks .* 0
  pvs = [(pos = pos + 1., vel = 1e0 * 0.99f0 - 0.03f0) for _ ∈ ticks]
  _get_slime_bounces_new(pvs, prevI, CartesianIndices(ticks), threshold, entitytype, explosionheight, minpos, maxpos, maxticks)
end

function _get_slime_bounces_new(pvs::Vector{NamedTuple{(:pos, :vel), Tuple{Float64, Float64}}}, prevI::NTuple{N, Int}, iter::CartesianIndices, threshold::Float64, entitytype::Type, explosionheight::Float64, minpos::Float64, maxpos::Float64, maxticks::Int) where N
  outpos = Float64[]
  outticks = Tuple[]
  outdelta = Float64[]
  outblock = Float64[]

  function npvbr!(pvs::Vector{NamedTuple{(:pos, :vel), Tuple{Float64, Float64}}}, ds::NTuple{N, Int}, entitytype::Type)
    pv = last(pvs)
    @inbounds(for i ∈ N:-1:1
      pvs[i] = ds[i] == 0 ? pv = pvs[i] :
               ds[i] < 0 ? (pos = _slime_bounce_pos(pv[1]) + pv[1], vel = 1.) :
               pv = tick_y(entitytype, pvs[i][1], pvs[i][2])
    end)
  end

  for index ∈ iter
    I = Tuple(index)
    npvbr!(pvs, I .- prevI, entitytype)
    prevI = I
    if sum(I) > maxticks continue end

    pos = first(first(pvs))
    if pos < minpos || pos > maxpos continue end
    block, delta = _nearestblockpos(entitytype, pos, explosionheight)
    if abs(delta) < threshold
      push!(outpos, pos)
      push!(outticks, (last(I) + 1, reverse(Base.front(I) .- 1)...))
      push!(outdelta, delta)
      push!(outblock, block)
    end
  end
  (outpos, outticks, outdelta, outblock)
end

struct BounceIndices{N}
  indices::NTuple{N, Int}
  pos::Float64
  vel::Float64
  entity_type::Type
end

struct BounceIndex{N}
  I::NTuple{N, Int}
  P::NTuple{N, Float64}
  V::NTuple{N, Float64}
end

struct BounceIndices2{N}
  indices::NTuple{N, Int}
  pos::Float64
  vels::Vector{Float64}
end

struct BounceIndex2{N}
  I::NTuple{N, Int}
  P::NTuple{N, Float64}
end

function BounceIndices2(pos::Float64, vel::Float64, ticks::NTuple{N, Int}) where N
  vels = Vector{Float64}(undef, maximum(ticks) - 1)
  for i ∈ 1:maximum(ticks) - 1
    vels[i] = vel
    vel = vel * 0.99f0 - 0.03f0
  end
  return BounceIndices2(ticks, pos, vels)
end

@inline function __inc(state::NTuple{N, Int}, pos::NTuple{N, Float64}, vel::NTuple{N, Float64}, indices::NTuple{N, Int}, entity_type::Type) where N
  I = first(state)
  ts, tp, tv = Base.tail(state), Base.tail(pos), Base.tail(vel)
  if I < first(indices)
    NP, NV = tick_y(entity_type, first(pos), first(vel))
    return true, (I + 1, ts...), (NP, tp...), (NV, tv...)
  end
  first_zero, I, P, V = __inc(ts, tp, tv, Base.tail(indices), entity_type)
  FP = first(P)
  NP = first_zero ? FP + _slime_bounce_pos(FP) : FP
  return false, (1, I...), (NP, P...), (1., V...)
end

@inline function Base.iterate(iter::BounceIndices{N}) where N
  BI = BounceIndex(ntuple(i -> 1, Val(N)), ntuple(i -> iter.pos, Val(N)), ntuple(i -> iter.vel, Val(N)))
  return BI, BI
end

@inline function Base.iterate(iter::BounceIndices{N}, state::BounceIndex{N}) where N
  state.I == iter.indices && return nothing
  _, I, P, V = __inc(state.I, state.P, state.V, iter.indices, iter.entity_type)
  next = BounceIndex(I, P, V)
  return next, next
end

@inline Base.eltype(::Type{BounceIndices{N}}) where N = BounceIndex{N}
@inline Base.length(iter::BounceIndices{N}) where N = prod(iter.indices)

@inline function __inc(state::NTuple{N, Int}, pos::NTuple{N, Float64}, vel::Vector{Float64}, indices::NTuple{N, Int}) where N
  I = first(state)
  ts, tp = Base.tail(state), Base.tail(pos)
  if I < first(indices)
    return true, (I + 1, ts...), (first(pos) + vel[I], tp...)
  end
  first_zero, I, P = __inc(ts, tp, vel, Base.tail(indices))
  FP = first(P)
  NP = first_zero ? FP + _slime_bounce_pos(FP) : FP
  return false, (1, I...), (NP, P...)
end

@inline function Base.iterate(iter::BounceIndices2{N}) where N
  BI = BounceIndex2(ntuple(i -> 1, Val(N)), ntuple(i -> iter.pos, Val(N)))
  return BI, BI
end

@inline function Base.iterate(iter::BounceIndices2{N}, state::BounceIndex2{N}) where N
  state.I == iter.indices && return nothing
  _, I, P = __inc(state.I, state.P, iter.vels, iter.indices)
  next = BounceIndex2(I, P)
  return next, next
end

"""
    get_slime_bounces(pos, vel, ticks, threshold; entity_type, explosion_height, min_pos, max_pos, max_ticks)

Simulates all possible slime bounces within the limits specified, and returns only those which are closer to the perfect alignment than `threshold`.

Returns (pos::Vector{Float64}, vel::Vector{Float64}, ticks::Vector{Tuple}, delta::Vector{Float64}, block::Vector{Float64}).

See also [`sim_bounces`](@ref).

# Arguments
- `pos::Float64`: The starting Y position.
- `vel::Float64`: The starting Y velocity.
- `ticks::Tuple`: A list of maximum bounce lengths, in ticks, written (bounce 1 length, bounce 2 length, ..., bounce N length).
- `threshold::Float64`: The biggest acceptable distance from the perfect alignment.
- `entity_type::Type`: The entity type to bounce. (see also [`Entity`](@ref))
- `explosion_height::Float64=0.061250001192092896`: The explosion height of the explosive.
- `min_pos::Float64=pos`: The minimum projectile arrival position.
- `max_pos::Float64=256.0`: The minimum projectile arrival position.
- `max_ticks::Int=sum(ticks)`: The maximum total ticks.

# Examples
```jldoctest
julia> pos, vel, ticks, delta, block = get_slime_bounces(0.0, 1.0, (100, 100), 1e-4);

julia> [pos vel ticks delta block][sortperm(delta, by=abs), :]
42×5 Matrix{Any}:
 -36.3387    (-36.3387, 8.55311)   (48, 99)   2.72403e-6  -36.1875
 -35.3387    (-35.3387, -44.8919)  (98, 49)   2.72403e-6  -35.1875
   4.91125   (4.91125, -9.99654)   (71, 33)  -3.43297e-6    5.0625
   3.59876   (3.59876, -5.95435)   (67, 49)   8.31251e-6    3.75
[...]

julia> sim_bounces(0.0, 1.0, (48, 99))
(2, 1.9600000102072954, 0.9204000199310483)
(3, 2.8804000301383437, 0.8811960291799086)
(4, 3.7615960593182525, 0.842384077962402)
[...]
```
"""
function get_slime_bounces(pos::Float64, vel::Float64, ticks::NTuple{N, Int}, threshold::Float64; entity_type=ThrownEntity, explosion_height=0.061250001192092896, min_pos=pos, max_pos=256., max_ticks=sum(ticks)) where N
  return _get_slime_bounces(BounceIndices(ticks, pos + vel, vel * 0.99f0 - 0.03f0, entity_type), threshold, entity_type, explosion_height, min_pos, max_pos, max_ticks)
end

function _get_slime_bounces(iter::BounceIndices{N}, threshold::Float64, entity_type::Type, explosion_height::Float64, min_pos::Float64, max_pos::Float64, max_ticks::Int) where N
  outpos = Float64[]
  outvel = Float64[]
  outticks = NTuple{N, Int}[]
  outdelta = Float64[]
  outblock = Float64[]
  for BI ∈ iter
    pos = first(BI.P)
    if pos < min_pos || pos > max_pos || sum(BI.I) > max_ticks
      continue
    end
    block, delta = _nearestblockpos(entity_type, pos, explosion_height)
    if abs(delta) < threshold
      push!(outpos, pos)
      push!(outvel, first(BI.V))
      push!(outticks, reverse(BI.I) .- 1)
      push!(outdelta, delta)
      push!(outblock, block)
    end
  end
  (outpos, outvel, outticks, outdelta, outblock)
end

"""
    sim_bounces(pos::Float64, indices::Tuple; entitytype::Type)

Simulates a set of bounces listed in `indices`, starting at `pos`;

Useful for reading the output of [`get_slime_bounces`](@ref).

See also [`get_slime_bounces`](@ref).

# Examples
```jldoctest
julia> sim_bounces(0.0, (45, 57))
Bounced: +0.51, long pulse
(0, 0.51, 1.0)
(1, 1.51, 0.9600000102072954)
(2, 2.470000010207295, 0.9204000199310483)
(3, 3.3904000301383435, 0.8811960291799086)
[...]
(45, 11.035809967136265, -0.45525796911456506)
Bounced: +1.0, long pulse
(45, 12.035809967136265, 1.0)
(46, 13.035809967136265, 0.9600000102072954)
[...]
(101, 16.1953340411928, -0.7215951635213262)
(102, 15.473738877671474, -0.7443792180972285)
```
"""
function sim_bounces(pos::Float64, vel::Float64, indices::Tuple; entitytype::Type=ThrownEntity)
  pos, vel, tick = pos + vel, vel * 0.99f0 - 0.03f0, 1
  pos, vel = tick_y(entitytype, pos, vel)
  for (i, bounce) ∈ enumerate(indices)
    for _ ∈ 2:bounce
      println((tick += 1, pos, vel))
      pos, vel = tick_y(entitytype, pos, vel)
    end
    if i != length(indices)
      println("Piston extension on tick $(tick + 1)")
      for d ∈ _slime_bounce(pos)
        println((tick += 1, pos += d, 1e0))
        pos, vel = tick_y(entitytype, pos, 1e0)
      end
      println("Piston retraction on tick $(tick + 1)")
    end
  end
  println((tick += 1, pos, vel))
end

function _slime_bounce(pos::Float64)
  pos -= floor(pos)
  return 0.75 <= pos ? (0.51,) :
    0.5 <= pos ? (0.,) :
    0.25 <= pos ? (0.51, 0.) :
     (0., 0.)
end

end
