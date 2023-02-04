"A Julia library for Minecraft simulations useful for making TNT cannons"
module CannonLib

include("Entity.jl")
import .Entity.tick_projectile_y
import .Entity.tick_projectile_y_honey
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

Number of possible angles given a maximum TNT count for 1 side of the cannon.
"""
function cannon_angles(tnt::Int)
  count::Int = 0
  for x::Int ∈ 2:(tnt-1)
    for z::Int ∈ x:tnt
      gcd(x, z) == 1 && (count += 1) end end
  count*8 + tnt*8
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
function cannon_destinations(efficiencyaxis::Float64, startx::Float64, startz::Float64, endx::Float64, endz::Float64, ticks::Int64)
  dx = endx - startx
  dz = endz - startz
  tnt_northwest = []
  tnt_southwest = []
  pos_x = []
  pos_z = []
  dist = []
  for t ∈ 1:ticks
    deltap_over_v = 100*(1 - 0.99^t)
    vx, vz = (dx, dz)./deltap_over_v
    tntnw, tntsw = round.(Int64, (vx + vz, vx - vz)./(2*efficiencyaxis))
    finalx, finalz = (tntnw + tntsw, tntnw - tntsw).*(deltap_over_v*efficiencyaxis) .+ (startx, startz)
    push!(tnt_northwest, tntnw)
    push!(tnt_southwest, tntsw)
    push!(pos_x, finalx)
    push!(pos_z, finalz)
    push!(dist, sqrt((finalx - endx)^2 + (finalz - endz)^2))
  end
  println("ticks, tnt_northwest, tnt_southwest, pos_x, pos_z, distance:")
  display([[1:ticks...] tnt_northwest tnt_southwest pos_x pos_z dist])
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
cannon_destinations(cannon::String, startx::Float64, startz::Float64, endx::Float64, endz::Float64, ticks::Int64) = cannon_destinations(cannondict[cannon][:efficiencyaxis], startx, startz, endx, endz, ticks)

function cannon_encode(tnt_northwest::Int64, tnt_southwest::Int64, parity::Bool, bits::Matrix{Int64})
  case = (tnt_northwest > 0) << 0x1 + (tnt_southwest > 0)
  swap = Bool(case & 0x2 == 0x2) ⊻ parity
  dir = ('w', 'n', 's', 'e')[case+1]
  bits_west = []
  bits_east = []
  tnt_northwest = abs(tnt_northwest)
  tnt_southwest = abs(tnt_southwest)
  for i ∈ eachindex(bits[1,:])
    bit_west = min(tnt_northwest ÷ bits[1,i], bits[2,i])
    bit_east = min(tnt_southwest ÷ bits[1,i], bits[2,i])
    tnt_northwest = tnt_northwest - bit_west * bits[1,i]
    tnt_southwest = tnt_southwest - bit_east * bits[1,i]
    push!(bits_west, bit_west)
    push!(bits_east, bit_east)
  end
  swap ? (dir=dir, bits_west=bits_west, bits_east=bits_east) : (dir=dir, bits_west=bits_east, bits_east=bits_west)
end

cannon_encode(tnt_northwest::Int64, tnt_southwest::Int64, parity::Bool, bits::Vector{Int64}) = cannon_encode(tnt_northwest, tnt_southwest, parity, [bits; ones(length(bits))'])

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
function cannon_encode(cannon::String, tnt_northwest::Int64, tnt_southwest::Int64, rotation::String, mirror::String)
  parity = mirror == "FRONT_BACK" ? 0x3 : # Generate transformation from mirror
    mirror == "LEFT_RIGHT" ? 0x1 : 0x0
  case = rotation == "CW_90" ? 0x3 : # Generate transformation from rotation
    rotation == "CW_180" ? 0x2 :
    rotation == "CCW_90" ? 0x1 : 0x0

  c = cannondict[cannon]
  parity = c[:basketparity] ⊻ parity # Apply mirror transformation (litematica applies mirror before rotation)
  parity ⊻= Bool(parity & 1) && Bool(case & 1) ? case ⊻ 0x2 : case # Apply rotation transformation

  dir, bits_west, bits_east = cannon_encode(tnt_northwest, tnt_southwest, Bool(parity & 0x1), c[:bits])
  e = [c[:bits][1,:] bits_east]'
  w = [c[:bits][1,:] bits_west]'
  println("Direction: $dir")
  println("$(c[:sides][1]) side:")
  display([c[:bits][1,:] bits_west]')
  println("$(c[:sides][2]) side:")
  display([c[:bits][1,:] bits_east]')

  parity & 0x2 == 0x2 ? (dir, c[:sides], e, w) : (dir, c[:sides], w, e)
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

function get_slime_bounces_with_honey(pos_::Float64, honeyticks::Int64, ticks::Tuple, threshold::Float64; eyeheight::Float32=eyeheightdict["ender_pearl"], explosionheight=Float64(0.98f0*0.0625f0), minpos=-256.0, maxpos=256.0, maxticks=999, prehoneyticks=0)
  prevI = Tuple(0x0 for i ∈ ticks)
  outdelta = Float64[]
  outpos = Float64[]
  outblock = Float64[]
  outbounces = Tuple[]
  outhoney = Int64[]
  pos_ = tick_projectile_y(pos_, 0.0, prehoneyticks)
  starting_positions = tick_projectile_y_honey.(pos_, 0.0, 1:honeyticks)
  for i ∈ eachindex(starting_positions)
    pos, bounces, delta, block = _get_slime_bounces((starting_positions[i], Float64.(prevI)...), prevI, CartesianIndices(ticks), threshold, eyeheight, explosionheight, minpos, maxpos, maxticks)
    outpos = vcat(outpos, pos)
    for _ ∈ 1:length(pos) push!(outhoney, i) end
    outbounces = vcat(outbounces, bounces)
    outdelta = vcat(outdelta, delta)
    outblock = vcat(outblock, block)
  end
  outpos, outhoney, outbounces, outdelta, outblock
end

function _slime_bounce_pos(pos::Float64)
  pos = pos - floor(pos)
  delta = 0.51
  if     0.5  < pos < 0.75 delta = 0.
  elseif 0.25 < pos < 0.5  delta = 1.51
  elseif 0.   < pos < 0.25 delta = 1. end
  return delta
end

function _slime_bounce_pulse(pos::Float64)
  pos = pos - floor(pos)
  pulse = 2
  if 0.5 < pos < 0.75 pulse = 1 end
  return pulse
end

function _nextpostuple(index::Tuple, previndex::Tuple, positions::NTuple{N1, Float64}) where N1
  ntuple(i -> i != 1 && index[N1-i+1] != previndex[N1-i+1] ? tick_projectile_y(positions[i-1] + _slime_bounce_pos(positions[i-1]), 1e0, index[N1-i+1]) : positions[i], Val(N1))
end

"""
    get_slime_bounces(pos, ticks, threshold; eyeheight, explosionheight, minpos, maxpos)

Simulates all possible slime bounces within the limits specified, and returns only those which are closer to the perfect alignment (where 0 vertical velocity added) than `threshold`.

Returns (position::Vector{Float64}, ticks::Vector{Tuple}, delta::Vector{Float64}, block::Vector{Float64}). A negative tick number means that bounce requires a 1gt pulse.

See also [`sim_bounces`](@ref).

# Arguments
- `pos::Float64`: The starting Y position.
- `ticks::Tuple`: A list of maximum bounce lengths, in ticks, in the form (bounce 1 length, bounce 2 length, ..., bounce N length).
- `threshold::Float64`: The biggest acceptable distance from the perfect alignment.
- `eyeheight::Float32`: The eye height of the projectile being aligned (optional, default 0.2125f0).
- `explosionheight::Float64`: The explosion height of the explosive (optional, default 0.061250001192092896).
- `minpos::Float64`: The minimum acceptable projectile arrival position (optional, default -256.0).
- `maxpos::Float64`: The minimum acceptable projectile arrival position (optional, default 256.0).

# Examples
```jldoctest
julia> pos, ticks, delta, block = get_slime_bounces(0.0, (100, 100), 1e-4);

julia> [pos ticks delta block][sortperm(delta, by=abs), :]
34×4 Matrix{Any}:
  15.4737    (45, 57)   -1.11176e-5   15.625
 -34.4013    (100, 14)  -1.23203e-5  -34.25
 -34.4013    (14, 100)  -1.23203e-5  -34.25
[...]

julia> sim_bounces(0.0, (45, 57))
Bounced: +0.51, long pulse
(0, 0.51, 1.0)
(1, 1.51, 0.9600000102072954)
[...]
```
"""
function get_slime_bounces(pos::Float64, ticks::Tuple, threshold::Float64; eyeheight::Float32=eyeheightdict["ender_pearl"], explosionheight=Float64(0.98f0*0.0625f0), minpos=-256.0, maxpos=256.0, maxticks=999)
  prevI = Tuple(0x0 for i ∈ ticks)
  poses = (pos, Float64.(prevI)...)
  _get_slime_bounces(poses, prevI, CartesianIndices(ticks), threshold, eyeheight, explosionheight, minpos, maxpos, maxticks)
end

function _get_slime_bounces(positions::NTuple{N, Float64}, previndex::NTuple{N1, UInt8}, iter::CartesianIndices, threshold::Float64, eyeheight::Float32, explosionheight::Float64, minpos::Float64, maxpos::Float64, maxticks::Int64) where {N, N1}
  outpos = Float64[]
  outticks = Tuple[]
  outdelta = Float64[]
  outblock = Float64[]
  yoffset = Float64(eyeheight) - explosionheight
  changeindexsign = (index, positions) -> ntuple(i -> _slime_bounce_pulse(positions[i]) == 1 ? -index[end-i+1] : index[end-i+1], length(index))
  for index ∈ iter
    index = Tuple(index)
    if sum(index) > maxticks continue end
    positions = _nextpostuple(index, previndex, positions) #Calculate next pos
    previndex = index
    if(minpos < last(positions) < maxpos)
      block = round(last(positions) + yoffset; base = 16, digits = 1);
      delta = last(positions) + yoffset - block;
      if abs(delta) < threshold
        push!(outpos, last(positions))
        push!(outticks, changeindexsign(index, positions))
        push!(outdelta, delta)
        push!(outblock, block)
      end
    end
  end
  (outpos, outticks, outdelta, outblock)
end

"""
    sim_bounces(pos::Float64, indices::Tuple)

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
function sim_bounces(pos::Float64, indices::Tuple)
  vel = 1e0
  tick = 0
  for bounce ∈ indices
    slimepos, slimepulse = _slime_bounce_pos(pos), _slime_bounce_pulse(pos)
    println("Bounced: +$slimepos, $(slimepulse == 1 ? "1gt" : "long") pulse")
    pos += slimepos
    vel = 1e0
    println((tick, pos, vel))
    for _ ∈ 1:bounce
      pos += vel
      vel *= 0.99f0
      vel -= 0.03f0
      tick += 1
      println((tick, pos, vel))
    end
  end
end

end
