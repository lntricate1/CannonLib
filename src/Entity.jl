"A Julia library for simulating Minecraft entities"
module Entity

export simprojectile, tickprojectile, tickprojectilelist, tickprojectile!,
       simtnt, ticktnt, ticktnt!,
       projectilecalcvel
const projectile_pos_x_table = Float64[1.0, 1.9900000095367432, 2.9701000284194947, 3.940399056460381, 4.900995103474351, 5.851985199179139, 6.793465402996228, 7.7255308137538, 8.648275579292665, 9.561792905976121, 10.466175068104723, 11.3615134172369, 12.247898391416367, 13.125419524307265, 13.994165454237947, 14.85422393315433, 15.705681835483706, 16.548625166909922, 17.38313907306081, 18.209307848108736, 19.02721494328514, 19.836942975309952, 20.638573734736685, 21.432188194214095, 22.21786651666523, 22.995688063384662, 23.765731402054787, 24.528074314681916, 25.282793805453043, 26.029966108514024, 26.769666695669986, 27.501970284008724, 28.226950843447863, 28.944681604206565, 29.655235064202493, 30.358682996374828, 31.05509645593404, 31.74454578753918, 32.427100632403366, 33.102829935328266, 33.771801951668174, 34.434084254224494, 35.08974374007127, 35.738846637312434, 36.38145851177151, 37.01764427361442, 37.64746818390604, 38.27099386110122, 38.888284287470846, 39.49940181546372, 40.10440817400473, 40.70336447473012, 41.296331218160354, 41.88336829981125, 42.464535016244064, 43.039890071054984, 43.609491580804814, 44.17339708088929, 44.73166353135074, 45.28434732263162, 45.83150428127049, 46.37318967554107, 46.90945822103486, 47.44036408618797, 47.96596089775266, 48.48630174621418, 49.001439191153445, 49.51142526655605, 50.01631148606824, 50.516148848200274, 51.01098784147781, 51.500878449541716, 51.98587015619695, 52.46601195041087, 52.94135233126165, 53.41193931283711, 53.877820429084686, 54.339042738612775, 54.795652829444144, 55.24769682372177, 55.695220382367644, 56.138268709694984, 56.576886557974284, 57.01111823195378, 57.44100759333463, 57.86659806520142, 58.28793263640829, 58.70505386592125, 59.11800388711706, 59.52682441203911, 59.931556735610755, 60.33224173980651, 60.728919897781545, 61.12163127795984, 61.510415548081546, 61.895311979209765, 62.27635944969736, 62.65359644911403, 63.02706108213415, 63.3967910723857]
const projectile_vel_x_table = Float64[0.9900000095367432, 0.9801000188827516, 0.9702990280408862, 0.96059604701397, 0.950990095704788, 0.9414802038170884, 0.9320654107575724, 0.9227447655388652, 0.9135173266834564, 0.904382162128602, 0.8953383491321764, 0.8863849741794665, 0.8775211328908976, 0.8687459299306823, 0.8600584789163823, 0.8514579023293753, 0.8429433314262169, 0.8345139061508888, 0.8261687750479247, 0.8179070951764049, 0.8097280320248108, 0.801630759426731, 0.7936144594774104, 0.7856783224511336, 0.7778215467194346, 0.7700433386701245, 0.7623429126271288, 0.7547194907711261, 0.7471723030609808, 0.7397005871559613, 0.7323035883387362, 0.7249805594391401, 0.7177307607587021, 0.710553459995929, 0.7034479321723356, 0.6964134595592145, 0.6894493316051387, 0.6825548448641885, 0.6757293029248969, 0.6689720163399047, 0.6622823025563199, 0.6556594858467729, 0.6491028972411613, 0.6426118744590773, 0.636185761842911, 0.6298239102916221, 0.6235256771951747, 0.6172904263696272, 0.6111175279928711, 0.6050063585410133, 0.5989563007253934, 0.5929667434302319, 0.5870370816509012, 0.581166716432814, 0.5753550548109236, 0.5696015097498278, 0.5639055000844728, 0.5582664504614501, 0.5526837912808793, 0.5471569586388739, 0.5416853942705806, 0.5362685454937893, 0.5309058651531068, 0.5255968115646886, 0.5203408484615235, 0.5151374449392653, 0.5099860754026061, 0.5048862195121863, 0.4998373621320346, 0.49483899327753483, 0.4898906080639119, 0.48499170665523367, 0.48014179421392267, 0.4753403808507724, 0.4705869815754638, 0.46588111624757633, 0.4612223095280891, 0.45661009083136694, 0.4520439942776264, 0.4475235586458776, 0.44304832732733607, 0.43861784827930084, 0.4342316739794936, 0.4298893613808546, 0.4255904718667905, 0.4213345712068696, 0.41712122951296055, 0.412950021195809, 0.4088205249220492, 0.404732323571645, 0.4006850041957568, 0.3966781579750292, 0.3927113801782966, 0.3887842701217012, 0.3848964311282199, 0.38104747048759613, 0.37723699941667205, 0.3734646330201177, 0.3697299902515528, 0.3660326938750572]
const projectile_pos_y_table = Float64[1.0, 1.9600000102072954, 2.8804000301383437, 3.7615960593182525, 4.603980137280654, 5.407940383167585, 6.173861034933361, 6.90212248815642, 7.593101334463043, 8.24717039956684, 8.86469878092784, 9.446051885034992, 9.991591464315839, 10.5016756536771, 10.976659006679844, 11.416892531352907, 11.822723725648187, 12.194496612541373, 12.532551774781682, 12.837226389294086, 13.108854261237521, 13.34776585772252, 13.55428834119166, 13.728745602466212, 13.871458293462325, 13.982743859580044, 14.062916571768438, 14.112287558270088, 14.131164836048113, 14.119853341898937, 14.07865496325393, 14.007868568673027, 13.907790038033413, 13.778712292416325, 13.620925323694978, 13.434716223826623, 13.220369213851676, 12.978165672602858, 12.708384165127248, 12.41130047082411, 12.087187611301344, 11.736315877953377, 11.358952859263267, 10.955363467831797, 10.525809967136265, 10.070551998021699, 9.589846604927153, 9.083948261849741, 8.553108898049034, 7.997577923494407, 7.417602254057922, 6.813426336455275, 6.185292172937336, 5.533439345734775, 4.858105041258239, 4.15952407405653, 3.437928910535204, 2.6935496924379754, 1.9266142600933178, 1.137348175428593, 0.32597474475403976, -0.5072849586810757, -1.3622120723578717, -2.2385899223805676, -3.136204001590275, -4.054841947897648, -4.994293522832209, -5.954350590306181, -6.934807095590678, -7.93545904450214, -8.956104482796896, -9.996543475771785, -11.056578088068774, -12.136012363681516, -13.234652306161866, -14.352305859024309, -15.48878288634635, -16.64389515356291, -17.81745630845276, -19.009281862315106, -20.219189171334413, -21.446997418131552, -22.69252759349946, -23.955602478321435, -25.23604662567026, -26.53368634308631, -27.848349675032903, -29.179866385527088, -30.528067940944112, -31.892787492993865, -33.27385986186755, -34.67112151955288, -36.08441057331612, -37.51356674934936, -38.9584313765812, -40.41884737064948, -41.894659218034136, -43.38571296034883, -44.891856178789624, -46.41293797873916]
const projectile_vel_y_table = Float64[0.9600000102072954, 0.9204000199310483, 0.8811960291799086, 0.842384077962402, 0.8039602458869308, 0.7659206517657762, 0.7282614532230592, 0.6909788463066233, 0.6540690651037971, 0.6175283813610001, 0.5813531041071519, 0.5455395792808478, 0.5100841893612624, 0.47498335300274397, 0.44023352467306304, 0.4058311942952787, 0.37177288689318605, 0.33805516224030896, 0.3046746145124034, 0.2716278719434352, 0.23891159648499832, 0.20652248346913912, 0.17445726127455186, 0.1427126909961127, 0.11128556611771812, 0.08017271218839506, 0.04937098650164992, 0.018877277778024094, -0.011311494149176142, -0.041198378645006944, -0.07078639458090297, -0.10007853063961335, -0.1290777456170882, -0.15778696872134637, -0.18620909986835446, -0.21434700997494702, -0.24220354124881768, -0.26978150747561025, -0.29708369430313886, -0.3241128595227661, -0.3508717333479673, -0.37736301869010896, -0.4035893914314698, -0.4295535006955312, -0.45525796911456506, -0.48070539309454546, -0.5058983430774115, -0.5308393638007077, -0.5555309745546271, -0.5799756694364848, -0.6041759176026468, -0.6281341635179386, -0.6518528272025611, -0.6753343044765362, -0.6985809672017084, -0.7215951635213262, -0.7443792180972285, -0.7669354323446574, -0.7892660846647248, -0.8113734306745533, -0.8332597034351155, -0.8549271136767959, -0.876377850022696, -0.8976140792097073, -0.9186379463073728, -0.939451574934561, -0.9600570674739716, -0.9804565052844972, -1.0006519489114618, -1.0206454382947556, -1.0404389929748892, -1.0600346122969875, -1.0794342756127433, -1.098639942480351, -1.1176535528624423, -1.1364770273220406, -1.1551122672165575, -1.1735611548898486, -1.1918255538623492, -1.2099073090193078, -1.2278082467971376, -1.245530175367906, -1.2630748848219762, -1.2804441473488248, -1.2976397174160512, -1.3146633319465952, -1.3315167104941836, -1.3482015554170224, -1.364719552049752, -1.3810723688736821, -1.3972616576853256, -1.4132890537632457, -1.4291561760332356, -1.4448646272318464, -1.4604159940682786, -1.4758118473846558, -1.4910537423146955, -1.5061432184407928, -1.5210817999495336, -1.5358709957856524]
export projectile_pos_x_table, projectile_vel_x_table, projectile_pos_y_table, projectile_vel_y_table

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
    for _ ∈ 1:ticks
      pos1 += vel1
      vel1 *= D
      vel1 += g
    end
  else
    for _ ∈ 1:(-ticks)
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
  posout::Vector{Vector{Float64}} = []
  velout::Vector{Vector{Float64}} = []
  for _ ∈ 1:(ticks[1]+length(ticks)-1)
    pos1 += vel1
    vel1 *= 0.99f0
    vel1 += [0f0, -0.03f0, 0f0]
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
    for _ ∈ 1:ticks
      pos += vel
      vel *= D
      vel += g
    end
  else
    for _ ∈ 1:(-ticks)
      vel -= g
      vel /= D
      pos -= vel
    end
  end
  (pos=pos, vel=vel)
end

"""
    tickprojectiley(pos, vel, ticks)

Accurate state of projectile Y position after `ticks` ticks.

Returns `pos::Float64`.
"""
function tickprojectiley(pos::Float64, vel::Float64, ticks::Int=1)::Float64
  pos1 = pos
  vel1 = vel
  for _ ∈ 1:ticks
    pos1 += vel1
    vel1 *= 0.99f0
    vel1 -= 0.03f0
  end
  pos1
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
    for _ ∈ 1:ticks
      vel1 += g
      pos1 += vel1
      vel1 *= D
    end
  else
    for _ ∈ 1:(-ticks)
      vel1 -= g
      pos1 -= vel1
      vel1 /= D
    end
  end
  (pos=pos1, vel=vel1)
end

function ticktnt(pos::Float64, vel::Float64, ticks::Int=1)::NamedTuple{(:pos, :vel), Tuple{Float64, Float64}}
  pos1 = pos
  vel1 = vel
  if ticks > 0
    for _ ∈ 1:ticks
      vel1 -= 0.04e0
      pos1 += vel1
      vel1 *= 0.98e0
    end
  else
    for _ ∈ 1:(-ticks)
      vel1 -= 0.04e0
      pos1 -= vel1
      vel1 /= 0.98e0
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
    for _ ∈ 1:ticks
      vel += g
      pos += vel
      vel *= D
    end
  else
    for _ ∈ 1:(-ticks)
      vel -= g
      pos -= vel
      vel /= D
    end
  end
  (pos=pos, vel=vel)
end


end
