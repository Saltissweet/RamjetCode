
# Helpful Equations
import RequiredEquations as he
import numpy as np
import math
import matplotlib.pyplot as plt
import cea


mach0 = 3.2
static_pressure_10k = 69680  #{Pascals @ 10k}
temp_0 = 268.3  #{Kelvin @ 10k}
tot_press_0 = he.tp(mach0,static_pressure_10k)
density_0 = .905 # {kg/m^3}



"""Inlet"""
initangle = math.ceil(math.degrees(math.asin(1/mach0)))
press_tot_ideal = 0
angles = [0,0]
for i in np.arange(initangle+3,40,.5):
    output1 = he.osr(mach0, i, tot_press_0, temp_0, density_0)
    initangle2 = math.ceil(math.degrees(math.asin(1/output1[0])))
    for j in np.arange(initangle2+3,50,.5):
        output2 = he.osr(output1[0], j, output1[2], output1[3], output1[4])
        output3 = he.nsr(output2[0], output2[2], output2[3], output2[4])
        press_tot_final = output3[1]
        if press_tot_final > press_tot_ideal:
            press_tot_ideal = press_tot_final
            angles = [i,j]

# Recalc Final Outputs so you can print
output1_final = he.osr(mach0, angles[0], tot_press_0, temp_0, density_0)
output2_final = he.osr(output1_final[0], angles[1], output1_final[2], output1_final[3], output1_final[4])
output3_final = he.nsr(output2_final[0], output2_final[2], output2_final[3], output2_final[4])
static_pressure_final = he.sp(output3_final[0], output3_final[1])


# Inlet Outputs
# print("angles: " + str(angles))
# print("After Shockwave 1: " + str(output1_final))
# print("After Shockwave 2: " + str(output2_final))
# print("After Shockwave 3: " + str(output3_final))
# print("Inlet Pressure Recovery: " + str(output3_final[1]/tot_press_0))

"""Subsonic Diffuser"""
inlet_sub_mach = output3_final[0]
sub_press_tot = output3_final[1]
inlet_sub_temp = output3_final[2]
inlet_sub_dens = output3_final[3]
inlet_sub_area = 0.0182414692  #{m^2}
exit_sub_mach = .15


sub_temp_tot = he.tt(inlet_sub_mach, inlet_sub_temp)
sub_dens_tot = he.td(inlet_sub_mach, inlet_sub_dens)
inlet_sub_velocity = he.speed(inlet_sub_mach, inlet_sub_temp)[0]
inlet_sub_massflow = he.massflow(inlet_sub_velocity, inlet_sub_dens, inlet_sub_area)

theor_throat_area_sub = he.fta(inlet_sub_mach, inlet_sub_area)
# e for exit
e_sub_area = he.fwa(exit_sub_mach, theor_throat_area_sub)
e_sub_press_stat = he.sp(exit_sub_mach, sub_press_tot)
e_sub_temp_stat = he.st(exit_sub_mach, sub_temp_tot)
e_sub_dens_stat = he.sd(exit_sub_mach, sub_dens_tot)

# print("Mass Flow Inlet: (kg/s?) " + str(inlet_sub_massflow))
# print("Subsonic Diffuser Exit Area: " + str(exit_sub_area))
# print("Subsonic Diffuser Exit Pressure S: " + str(exit_sub_press_static))
# print("Subsonic Diffuser Exit Temperature S: " + str(exit_sub_temp_static))
# print("Subsonic Diffuser Exit Density S: " + str(exit_sub_dens_static))
# print("Subsonic Diffuser Inlet Mach: " + str(inlet_sub_mach))
# print("Subsonic Diffuser Exit Mach: " + str(exit_sub_mach))


"""Combustor"""
# Part 1: Short distance (only Fanno Flow), I guess a little break and mixing area
# Part 2: Combustor (Rayliegh + Fanno Flow)
# Part 3: IDK another short distance (Fanno Only)
# IG I should change prev stuff to include Fanno, im too lazy for that rn tho

# Part 2 Method: Divide up the combustor into x steps (distance)
# At each distance, apply Fanno Flow, Rayliegh with totalheat/x heat added,
#       and then work the isentropic backwords to make it const P

# From now it is slop numbers lol
# Start off with x = 1e4 steps?
# Total Dist: Allot 1 meter (big), [Part 1: .2m, Part 2: .7m, Part3: .1m]
# Assume Stainless Steel For Friction Factor = .02
# Assume Heat Additions is 1000 kj/kg
heatadd = 1e6
frictfact = .02

# Get "Diameter", but I actually made a 2-D ramjet, so it is actually cuboidal shhhh
diam_comb1 = 2*math.sqrt(e_sub_area/math.pi)
dist1 = .2

# Part 1
exit_comb_pt1 = he.ff(exit_sub_mach, e_sub_press_stat, e_sub_temp_stat, e_sub_dens_stat, frictfact, dist1, diam_comb1)
e1_mach = exit_comb_pt1[0]
e1_press = exit_comb_pt1[1]
e1_temp = exit_comb_pt1[2]
e1_dens = exit_comb_pt1[3]

# Part2
numsteps = int(1e4)
locs = np.linspace(0, 0.7, numsteps)
d_dist = .7/numsteps
# if you want to in the future you could make a list of vals of heatadd which all add up to total heatadd
# this gives more control on how much heat is added in different places
# just make a function which models an equation, like a quadratic, and normalize it or something (add em all up and divide by it
# to get it to 1 and multiply by total heatadd), this could give you a varied profile of heat addition (which I believe a combustor
# prolly has) (I chose d_heatadd bc its like a small value (dx))
d_heatadd = heatadd/numsteps

# Initialize a list of these vals so you can plot em later
c_areas = [0] * numsteps
c_machs = [0] * numsteps
c_pressures = [0] * numsteps
c_temps = [0] * numsteps
c_densities = [0] * numsteps

# Initialize first value of list 
c_areas[0] = e_sub_area
c_machs[0] = e1_mach
c_pressures[0] = e1_press
c_temps[0] = e1_temp
c_densities[0] = e1_dens


# First run a Rayliegh with that amount of heat added
# Then run a Fanno with that amount of heat added
# Then work back isentropic to figure out how much area needs to change to keep static P
for i in range(1, numsteps):
    diam = 2*math.sqrt(c_areas[i-1]/math.pi)
    valsR = he.rf(c_machs[i-1], c_pressures[i-1], c_temps[i-1], c_densities[i-1], d_heatadd)
    valsF = he.ff(valsR[0], valsR[1], valsR[2], valsR[3], frictfact, d_dist, diam)
    
    # Post Fanno Totals
    pf_mach = valsF[0]
    pf_presstot = he.tp(valsF[0], valsF[1])
    pf_temptot = he.tt(valsF[0], valsF[2])
    pf_denstot = he.td(valsF[0], valsF[3])

    current_area = c_areas[i-1]
    if pf_presstot < e1_press:
        print(i)
        raise ValueError("why is total pressure so low gng")
    
    post_isen_mach = he.fmp(pf_presstot, e1_press)
    c_throat_area = he.fta(pf_mach, current_area)
    c_new_area = he.fwa(post_isen_mach, c_throat_area)
    c_new_temp = he.st(post_isen_mach, pf_temptot)
    c_new_dens = he.sd(post_isen_mach, pf_denstot)
    # Should be unnecessary, just for checking for now
    c_new_press = he.sp(post_isen_mach, pf_presstot)

    # Update lists
    c_areas[i] = c_new_area
    c_machs[i] = post_isen_mach
    c_pressures[i] = c_new_press
    c_temps[i] = c_new_temp
    c_densities[i] = c_new_dens

# Combustor Radius
result = [x / math.pi for x in c_areas]
c_rads = [math.sqrt(x) for x in result]


# Exit Part 2 vals
e2_area = c_areas[numsteps-1]
e2_mach = c_machs[numsteps-1]
e2_press = c_pressures[numsteps-1]
e2_temp = c_temps[numsteps-1]
e2_dens = c_densities[numsteps-1]

# Part 3
# Get "Diameter"
diam_comb3 = 2*math.sqrt(e2_area / math.pi)
dist3 = .1

exit_comb_pt1 = he.ff(e2_mach, e2_press, e2_temp, e2_dens, frictfact, dist3, diam_comb3)
e3_mach = exit_comb_pt1[0]
e3_press = exit_comb_pt1[1]
e3_temp = exit_comb_pt1[2]
e3_dens = exit_comb_pt1[3]


# Mach, Press, Temp, Density,
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows=5, ncols=1, sharex=True)

ax1.plot(locs, c_machs)
ax1.set_ylabel('Mach # Across Combustor')

ax2.plot(locs, c_pressures)
ax2.set_ylabel('Pressure Across Combustor')

ax3.plot(locs, c_temps)
ax3.set_ylabel('Temperature Across Combustor')

ax4.plot(locs, c_densities)
ax4.set_ylabel('Density Across Combustor')

ax5.plot(locs, c_rads)
ax5.set_ylabel('Radius Across Combustor')

ax5.set_xlabel('Location Across Combustor')

plt.tight_layout()
plt.grid()
plt.show()