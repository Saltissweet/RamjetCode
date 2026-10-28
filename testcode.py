import RequiredEquations as he
import numpy as np # type: ignore
import math
import matplotlib.pyplot as plt # type: ignore
import cea # type: ignore
import CoolProp.CoolProp as cp # type: ignore
from scipy.optimize import brentq # type: ignore

mach0 = 3.2
static_pressure_10k = 11597  #{Pascals @ 50k}
temp_0 = 268.3  #{Kelvin @ 10k}
tot_press_0 = he.tp(mach0,static_pressure_10k)
density_0 = .905 # {kg/m^3}
specgasconst = 287

## To Help Sai
# Virtual Environments:
# python -m venv env
# env\Scripts\activate
# deactivate
# pip freeze > requirements.txt (makes requirements file)
# pip install -r requirements.txt


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

# # gamma = k
"""Mach Number Diff eq"""
"""




REDO WITH NEG NUMBERS




"""
# Diff eq
# dM^2/M^2 = -M_Fa*dA + M_Fb*dQ + M_Fc*dx - M_Fd*dk

def M_Fa(Mach, k, Area):
    a = 2*(1+ .5*(k-1)*Mach**2)
    b = (1-Mach**2)
    return (a/(b*Area))

def M_Fb(Mach, k, C_p, Temp):
    a = (1 + k*Mach**2)/(1 - Mach**2)
    return (a/(C_p * Temp))

def M_Fc(Mach, k, Darcy, diameter):
    a = k*Mach**2 * (1 + (k-16/2)*Mach**2)
    b = (1-Mach**2)
    return (a*Darcy/(b*diameter))

def M_Fd(k):
    return 1/k

"""Press Diff Eq"""
# dP/P = P_Fa*dA - P_Fb*dQ - P_Fc*dx

def P_Fa(Mach, k, Area):
    a = k*Mach**2 / (1-Mach**2)
    return (a/Area)

def P_Fb(Mach, k, C_p, Temp):
    a = 1*k*Mach**2 / (1-Mach**2)
    return (a/(C_p*Temp))

def P_Fc(Mach, k, Darcy, diameter):
    a = 1*k*Mach**2*(1 + (k-1)*Mach**2)
    b = 2*(1-Mach**2)
    return (a*Darcy/(b*diameter))


"""Temp Diff Eq"""
# dP/P = T_Fa*dA + T_Fb*dQ - T_Fc*dx

def T_Fa(Mach, k, Area):
    a = (k-1)*Mach**2 / (1-Mach**2)
    return (a/Area)

def T_Fb(Mach, k, C_p, Temp):
    a = (1-k*Mach**2) / (1-Mach**2)
    return (a/(C_p * Temp))

def T_Fc(Mach, k, Darcy, diameter):
    a = 1*k*(k-1)*Mach**4
    b = 2*(1-Mach**2)
    return(a*Darcy/(b*diameter))

numsteps = 1e4
tot_dist = .7
tot_heat = 43.5e6
locs = np.linspace(0,tot_dist,numsteps)
machs = [0]*numsteps
machsquared = [0]*numsteps
pressures = [0]*numsteps
temps = [0]*numsteps
densities = [0]*numsteps
areas = [0]*numsteps
radiuses = [0]*numsteps

machs[0] = exit_sub_mach
machsquared[0] = exit_sub_mach**2
pressures[0] = e_sub_press_stat
temps[0] = e_sub_temp_stat
densities[0] = e_sub_dens_stat
areas[0] = e_sub_area
radiuses[0] = math.sqrt(e_sub_area/math.pi)

"""CHANGE numsteps and totvals are up^"""
darcy = .02
dx = tot_dist/numsteps
dQ = tot_heat/numsteps 
dP = 0
past_cp = 0
past_cv = 0
for i in range(1, numsteps):
    c_p = cp.PropsSI('CPMASS', 'P', pressures[i-1], 'D', densities[i-1], 'Air')
    c_v = cp.PropsSI('CVMASS', 'P', pressures[i-1], 'D', densities[i-1], 'Air')
    k = c_p/c_v

    mfa = M_Fa(machs[i-1], k, areas[i-1])
    mfb = M_Fb(machs[i-1], k, c_p, temps[i-1])
    mfc = M_Fc(machs[i-1], k, darcy, 2*radiuses[i-1])
    mfd = M_Fd(k)

    pfa = P_Fa(machs[i-1], k, areas[i-1])
    pfb = P_Fb(machs[i-1], k, c_p, temps[i-1])
    pfc = P_Fc(machs[i-1], k, darcy, 2*radiuses[i-1])

    tfa = T_Fa(machs[i-1], k, areas[i-1])
    tfb = T_Fb(machs[i-1], k, c_p, temps[i-1])
    tfc = T_Fc(machs[i-1], k, darcy, 2*radiuses[i-1])

    """What is dk
    Prolly not right"""
    dk = c_p/c_v - past_cp/past_cv

    # Start from dP/p = 0
    dA = (dP + pfb*dQ + pfc*dx)/pfa
    dT = (tfa*dA + tfb*dQ + tfc*dx)*temps[i-1]
    dmachsquared = (mfa*dA + mfb*dQ + mfc*dx + mfd*dk)

    areas[i] = areas[i-1] + dA
    temps[i] = temps[i-1] + dT
    machsquared[i] = machsquared[i-1] + dmachsquared
    if dmachsquared >= 0:
        machs[i] = machs[i-1] + math.sqrt(dmachsquared)
    else:
        machs[i] = machs[i-1] - math.sqrt(abs(dmachsquared))
    # Unnecessary
    pressures[i] = pressures[i-1] + dP
    # Ideal gas law is used for influence coeffients
    densities[i] = pressures[i]/(specgasconst*temps[i])
    if dA >= 0:
        radiuses[i] = radiuses[i-1] + math.sqrt(dA/math.pi)
    else:
        radiuses[i] = radiuses[i-1] - math.sqrt(dA/math.pi)


# Plotting Combustor Rayleigh + Fanno
fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True)
ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()

ax1.plot(locs, machs)
ax1.set_ylabel('Mach #')

ax2.plot(locs, pressures)
ax2.set_ylabel('Pressure')

ax3.plot(locs, temps)
ax3.set_ylabel('Temperature')

ax4.plot(locs, densities)
ax4.set_ylabel('Density')

ax5.plot(locs, radiuses)
ax5.set_ylabel('Radius')
ax5.set_xlabel('Location')

ax6.plot(locs, machsquared)
ax6.set_ylabel('Speed')
ax6.set_xlabel('Machsquared')

plt.tight_layout()
plt.show(block = False)

# close it?
close = 0
close = input('type literally anything and hit enter to close the window \n')
if close != 0:
    plt.close()