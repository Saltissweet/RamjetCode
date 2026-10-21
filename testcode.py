# Helpful Functions Calculator
# Currently Includes: Oblique Shock Relations, Total Pressure Calculator, Static Pressure Calculator, Normal Shock Relations

import math
from scipy.optimize import brentq # type: ignore
import numpy as np # type: ignore

# global variables
# Canino said gamma is like 1.3 when you increase pressure
gamma = 1.4
specgasconst = 287  #{J/(kg*K)}

"""ISENTROPIC EQUATION CALCULATORS"""

# Isentropic Total Pressure Calculator
# Inputs: Mach Number, Static Pressure
# Outputs: Total Pressure
def TotalPressure(Mach, static_pressure):
    equation = (1 + ((gamma-1)/2) * Mach**2)
    top = (gamma/(gamma-1))
    totalpressure = equation**top * static_pressure
    return (totalpressure)

# Isentropic Static Pressure Calculator
# Inputs: Mach Number, Total Pressure
# Outputs: Static Pressure
def StaticPressure(Mach, total_pressure):
    equation = (1 + ((gamma-1)/2) * Mach**2)
    top = (gamma/(gamma-1))
    static_pressure = total_pressure / equation**top
    return (static_pressure)

# Isentropic Throat Area Calculator
# Inputs: Mach Number, Current Area
# Outputs: Throat Area
def FindThroatArea(MachIn, Area):
    eqbase = (2/(gamma+1)*(1+(gamma-1)/2 * MachIn**2))
    eqexponent = .5*(gamma+1)/(gamma-1)
    arearatioT = (1/MachIn) * eqbase**eqexponent
    throatarea = Area/arearatioT
    return throatarea

# Isentropic Wanted Area Calculator
# Inputs: Wanted Mach Number Out, Throat Area
# Outputs: Wanted Area
def FindWantedArea(MachOut, ThroatArea):
    eqbase = (2/(gamma+1)) * (1 + ((gamma-1)/2)*MachOut**2)
    eqexponent = .5*((gamma+1)/(gamma-1))
    arearatioW = (1/MachOut) * eqbase ** eqexponent
    wantedarea = ThroatArea*arearatioW
    return wantedarea

# Isentropic Total Temp Calculator
# Inputs: Mach Number, Static Temp
# Outputs: Total Temp
def TotalTemp(Mach, statictemp):
    equation = (1 + ((gamma-1)/2) * Mach**2)
    totaltemp = statictemp * equation
    return (totaltemp)

# Isentropic Static Temp Calculator
# Inputs: Mach Number, Total Temp
# Outputs: Static Temp
def StaticTemp(Mach, totaltemp):
    equation = (1 + ((gamma-1)/2) * Mach**2)
    statictemp = totaltemp / equation
    return (statictemp)

# Isentropic Total Density Calculator
# Inputs: Mach Number, Static Density
# Outputs: Total Density
def TotalDensity(Mach, static_density):
    equation = (1 + ((gamma-1)/2) * Mach**2)
    top = (1/(gamma-1))
    totaldensity = equation**top * static_density
    return (totaldensity)

# Isentropic Static Density Calculator
# Inputs: Mach Number, Total Pressure
# Outputs: Static Density
def StaticDensity(Mach, total_density):
    equation = (1 + ((gamma-1)/2) * Mach**2)
    top = (1/(gamma-1))
    static_density = total_density / equation**top
    return (static_density)

def FindMachPress(presstot, press):
    frac = (presstot/press)**((gamma-1)/gamma)
    Mach = math.sqrt((frac-1)*2/(gamma-1))
    return Mach

# Calculate Speed of Flow, and speed of sound
def speed(Mach, Temperature):
    a = math.sqrt(Temperature * gamma * specgasconst)
    v = a * Mach                  
    return (v, a) 

def massflow(speed, density, area):
    return(speed * density * area)

def machnumber(Speed, Temperature):
    a = math.sqrt(Temperature * gamma * specgasconst)
    mach = Speed/a
    return(mach, a)



"""COMPLICATED EQUATIONS"""

# Oblique Shock Relation Function
# Inputs: Mach in, Beta, Pressure Total in, Temperature in
# Outputs: Mach Out, Theta, Pressure Total out, Temperature out
def ObliqueShock(MachIn, BetaRelative, PressTotIn, TempIn, DensIn) -> tuple:
    BetaR = math.radians(BetaRelative)

    # Calculate Wedge Angle
    part1eq = (MachIn**2 * (math.sin(BetaR)**2)-1)/(MachIn**2 * (gamma + math.cos(2*BetaR))+2)
    tantheta = (2/math.tan(BetaR))*part1eq
    theta = math.atan(tantheta)

    # Calculate Mach Out
    part2eq = ((gamma-1) * MachIn**2 * math.sin(BetaR)**2 + 2)/(2*gamma*MachIn**2*math.sin(BetaR)**2 - (gamma-1))
    MachOut = math.sqrt(part2eq/math.sin(BetaR-theta)**2)

    # Calculate Pto out
    part3_1eq = (((gamma+1)*MachIn**2*math.sin(BetaR)**2)/((gamma-1)*MachIn**2*math.sin(BetaR)**2 + 2))**(gamma/(gamma-1))
    part3_2eq = ((gamma+1)/(2*gamma*MachIn**2*math.sin(BetaR)**2 - (gamma-1)))**(1/(gamma-1))
    pt1_pt2 = part3_1eq * part3_2eq
    PressTotOut = pt1_pt2 * PressTotIn

    # Calculate T out
    part4_topeq = (2*gamma*MachIn**2 * math.sin(BetaR)**2 - (gamma-1))*((gamma-1)*MachIn**2*math.sin(BetaR)**2 + 2)
    part4_botteq = (gamma+1)**2*MachIn**2*math.sin(BetaR)**2
    t1_t2 = part4_topeq/part4_botteq
    TempOut = t1_t2 * TempIn

    # Calculate Rho out
    part5_eq = ((gamma+1)*MachIn**2 * math.sin(BetaR)**2)/((gamma-1)*MachIn**2 * math.sin(BetaR)**2 + 2)
    DensOut = DensIn*part5_eq

    return (MachOut, math.degrees(theta), PressTotOut, TempOut, DensOut)


# Normal Shock Calculator
# Inputs: Mach in, Pressure Total in, Temperature in
# Outputs: Mach Out, Pressure Total out, Temperature out
def NormalShock(MachIn, PressTotIn, TempIn, DensIn) -> tuple:
    # Calculate Mach Out
    MachOutSq = ((gamma-1)*MachIn**2+2)/(2*gamma*MachIn**2 -(gamma-1))
    MachOut = math.sqrt(MachOutSq)

    # Calculate P tot out
    p_eq_1 = (((gamma+1)*MachIn**2)/((gamma-1)*MachIn**2+2))**(gamma/(gamma-1))
    p_eq_2 = ((gamma+1)/(2*gamma*MachIn**2-(gamma-1)))**(1/(gamma-1))
    PressTotOut = p_eq_1*p_eq_2*PressTotIn

    # Calculate T out
    t_eq_top = (2*gamma*MachIn**2 - (gamma-1))*((gamma-1)*MachIn**2+2)
    t_eq_bott = (gamma+1)**2*MachIn**2
    TempOut = TempIn*t_eq_top/t_eq_bott

    d_eq = ((gamma+1)*MachIn**2)/((gamma-1)*MachIn**2 + 2)
    DensOut = DensIn*d_eq

    return (MachOut, PressTotOut, TempOut, DensOut)


# SUBSONIC RAYLIEGH FLOW, for supersonic change the bounds of the if statement
# Inputs: Mach In, Static Pressure In, Static Temperature In, Static Density In, Heat Added
# Outputs: Mach Out, Static Pressure Out, Static Temperature Out, Static Density Out
def RayleighFlow(MachIn, PressStatIn, TempStatIn, DensStatIn, HeatIn) -> tuple:
    # Canino said gamma is like 1.3 when you increase pressure
    # Fanno flow is basically only used in combustor
    gamma = 1.3
    # All this is to solve for MachOut
    c_p = gamma*specgasconst/(gamma-1)
    totaltemp1 = TotalTemp(MachIn, TempStatIn)

    value = (1+gamma)/(1+gamma*MachIn**2)**2

    temp_tot_star = totaltemp1 / (value*MachIn**2 * (2+(gamma-1)*MachIn**2))
    tempadd = HeatIn/c_p
    totaltemp2 = totaltemp1 + tempadd
    fv = totaltemp2/temp_tot_star # fracval
    p0 = (gamma+1)*(gamma-1) - fv*gamma**2
    p1 = 2*(gamma+1 - gamma*fv)
    p2 = -1*fv
    coeffs = [p0,p1,p2]
    roots = np.roots(coeffs)
    msquared = -1
    for x in roots:
        if 1e-5 <= x <= 1:
            msquared = x
    MachOut = math.sqrt(msquared)
    TempStatOut = StaticTemp(MachOut, totaltemp2)
    def value1(Mach):
        return((1+gamma)/(1+gamma*Mach**2))
    
    press_star = PressStatIn/value1(MachIn)
    PressStatOut = press_star * value1(MachOut)
    dens_star = DensStatIn * value1(MachIn)*MachIn**2
    DensStatOut = dens_star / (value1(MachOut)*MachOut**2)

    return (MachOut, PressStatOut, TempStatOut, DensStatOut)



# SUBSONIC FANNO FLOW, for supersonic change the bounds of Brents stuff
# Inputs: Mach In, Static Pressure In, Static Temperature In, Static Density In, Friction Factor (fbar), Distance, Diameter
# Outputs: Mach Out, Static Pressure Out, Static Temperature Out, Static Density Out
def FannoFlow(MachIn, PressStatIn, TempStatIn, DensStatIn, FrictFact, Dist, Diam) -> tuple:
    # Canino said gamma is like 1.3 when you increase pressure
    # Fanno flow is basically only used in combustor/past the inlet so its hot flow
    gamma = 1.3
    # Simplify Calcs
    def value(Mach):
        return ((gamma+1)/(2+(gamma-1)*Mach**2))
    
    p_star = PressStatIn * MachIn / math.sqrt(value(MachIn))
    t_star = TempStatIn / value(MachIn)
    d_star = DensStatIn * MachIn * math.sqrt(value(MachIn))\

    # Solve For l star
    a = (1 - MachIn**2)/(gamma*MachIn**2)
    b = ((gamma+1)/(2*gamma))*math.log(value(MachIn)*MachIn**2)

    # Solve for final Right Hand Side Value
    val1 = a+b
    diststar1 = val1 * Diam/FrictFact
    diststar2 = diststar1 - Dist
    val2 = diststar2*FrictFact/Diam

    # Find MachOut Using Brents
    def BrentFuncResidual(Mach):
        a1 = (1 - Mach**2)/(gamma*Mach**2)
        b1 = ((gamma+1)/(2*gamma))*math.log(value(Mach)*Mach**2)
        return (a1+b1-val2)
    
    MachOut, result = brentq(BrentFuncResidual, 1e-5, 1.0, full_output=True)
    if not result.converged:
        raise RuntimeError("Root finder failed to converge.")

    PressStatOut = p_star * math.sqrt(value(MachOut)) / MachOut
    TempStatOut = t_star * value(MachOut)
    DensStatOut = d_star/(math.sqrt(value(MachOut))*MachOut)

    return (MachOut, PressStatOut, TempStatOut, DensStatOut)




# Faster Call equations lmao
def osr(MachIn, BetaRelative, PressTotIn, TempTotIn, DensIn):
    return ObliqueShock(MachIn, BetaRelative, PressTotIn, TempTotIn, DensIn)

def nsr(MachIn, PressTotIn, TempIn, DensIn):
    return NormalShock(MachIn, PressTotIn, TempIn, DensIn)

def tp(Mach, static_pressure):
    return TotalPressure(Mach, static_pressure)

def sp(Mach, total_pressure):
    return StaticPressure(Mach, total_pressure)

def fta(MachIn, Area):
    return FindThroatArea(MachIn, Area)

def fwa(MachOut, ThroatArea):
    return FindWantedArea(MachOut, ThroatArea)

def tt(Mach, static_temp):
    return TotalTemp(Mach, static_temp)

def st(Mach, total_temp):
    return StaticTemp(Mach, total_temp)

def td(Mach, static_density):
    return TotalDensity(Mach, static_density)

def sd(Mach, total_density):
    return StaticDensity(Mach, total_density)

def ff(MachIn, PressStatIn, TempStatIn, DensStatIn, FrictFact, Dist, Diam):
    return FannoFlow(MachIn, PressStatIn, TempStatIn, DensStatIn, FrictFact, Dist, Diam)

def rf(MachIn, PressStatIn, TempStatIn, DensStatIn, HeatIn):
    return RayleighFlow(MachIn, PressStatIn, TempStatIn, DensStatIn, HeatIn)

def fmp(presstot, press):
    return FindMachPress(presstot, press)