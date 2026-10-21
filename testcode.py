gamma = 1.4
specgasconst = 287  #{J/(kg*K)}
import numpy as np # type: ignore
import math
import RequiredEquations as he

def RayleighFlow1(MachIn, PressStatIn, TempStatIn, DensStatIn, HeatIn) -> tuple:
    # All this is to solve for MachOut
    c_p = gamma*specgasconst/(gamma-1)
    totaltemp1 = he.tt(MachIn, TempStatIn)

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
    TempStatOut = he.st(MachOut, totaltemp2)
    def value1(Mach):
        return((1+gamma)/(1+gamma*Mach**2))
    
    press_star = PressStatIn/value1(MachIn)
    PressStatOut = press_star * value1(MachOut)
    dens_star = DensStatIn * value1(MachIn)*MachIn**2
    DensStatOut = dens_star / (value1(MachOut)*MachOut**2)

    return (MachOut, TempStatOut, PressStatOut, DensStatOut)


speed = 76.2 #m/s
press = 20 #psia
temp  = 294.261 #k
mach = he.machnumber(speed, temp)[0]
heatadd = 2326*400 #j/kg
density = 16.01846 * .1019

print(heatadd)
a = RayleighFlow1(mach, press, temp, density, heatadd)
print(a)



# def RayleighFlow(MachIn, PressStatIn, TempStatIn, DensStatIn, HeatIn, MassFlow) -> tuple:
#     # Simplify Stuff
#     def value(Mach):
#         return((1+gamma)/(1+gamma*Mach**2))
        
#     # Calculate Star Values
#     p_star = PressStatIn / value(MachIn)
#     t_star = TempStatIn / (value(MachIn)*MachIn)**2
#     d_star = DensStatIn * value(MachIn) * MachIn**2

#     # Calculate MachOut
#     # Turn it into a quadratic equation
#     c = HeatIn/MassFlow
#     c1 = c*(gamma-1)/(gamma*specgasconst*TempStatIn)
#     c2 = c1 + 1 + .5*(gamma-1)*MachIn**2
#     d = ((1+gamma*MachIn**2)/MachIn)**2
#     p0 = (.5*(gamma-1)*d - gamma**2)
#     p1 = (d - 2*gamma)
#     p2 = -1*c2
#     coeffs = [p0,p1,p2]
#     roots = np.roots(coeffs)
#     msquared = 0
#     for x in roots:
#         if 1e-5 <= x <= 1:
#             msquared = x
#     MachOut = math.sqrt(msquared)
#     return MachOut