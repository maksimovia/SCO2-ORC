import pandas as pd
import prop
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
streams = pd.read_excel('streams.xlsx', index_col=0, sheet_name='reg-comp')

Tmax = 600
Pk = 7.8
pi=3
KPDcomp = 0.9
KPDturb = 0.9
Tmin = 32
dTreg = 80
Gas = 'CO2'
G0 = 1
#cool
streams.at['Cool-Comp1', 'T'] = Tmin
streams.at['Cool-Comp1', 'P'] = Pk
streams.at['Cool-Comp1', 'X'] = "CO2"
streams.at['Cool-Comp1', 'H'] = prop.t_p(streams.at['Cool-Comp1', 'T'],streams.at['Cool-Comp1', 'P'],streams.at['Cool-Comp1', 'X'])['H']
streams.at['Cool-Comp1', 'G'] = G0
streams.at['Cool-Comp1', 'Q'] = 1
streams.at['Cool-Comp1', 'S'] = prop.h_p(streams.at['Cool-Comp1', 'H'],streams.at['Cool-Comp1', 'P'],streams.at['Cool-Comp1', 'X'])['S']

#heat
streams.at['Heat-Turb', 'T'] = Tmax
streams.at['Heat-Turb', 'P'] = Pk*pi
streams.at['Heat-Turb', 'X'] = Gas
streams.at['Heat-Turb', 'G'] = G0
streams.at['Heat-Turb', 'H'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['H']
streams.at['Heat-Turb', 'Q'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['Q']
streams.at['Heat-Turb', 'S'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['S']

#turb
streams.at['Turb-Reg', 'P'] = Pk
streams.at['Turb-Reg', 'X'] = Gas
streams.at['Turb-Reg', 'G'] = G0
streams.at['Turb-Reg', 'T'] = Tmax
Ht = prop.p_s(streams.at['Turb-Reg', 'P'], streams.at['Heat-Turb', 'S'], streams.at['Turb-Reg', 'X'])['H']
streams.at['Turb-Reg', 'H'] = streams.at['Heat-Turb', 'H'] - (streams.at['Heat-Turb', 'H']-Ht)*KPDturb
streams.at['Turb-Reg', 'T'] = prop.h_p(streams.at['Turb-Reg', 'H'],streams.at['Turb-Reg', 'P'],streams.at['Turb-Reg', 'X'])["T"]
streams.at['Turb-Reg', 'Q'] = prop.h_p(streams.at['Turb-Reg', 'H'],streams.at['Turb-Reg', 'P'],streams.at['Turb-Reg', 'X'])["Q"]
streams.at['Turb-Reg', 'S'] = prop.h_p(streams.at['Turb-Reg', 'H'],streams.at['Turb-Reg', 'P'],streams.at['Turb-Reg', 'X'])["S"]

#comp1
streams.at['Comp1-Evap1', 'P'] = Pk*pi**0.5
streams.at['Comp1-Evap1', 'X'] = Gas
streams.at['Comp1-Evap1', 'G'] = G0
Ht = prop.p_s(streams.at['Comp1-Evap1', 'P'], streams.at['Cool-Comp1', 'S'], streams.at['Comp1-Evap1', 'X'])['H']
streams.at['Comp1-Evap1', 'H'] = streams.at['Cool-Comp1', 'H'] + (Ht - streams.at['Cool-Comp1', 'H'])/KPDcomp
streams.at['Comp1-Evap1', 'T'] = prop.h_p(streams.at['Comp1-Evap1', 'H'],streams.at['Comp1-Evap1', 'P'],streams.at['Comp1-Evap1', 'X'])["T"]
streams.at['Comp1-Evap1', 'Q'] = prop.h_p(streams.at['Comp1-Evap1', 'H'],streams.at['Comp1-Evap1', 'P'],streams.at['Comp1-Evap1', 'X'])["Q"]
streams.at['Comp1-Evap1', 'S'] = prop.h_p(streams.at['Comp1-Evap1', 'H'],streams.at['Comp1-Evap1', 'P'],streams.at['Comp1-Evap1', 'X'])["S"]

Tcomp1 = 44

#evap1
streams.at['Evap1-Comp2', 'P'] = Pk*pi**0.5
streams.at['Evap1-Comp2', 'X'] = Gas
streams.at['Evap1-Comp2', 'G'] = G0
streams.at['Evap1-Comp2', 'T'] = Tcomp1
streams.at['Evap1-Comp2', 'H'] = prop.t_p(streams.at['Evap1-Comp2', 'T'],streams.at['Evap1-Comp2', 'P'],streams.at['Evap1-Comp2', 'X'])['H']
streams.at['Evap1-Comp2', 'Q'] = prop.t_p(streams.at['Evap1-Comp2', 'T'],streams.at['Evap1-Comp2', 'P'],streams.at['Evap1-Comp2', 'X'])['Q']
streams.at['Evap1-Comp2', 'S'] = prop.t_p(streams.at['Evap1-Comp2', 'T'],streams.at['Evap1-Comp2', 'P'],streams.at['Evap1-Comp2', 'X'])['S']

#comp2
streams.at['Comp2-Reg', 'P'] = Pk*pi
streams.at['Comp2-Reg', 'X'] = Gas
streams.at['Comp2-Reg', 'G'] = G0
Ht = prop.p_s(streams.at['Comp2-Reg', 'P'], streams.at['Evap1-Comp2', 'S'], streams.at['Comp2-Reg', 'X'])['H']
streams.at['Comp2-Reg', 'H'] = streams.at['Evap1-Comp2', 'H'] + (Ht - streams.at['Evap1-Comp2', 'H'])/KPDcomp
streams.at['Comp2-Reg', 'T'] = prop.h_p(streams.at['Comp2-Reg', 'H'],streams.at['Comp2-Reg', 'P'],streams.at['Comp2-Reg', 'X'])["T"]
streams.at['Comp2-Reg', 'Q'] = prop.h_p(streams.at['Comp2-Reg', 'H'],streams.at['Comp2-Reg', 'P'],streams.at['Comp2-Reg', 'X'])["Q"]
streams.at['Comp2-Reg', 'S'] = prop.h_p(streams.at['Comp2-Reg', 'H'],streams.at['Comp2-Reg', 'P'],streams.at['Comp2-Reg', 'X'])["S"]

#reg hot
streams.at['Reg-Evap2', 'T'] = streams.at['Comp2-Reg', 'T']+dTreg
streams.at['Reg-Evap2', 'P'] = streams.at['Turb-Reg', 'P']
streams.at['Reg-Evap2', 'X'] = streams.at['Turb-Reg', 'X']
streams.at['Reg-Evap2', 'G'] = streams.at['Turb-Reg', 'G']
streams.at['Reg-Evap2', 'H'] = prop.t_p(streams.at['Reg-Evap2', 'T'],streams.at['Reg-Evap2', 'P'],streams.at['Reg-Evap2', 'X'])['H']
streams.at['Reg-Evap2', 'Q'] = prop.t_p(streams.at['Reg-Evap2', 'T'],streams.at['Reg-Evap2', 'P'],streams.at['Reg-Evap2', 'X'])['Q']
streams.at['Reg-Evap2', 'S'] = prop.t_p(streams.at['Reg-Evap2', 'T'],streams.at['Reg-Evap2', 'P'],streams.at['Reg-Evap2', 'X'])['S']

#reg cold
streams.at['Reg-Heat', 'P'] = streams.at['Comp2-Reg', 'P']
streams.at['Reg-Heat', 'X'] = streams.at['Comp2-Reg', 'X']
streams.at['Reg-Heat', 'G'] = streams.at['Comp2-Reg', 'G']
streams.at['Reg-Heat', 'H'] = streams.at['Comp2-Reg', 'H']+(streams.at['Turb-Reg', 'H'] - streams.at['Reg-Evap2', 'H'])
streams.at['Reg-Heat', 'T'] = prop.h_p(streams.at['Reg-Heat', 'H'],streams.at['Reg-Heat', 'P'],streams.at['Reg-Heat', 'X'])['T']
streams.at['Reg-Heat', 'Q'] = prop.h_p(streams.at['Reg-Heat', 'H'],streams.at['Reg-Heat', 'P'],streams.at['Reg-Heat', 'X'])['Q']
streams.at['Reg-Heat', 'S'] = prop.h_p(streams.at['Reg-Heat', 'H'],streams.at['Reg-Heat', 'P'],streams.at['Reg-Heat', 'X'])['S']


Tmin_orc = 30
Fluid_orc = "R236ea"
Pk_orc = prop.t_q(Tmin_orc,0,Fluid_orc)["P"]

#ORC Cool
streams.at['OCool-OPump', 'P'] = Pk_orc
streams.at['OCool-OPump', 'T'] = Tmin_orc
streams.at['OCool-OPump', 'X'] = Fluid_orc
streams.at['OCool-OPump', 'H'] = prop.t_p(streams.at['OCool-OPump', 'T'],streams.at['OCool-OPump', 'P'],streams.at['OCool-OPump', 'X'])["H"]
streams.at['OCool-OPump', 'Q'] = prop.t_p(streams.at['OCool-OPump', 'T'],streams.at['OCool-OPump', 'P'],streams.at['OCool-OPump', 'X'])["Q"]
streams.at['OCool-OPump', 'S'] = prop.t_p(streams.at['OCool-OPump', 'T'],streams.at['OCool-OPump', 'P'],streams.at['OCool-OPump', 'X'])["S"]

pi_orc = 4
KPDpump = KPDcomp
#ORC Pump
streams.at['OPump-OEvap1', 'P'] = streams.at['OCool-OPump', 'P']*pi_orc
streams.at['OPump-OEvap1', 'X'] = streams.at['OCool-OPump', 'X']
Ht = prop.p_s(streams.at['OPump-OEvap1', 'P'], streams.at['OCool-OPump', 'S'], streams.at['OPump-OEvap1', 'X'])['H']
streams.at['OPump-OEvap1', 'H'] = streams.at['OCool-OPump', 'H'] + (Ht - streams.at['OCool-OPump', 'H'])/KPDpump
streams.at['OPump-OEvap1', 'T'] = prop.h_p(streams.at['OPump-OEvap1', 'H'],streams.at['OPump-OEvap1', 'P'],streams.at['OPump-OEvap1', 'X'])["T"]
streams.at['OPump-OEvap1', 'Q'] = prop.h_p(streams.at['OPump-OEvap1', 'H'],streams.at['OPump-OEvap1', 'P'],streams.at['OPump-OEvap1', 'X'])["Q"]
streams.at['OPump-OEvap1', 'S'] = prop.h_p(streams.at['OPump-OEvap1', 'H'],streams.at['OPump-OEvap1', 'P'],streams.at['OPump-OEvap1', 'X'])["S"]

H22 = prop.p_q(streams.at['OPump-OEvap1', 'P'],1,Fluid_orc)['H']
H21 = prop.p_q(streams.at['OPump-OEvap1', 'P'],0,Fluid_orc)['H']
H11 = streams.at['Reg-Evap2', 'H']
T21 = prop.h_p(H21,streams.at['OPump-OEvap1', 'P'],Fluid_orc)['T']
T12 = T21 + 15
H12 = prop.t_p(T12,streams.at['Reg-Evap2', 'P'],Gas)['H']
G_orc = (H11-H12)/(H22-H21)*G0
print(T21)
#ORC Evap 1
streams.at['OEvap1-OEvap2', 'P'] = streams.at['OPump-OEvap1', 'P']
streams.at['OEvap1-OEvap2', 'H'] = (G0*(streams.at['Comp1-Evap1', 'H']-streams.at['Evap1-Comp2', 'H']))/G_orc + streams.at['OPump-OEvap1', 'H']
streams.at['OEvap1-OEvap2', 'X'] = "R236ea"
streams.at['OEvap1-OEvap2', 'G'] = G_orc
streams.at['OEvap1-OEvap2', 'T'] = prop.h_p(streams.at['OEvap1-OEvap2', 'H'],streams.at['OEvap1-OEvap2', 'P'],streams.at['OEvap1-OEvap2', 'X'])["T"]
streams.at['OEvap1-OEvap2', 'Q'] = prop.h_p(streams.at['OEvap1-OEvap2', 'H'],streams.at['OEvap1-OEvap2', 'P'],streams.at['OEvap1-OEvap2', 'X'])["Q"]
streams.at['OEvap1-OEvap2', 'S'] = prop.h_p(streams.at['OEvap1-OEvap2', 'H'],streams.at['OEvap1-OEvap2', 'P'],streams.at['OEvap1-OEvap2', 'X'])["S"]

H23 = streams.at['OEvap1-OEvap2', 'H']
H13 = H12 - (H21 - H23)*G_orc/G0

#evap2
streams.at['Evap2-Cool', 'H'] = H13
streams.at['Evap2-Cool', 'P'] = Pk
streams.at['Evap2-Cool', 'X'] = streams.at['Reg-Evap2', 'X']
streams.at['Evap2-Cool', 'G'] = streams.at['Reg-Evap2', 'G']
streams.at['Evap2-Cool', 'T'] = prop.h_p(streams.at['Evap2-Cool', 'H'],streams.at['Evap2-Cool', 'P'],streams.at['Evap2-Cool', 'X'])["T"]
streams.at['Evap2-Cool', 'Q'] = prop.h_p(streams.at['Evap2-Cool', 'H'],streams.at['Evap2-Cool', 'P'],streams.at['Evap2-Cool', 'X'])["Q"]
streams.at['Evap2-Cool', 'S'] = prop.h_p(streams.at['Evap2-Cool', 'H'],streams.at['Evap2-Cool', 'P'],streams.at['Evap2-Cool', 'X'])["S"]

#ORC Evap2
streams.at['OEvap2-OTurb', 'P'] = streams.at['OEvap1-OEvap2', 'P']
streams.at['OEvap2-OTurb', 'H'] = (G0*(streams.at['Reg-Evap2', 'H']-streams.at['Evap2-Cool', 'H']))/G_orc + streams.at['OEvap1-OEvap2', 'H']
streams.at['OEvap2-OTurb', 'X'] = "R236ea"
streams.at['OEvap2-OTurb', 'G'] = G_orc
streams.at['OEvap2-OTurb', 'T'] = prop.h_p(streams.at['OEvap2-OTurb', 'H'],streams.at['OEvap2-OTurb', 'P'],streams.at['OEvap2-OTurb', 'X'])["T"]
streams.at['OEvap2-OTurb', 'Q'] = prop.h_p(streams.at['OEvap2-OTurb', 'H'],streams.at['OEvap2-OTurb', 'P'],streams.at['OEvap2-OTurb', 'X'])["Q"]
streams.at['OEvap2-OTurb', 'S'] = prop.h_p(streams.at['OEvap2-OTurb', 'H'],streams.at['OEvap2-OTurb', 'P'],streams.at['OEvap2-OTurb', 'X'])["S"]

n=50
T_co2_evap1 = np.zeros(n)
H_co2_evap1 = np.linspace(streams.at['Comp1-Evap1', 'H'], streams.at['Evap1-Comp2', 'H'], n)
for i in range(0, n):
    T_co2_evap1[i] = prop.h_p(H_co2_evap1[i], streams.at['Evap1-Comp2', 'P'], Gas)["T"]
Q_co2_evap1 = -(H_co2_evap1-streams.at['Comp1-Evap1', 'H'])*streams.at['Evap1-Comp2', 'G']

T_orc_evap1 = np.zeros(n)
H_orc_evap1 = np.linspace(streams.at['OEvap1-OEvap2', 'H'], streams.at['OPump-OEvap1', 'H'], n)
for i in range(0, n):
    T_orc_evap1[i] = prop.h_p(H_orc_evap1[i], streams.at['OPump-OEvap1', 'P'], streams.at['OPump-OEvap1', 'X'])["T"]
Q_orc_evap1 = -(H_orc_evap1-streams.at['OEvap1-OEvap2', 'H'])*streams.at['OEvap1-OEvap2', 'G']
# plt.plot(Q_co2_evap1,T_co2_evap1)
# plt.plot(Q_orc_evap1,T_orc_evap1)
# plt.show()

n=50
T_co2_evap2 = np.zeros(n)
H_co2_evap2 = np.linspace(streams.at['Reg-Evap2', 'H'], streams.at['Evap2-Cool', 'H'], n)
for i in range(0, n):
    T_co2_evap2[i] = prop.h_p(H_co2_evap2[i], streams.at['Evap2-Cool', 'P'], Gas)["T"]
Q_co2_evap2 = -(H_co2_evap2-streams.at['Reg-Evap2', 'H'])*streams.at['Evap2-Cool', 'G']

T_orc_evap2 = np.zeros(n)
H_orc_evap2 = np.linspace(streams.at['OEvap2-OTurb', 'H'], streams.at['OEvap1-OEvap2', 'H'], n)
for i in range(0, n):
    T_orc_evap2[i] = prop.h_p(H_orc_evap2[i], streams.at['OEvap1-OEvap2', 'P'], streams.at['OEvap1-OEvap2', 'X'])["T"]
Q_orc_evap2 = -(H_orc_evap2-streams.at['OEvap2-OTurb', 'H'])*streams.at['OEvap2-OTurb', 'G']
plt.plot(Q_co2_evap2,T_co2_evap2,color="red")
plt.plot(Q_orc_evap2,T_orc_evap2,color="orange")
plt.plot(Q_co2_evap1+Q_co2_evap2[-1],T_co2_evap1,color="red")
plt.plot(Q_orc_evap1+Q_orc_evap2[-1],T_orc_evap1,color="orange")
plt.show()


print(streams)
print(G_orc,H11,H12,H21,H22, H13)
# def CalcPP(pi,pi_orc):
#     #pi = 2
#     #pi_orc = 10
#     def CalcGH(Input):
#         G_orc = Input[0]
#         H_percent = Input[1]
#         Tmax = 600
#         Pk = 7.8
#         #pi=12
#         KPDcomp = 0.9
#         KPDturb = 0.9
#         Tmin = 32
#         dTreg = 50
#
#         #cool
#         streams.at['Cool-Comp', 'T'] = Tmin
#         streams.at['Cool-Comp', 'P'] = Pk
#         streams.at['Cool-Comp', 'X'] = "CO2"
#         streams.at['Cool-Comp', 'H'] = prop.t_p(streams.at['Cool-Comp', 'T'],streams.at['Cool-Comp', 'P'],streams.at['Cool-Comp', 'X'])['H']
#         streams.at['Cool-Comp', 'G'] = 1
#         streams.at['Cool-Comp', 'Q'] = 1
#         streams.at['Cool-Comp', 'S'] = prop.h_p(streams.at['Cool-Comp', 'H'],streams.at['Cool-Comp', 'P'],streams.at['Cool-Comp', 'X'])['S']
#
#         #comp
#         streams.at['Comp-Reg', 'P'] = streams.at['Cool-Comp', 'P']*pi
#         streams.at['Comp-Reg', 'X'] = streams.at['Cool-Comp', 'X']
#         Ht = prop.p_s(streams.at['Comp-Reg', 'P'], streams.at['Cool-Comp', 'S'], streams.at['Comp-Reg', 'X'])['H']
#         streams.at['Comp-Reg', 'H'] = streams.at['Cool-Comp', 'H'] + (Ht - streams.at['Cool-Comp', 'H'])/KPDcomp
#         streams.at['Comp-Reg', 'G'] = streams.at['Cool-Comp', 'G']
#         streams.at['Comp-Reg', 'T'] = prop.h_p(streams.at['Comp-Reg', 'H'],streams.at['Comp-Reg', 'P'],streams.at['Comp-Reg', 'X'])["T"]
#         streams.at['Comp-Reg', 'Q'] = prop.h_p(streams.at['Comp-Reg', 'H'],streams.at['Comp-Reg', 'P'],streams.at['Comp-Reg', 'X'])["Q"]
#         streams.at['Comp-Reg', 'S'] = prop.h_p(streams.at['Comp-Reg', 'H'],streams.at['Comp-Reg', 'P'],streams.at['Comp-Reg', 'X'])["S"]
#
#         #reg_cold
#         streams.at['Reg-Heat', 'P'] = streams.at['Comp-Reg', 'P']
#         streams.at['Reg-Heat', 'X'] = streams.at['Comp-Reg', 'X']
#         streams.at['Reg-Heat', 'G'] = streams.at['Comp-Reg', 'G']
#
#
#         # #heat
#         streams.at['Heat-Turb', 'T'] = Tmax
#         streams.at['Heat-Turb', 'P'] = streams.at['Reg-Heat', 'P']
#         streams.at['Heat-Turb', 'X'] = streams.at['Reg-Heat', 'X']
#         streams.at['Heat-Turb', 'G'] = streams.at['Reg-Heat', 'G']
#         streams.at['Heat-Turb', 'H'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['H']
#         streams.at['Heat-Turb', 'Q'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['Q']
#         streams.at['Heat-Turb', 'S'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['S']
#
#         #turb
#         streams.at['Turb-Reg', 'P'] = Pk
#         streams.at['Turb-Reg', 'X'] = streams.at['Heat-Turb', 'X']
#         streams.at['Turb-Reg', 'G'] = streams.at['Heat-Turb', 'G']
#         Ht = prop.p_s(streams.at['Turb-Reg', 'P'], streams.at['Heat-Turb', 'S'], streams.at['Turb-Reg', 'X'])['H']
#         streams.at['Turb-Reg', 'H'] = streams.at['Heat-Turb', 'H'] - (streams.at['Heat-Turb', 'H']-Ht)*KPDturb
#         streams.at['Turb-Reg', 'T'] = prop.h_p(streams.at['Turb-Reg', 'H'],streams.at['Turb-Reg', 'P'],streams.at['Turb-Reg', 'X'])["T"]
#         streams.at['Turb-Reg', 'Q'] = prop.h_p(streams.at['Turb-Reg', 'H'],streams.at['Turb-Reg', 'P'],streams.at['Turb-Reg', 'X'])["Q"]
#         streams.at['Turb-Reg', 'S'] = prop.h_p(streams.at['Turb-Reg', 'H'],streams.at['Turb-Reg', 'P'],streams.at['Turb-Reg', 'X'])["S"]
#
#         #reg_hot
#         streams.at['Reg-Evap', 'T'] = streams.at['Comp-Reg', 'T']+dTreg
#         streams.at['Reg-Evap', 'P'] = streams.at['Turb-Reg', 'P']
#         streams.at['Reg-Evap', 'X'] = streams.at['Turb-Reg', 'X']
#         streams.at['Reg-Evap', 'G'] = streams.at['Turb-Reg', 'G']
#         streams.at['Reg-Evap', 'H'] = prop.t_p(streams.at['Reg-Evap', 'T'],streams.at['Reg-Evap', 'P'],streams.at['Reg-Evap', 'X'])['H']
#         streams.at['Reg-Evap', 'Q'] = prop.t_p(streams.at['Reg-Evap', 'T'],streams.at['Reg-Evap', 'P'],streams.at['Reg-Evap', 'X'])['Q']
#         streams.at['Reg-Evap', 'S'] = prop.t_p(streams.at['Reg-Evap', 'T'],streams.at['Reg-Evap', 'P'],streams.at['Reg-Evap', 'X'])['S']
#
#         #reg_cold
#         streams.at['Reg-Heat', 'H'] = streams.at['Comp-Reg', 'H']+(streams.at['Turb-Reg', 'H'] - streams.at['Reg-Evap', 'H'])
#         streams.at['Reg-Heat', 'T'] = prop.h_p(streams.at['Reg-Heat', 'H'],streams.at['Reg-Heat', 'P'],streams.at['Reg-Heat', 'X'])['T']
#         streams.at['Reg-Heat', 'Q'] = prop.h_p(streams.at['Reg-Heat', 'H'],streams.at['Reg-Heat', 'P'],streams.at['Reg-Heat', 'X'])['Q']
#         streams.at['Reg-Heat', 'S'] = prop.h_p(streams.at['Reg-Heat', 'H'],streams.at['Reg-Heat', 'P'],streams.at['Reg-Heat', 'X'])['S']
#
#         #H_percent = 0.6
#         #evap
#         streams.at['Evap-Cool', 'H'] = (streams.at['Reg-Evap', 'H']-streams.at['Cool-Comp', 'H'])*H_percent+streams.at['Cool-Comp', 'H']
#         streams.at['Evap-Cool', 'P'] = Pk
#         streams.at['Evap-Cool', 'X'] = streams.at['Reg-Evap', 'X']
#         streams.at['Evap-Cool', 'G'] = streams.at['Reg-Evap', 'G']
#         streams.at['Evap-Cool', 'T'] = prop.h_p(streams.at['Evap-Cool', 'H'],streams.at['Evap-Cool', 'P'],streams.at['Evap-Cool', 'X'])["T"]
#         streams.at['Evap-Cool', 'Q'] = prop.h_p(streams.at['Evap-Cool', 'H'],streams.at['Evap-Cool', 'P'],streams.at['Evap-Cool', 'X'])["Q"]
#         streams.at['Evap-Cool', 'S'] = prop.h_p(streams.at['Evap-Cool', 'H'],streams.at['Evap-Cool', 'P'],streams.at['Evap-Cool', 'X'])["S"]
#
#
#         Tmin_orc = 30
#         Fluid_orc = "R236ea"
#         Pk_orc = prop.t_q(Tmin_orc,0,Fluid_orc)["P"]
#         #G_orc = 0.66
#         #ORC Cool
#         streams.at['OCool-OPump', 'P'] = Pk_orc
#         streams.at['OCool-OPump', 'T'] = Tmin_orc
#         streams.at['OCool-OPump', 'X'] = Fluid_orc
#         streams.at['OCool-OPump', 'G'] = G_orc
#         streams.at['OCool-OPump', 'H'] = prop.t_p(streams.at['OCool-OPump', 'T'],streams.at['OCool-OPump', 'P'],streams.at['OCool-OPump', 'X'])["H"]
#         streams.at['OCool-OPump', 'Q'] = prop.t_p(streams.at['OCool-OPump', 'T'],streams.at['OCool-OPump', 'P'],streams.at['OCool-OPump', 'X'])["Q"]
#         streams.at['OCool-OPump', 'S'] = prop.t_p(streams.at['OCool-OPump', 'T'],streams.at['OCool-OPump', 'P'],streams.at['OCool-OPump', 'X'])["S"]
#
#         #pi_orc = 12
#         KPDpump = KPDcomp
#         #ORC Pump
#         streams.at['OPump-OEvap', 'P'] = streams.at['OCool-OPump', 'P']*pi_orc
#         streams.at['OPump-OEvap', 'X'] = streams.at['OCool-OPump', 'X']
#         Ht = prop.p_s(streams.at['OPump-OEvap', 'P'], streams.at['OCool-OPump', 'S'], streams.at['OPump-OEvap', 'X'])['H']
#         streams.at['OPump-OEvap', 'H'] = streams.at['OCool-OPump', 'H'] + (Ht - streams.at['OCool-OPump', 'H'])/KPDpump
#         streams.at['OPump-OEvap', 'G'] = streams.at['OCool-OPump', 'G']
#         streams.at['OPump-OEvap', 'T'] = prop.h_p(streams.at['OPump-OEvap', 'H'],streams.at['OPump-OEvap', 'P'],streams.at['OPump-OEvap', 'X'])["T"]
#         streams.at['OPump-OEvap', 'Q'] = prop.h_p(streams.at['OPump-OEvap', 'H'],streams.at['OPump-OEvap', 'P'],streams.at['OPump-OEvap', 'X'])["Q"]
#         streams.at['OPump-OEvap', 'S'] = prop.h_p(streams.at['OPump-OEvap', 'H'],streams.at['OPump-OEvap', 'P'],streams.at['OPump-OEvap', 'X'])["S"]
#
#         #ORC Evap
#         streams.at['OEvap-OTurb', 'P'] = streams.at['OPump-OEvap', 'P']
#         streams.at['OEvap-OTurb', 'H'] = (streams.at['Cool-Comp', 'G']*(streams.at['Reg-Evap', 'H']-streams.at['Evap-Cool', 'H']))/G_orc + streams.at['OPump-OEvap', 'H']
#         streams.at['OEvap-OTurb', 'X'] = "R236ea"
#         streams.at['OEvap-OTurb', 'G'] = G_orc
#         streams.at['OEvap-OTurb', 'T'] = prop.h_p(streams.at['OEvap-OTurb', 'H'],streams.at['OEvap-OTurb', 'P'],streams.at['OEvap-OTurb', 'X'])["T"]
#         streams.at['OEvap-OTurb', 'Q'] = prop.h_p(streams.at['OEvap-OTurb', 'H'],streams.at['OEvap-OTurb', 'P'],streams.at['OEvap-OTurb', 'X'])["Q"]
#         streams.at['OEvap-OTurb', 'S'] = prop.h_p(streams.at['OEvap-OTurb', 'H'],streams.at['OEvap-OTurb', 'P'],streams.at['OEvap-OTurb', 'X'])["S"]
#
#         KPDturb_orc = KPDturb
#         #ORC Turb
#         streams.at['OTurb-OCool', 'P'] = streams.at['OCool-OPump', 'P']
#         streams.at['OTurb-OCool', 'X'] = streams.at['OCool-OPump', 'X']
#         streams.at['OTurb-OCool', 'G'] = streams.at['OCool-OPump', 'G']
#         Ht = prop.p_s(streams.at['OTurb-OCool', 'P'], streams.at['OEvap-OTurb', 'S'], streams.at['OTurb-OCool', 'X'])['H']
#         streams.at['OTurb-OCool', 'H'] = streams.at['OEvap-OTurb', 'H'] - (streams.at['OEvap-OTurb', 'H']-Ht)*KPDturb_orc
#         streams.at['OTurb-OCool', 'T'] = prop.h_p(streams.at['OTurb-OCool', 'H'],streams.at['OTurb-OCool', 'P'],streams.at['OTurb-OCool', 'X'])["T"]
#         streams.at['OTurb-OCool', 'Q'] = prop.h_p(streams.at['OTurb-OCool', 'H'],streams.at['OTurb-OCool', 'P'],streams.at['OTurb-OCool', 'X'])["Q"]
#         streams.at['OTurb-OCool', 'S'] = prop.h_p(streams.at['OTurb-OCool', 'H'],streams.at['OTurb-OCool', 'P'],streams.at['OTurb-OCool', 'X'])["S"]
#
#
#
#
#
#         Q = streams.at['Reg-Heat', 'G']*(streams.at['Heat-Turb', 'H']-streams.at['Reg-Heat', 'H'])
#         Nco2 = 0.99*streams.at['Heat-Turb', 'G']*(streams.at['Heat-Turb', 'H']-streams.at['Turb-Reg', 'H']) - streams.at['Cool-Comp', 'G']*(streams.at['Comp-Reg', 'H']-streams.at['Cool-Comp', 'H'])/0.99
#         Norc = 0.99*streams.at['OTurb-OCool', 'G']*(streams.at['OEvap-OTurb', 'H']-streams.at['OTurb-OCool', 'H']) - streams.at['OEvap-OTurb', 'G']*(streams.at['OPump-OEvap', 'H']-streams.at['OCool-OPump', 'H'])/0.99
#         KPD = (Nco2+Norc)/Q*100
#         # Разбиение и проверка темп. напора
#         n = 50
#         T_co2_evap = np.zeros(n)
#         H_co2_evap = np.linspace(streams.at['Reg-Evap', 'H'], streams.at['Evap-Cool', 'H'], n)
#         for i in range(0, n):
#             T_co2_evap[i] = prop.h_p(H_co2_evap[i], streams.at['Reg-Evap', 'P'], streams.at['Reg-Evap', 'X'])["T"]
#         Q_co2_evap = -(H_co2_evap - streams.at['Reg-Evap', 'H']) * streams.at['Cool-Comp', 'G']
#
#         T_orc_evap = np.zeros(n)
#         H_orc_evap = np.linspace(streams.at['OEvap-OTurb', 'H'], streams.at['OPump-OEvap', 'H'], n)
#         for i in range(0, n):
#             T_orc_evap[i] = prop.h_p(H_orc_evap[i], streams.at['OPump-OEvap', 'P'], streams.at['OPump-OEvap', 'X'])["T"]
#         Q_orc_evap = -(H_orc_evap - streams.at['OEvap-OTurb', 'H']) * streams.at['OEvap-OTurb', 'G']
#         deltaT = T_co2_evap - T_orc_evap
#         # print(streams)
#         # print(min(deltaT))
#         # print(Pk_orc)
#         if min(deltaT) < 10:
#             KPD = float(0)
#         if streams.at['OEvap-OTurb', 'Q'] != 1:
#             KPD = float(0)
#         print(round(pi, 3),
#               round(G_orc, 3),
#               round(H_percent, 3),
#               round(pi_orc, 3),
#               round(KPD, 3))
#         # plt.plot(Q_co2_evap,T_co2_evap)
#         # plt.plot(Q_orc_evap,T_orc_evap)
#         # plt.show()

