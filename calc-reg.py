import pandas as pd
import prop

streams = pd.read_excel('streams.xlsx', index_col=0, sheet_name='reg')

Tmax = 600
Pk = 7.8
pi=2
KPDcomp = 0.9
KPDturb = 0.9
Tmin = 32
dTreg = 50

#cool
streams.at['Cool-Comp', 'T'] = Tmin
streams.at['Cool-Comp', 'P'] = Pk
streams.at['Cool-Comp', 'X'] = "CO2"
streams.at['Cool-Comp', 'H'] = prop.t_p(streams.at['Cool-Comp', 'T'],streams.at['Cool-Comp', 'P'],streams.at['Cool-Comp', 'X'])['H']
streams.at['Cool-Comp', 'G'] = 1
streams.at['Cool-Comp', 'Q'] = 1
streams.at['Cool-Comp', 'S'] = prop.h_p(streams.at['Cool-Comp', 'H'],streams.at['Cool-Comp', 'P'],streams.at['Cool-Comp', 'X'])['S']

#comp
streams.at['Comp-Reg', 'P'] = streams.at['Cool-Comp', 'P']*pi
streams.at['Comp-Reg', 'X'] = streams.at['Cool-Comp', 'X']
Ht = prop.p_s(streams.at['Comp-Reg', 'P'], streams.at['Cool-Comp', 'S'], streams.at['Comp-Reg', 'X'])['H']
streams.at['Comp-Reg', 'H'] = streams.at['Cool-Comp', 'H'] + (Ht - streams.at['Cool-Comp', 'H'])/KPDcomp
streams.at['Comp-Reg', 'G'] = streams.at['Cool-Comp', 'G']
streams.at['Comp-Reg', 'T'] = prop.h_p(streams.at['Comp-Reg', 'H'],streams.at['Comp-Reg', 'P'],streams.at['Comp-Reg', 'X'])["T"]
streams.at['Comp-Reg', 'Q'] = prop.h_p(streams.at['Comp-Reg', 'H'],streams.at['Comp-Reg', 'P'],streams.at['Comp-Reg', 'X'])["Q"]
streams.at['Comp-Reg', 'S'] = prop.h_p(streams.at['Comp-Reg', 'H'],streams.at['Comp-Reg', 'P'],streams.at['Comp-Reg', 'X'])["S"]

#reg_cold
streams.at['Reg-Heat', 'P'] = streams.at['Comp-Reg', 'P']
streams.at['Reg-Heat', 'X'] = streams.at['Comp-Reg', 'X']
streams.at['Reg-Heat', 'G'] = streams.at['Comp-Reg', 'G']


# #heat
streams.at['Heat-Turb', 'T'] = Tmax
streams.at['Heat-Turb', 'P'] = streams.at['Reg-Heat', 'P']
streams.at['Heat-Turb', 'X'] = streams.at['Reg-Heat', 'X']
streams.at['Heat-Turb', 'G'] = streams.at['Reg-Heat', 'G']
streams.at['Heat-Turb', 'H'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['H']
streams.at['Heat-Turb', 'Q'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['Q']
streams.at['Heat-Turb', 'S'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['S']

#turb
streams.at['Turb-Reg', 'P'] = Pk
streams.at['Turb-Reg', 'X'] = streams.at['Heat-Turb', 'X']
streams.at['Turb-Reg', 'G'] = streams.at['Heat-Turb', 'G']
Ht = prop.p_s(streams.at['Turb-Reg', 'P'], streams.at['Heat-Turb', 'S'], streams.at['Turb-Reg', 'X'])['H']
streams.at['Turb-Reg', 'H'] = streams.at['Heat-Turb', 'H'] - (streams.at['Heat-Turb', 'H']-Ht)*KPDturb
streams.at['Turb-Reg', 'T'] = prop.h_p(streams.at['Turb-Reg', 'H'],streams.at['Turb-Reg', 'P'],streams.at['Turb-Reg', 'X'])["T"]
streams.at['Turb-Reg', 'Q'] = prop.h_p(streams.at['Turb-Reg', 'H'],streams.at['Turb-Reg', 'P'],streams.at['Turb-Reg', 'X'])["Q"]
streams.at['Turb-Reg', 'S'] = prop.h_p(streams.at['Turb-Reg', 'H'],streams.at['Turb-Reg', 'P'],streams.at['Turb-Reg', 'X'])["S"]

#reg_hot
streams.at['Reg-Cool', 'T'] = streams.at['Comp-Reg', 'T']+dTreg
streams.at['Reg-Cool', 'P'] = streams.at['Turb-Reg', 'P']
streams.at['Reg-Cool', 'X'] = streams.at['Turb-Reg', 'X']
streams.at['Reg-Cool', 'G'] = streams.at['Turb-Reg', 'G']
streams.at['Reg-Cool', 'H'] = prop.t_p(streams.at['Reg-Cool', 'T'],streams.at['Reg-Cool', 'P'],streams.at['Reg-Cool', 'X'])['H']
streams.at['Reg-Cool', 'Q'] = prop.t_p(streams.at['Reg-Cool', 'T'],streams.at['Reg-Cool', 'P'],streams.at['Reg-Cool', 'X'])['Q']
streams.at['Reg-Cool', 'S'] = prop.t_p(streams.at['Reg-Cool', 'T'],streams.at['Reg-Cool', 'P'],streams.at['Reg-Cool', 'X'])['S']

#reg_cold
streams.at['Reg-Heat', 'H'] = streams.at['Comp-Reg', 'H']+(streams.at['Turb-Reg', 'H'] - streams.at['Reg-Cool', 'H'])
streams.at['Reg-Heat', 'T'] = prop.h_p(streams.at['Reg-Heat', 'H'],streams.at['Reg-Heat', 'P'],streams.at['Reg-Heat', 'X'])['T']
streams.at['Reg-Heat', 'Q'] = prop.h_p(streams.at['Reg-Heat', 'H'],streams.at['Reg-Heat', 'P'],streams.at['Reg-Heat', 'X'])['Q']
streams.at['Reg-Heat', 'S'] = prop.h_p(streams.at['Reg-Heat', 'H'],streams.at['Reg-Heat', 'P'],streams.at['Reg-Heat', 'X'])['S']




Q = streams.at['Reg-Heat', 'G']*(streams.at['Heat-Turb', 'H']-streams.at['Reg-Heat', 'H'])
Nco2 = 0.99*streams.at['Heat-Turb', 'G']*(streams.at['Heat-Turb', 'H']-streams.at['Turb-Reg', 'H']) - streams.at['Cool-Comp', 'G']*(streams.at['Comp-Reg', 'H']-streams.at['Cool-Comp', 'H'])/0.99
KPD = Nco2/Q*100


print(streams)
print(KPD)

