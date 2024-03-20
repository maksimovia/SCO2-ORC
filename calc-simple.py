import pandas as pd
import prop

streams = pd.read_excel('streams.xlsx', index_col=0, sheet_name='simple')

Tmax = 600
Pk = 7.8
pi = 2
KPDcomp = 0.9
KPDturb = 0.9
Tmin = 32

#cool
streams.at['Cool-Comp', 'T'] = Tmin
streams.at['Cool-Comp', 'P'] = Pk
streams.at['Cool-Comp', 'X'] = "CO2"
streams.at['Cool-Comp', 'H'] = prop.t_p(streams.at['Cool-Comp', 'T'],streams.at['Cool-Comp', 'P'],streams.at['Cool-Comp', 'X'])['H']
streams.at['Cool-Comp', 'G'] = 1
streams.at['Cool-Comp', 'Q'] = 1
streams.at['Cool-Comp', 'S'] = prop.h_p(streams.at['Cool-Comp', 'H'],streams.at['Cool-Comp', 'P'],streams.at['Cool-Comp', 'X'])['S']

#comp
streams.at['Comp-Heat', 'P'] = streams.at['Cool-Comp', 'P']*pi
streams.at['Comp-Heat', 'X'] = streams.at['Cool-Comp', 'X']
Ht = prop.p_s(streams.at['Comp-Heat', 'P'], streams.at['Cool-Comp', 'S'], streams.at['Comp-Heat', 'X'])['H']
streams.at['Comp-Heat', 'H'] = streams.at['Cool-Comp', 'H'] + (Ht - streams.at['Cool-Comp', 'H'])/KPDcomp
streams.at['Comp-Heat', 'G'] = streams.at['Cool-Comp', 'G']
streams.at['Comp-Heat', 'T'] = prop.h_p(streams.at['Comp-Heat', 'H'],streams.at['Comp-Heat', 'P'],streams.at['Comp-Heat', 'X'])["T"]
streams.at['Comp-Heat', 'Q'] = prop.h_p(streams.at['Comp-Heat', 'H'],streams.at['Comp-Heat', 'P'],streams.at['Comp-Heat', 'X'])["Q"]
streams.at['Comp-Heat', 'S'] = prop.h_p(streams.at['Comp-Heat', 'H'],streams.at['Comp-Heat', 'P'],streams.at['Comp-Heat', 'X'])["S"]

#heat
streams.at['Heat-Turb', 'T'] = Tmax
streams.at['Heat-Turb', 'P'] = streams.at['Comp-Heat', 'P']
streams.at['Heat-Turb', 'X'] = streams.at['Comp-Heat', 'X']
streams.at['Heat-Turb', 'H'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['H']
streams.at['Heat-Turb', 'G'] = streams.at['Comp-Heat', 'G']
streams.at['Heat-Turb', 'Q'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['Q']
streams.at['Heat-Turb', 'S'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['S']

#turb
streams.at['Turb-Cool', 'P'] = Pk
streams.at['Turb-Cool', 'X'] = streams.at['Heat-Turb', 'X']
streams.at['Turb-Cool', 'G'] = streams.at['Heat-Turb', 'G']
Ht = prop.p_s(streams.at['Turb-Cool', 'P'], streams.at['Heat-Turb', 'S'], streams.at['Turb-Cool', 'X'])['H']
streams.at['Turb-Cool', 'H'] = streams.at['Heat-Turb', 'H'] - (streams.at['Heat-Turb', 'H']-Ht)*KPDturb
streams.at['Turb-Cool', 'T'] = prop.h_p(streams.at['Turb-Cool', 'H'],streams.at['Turb-Cool', 'P'],streams.at['Turb-Cool', 'X'])["T"]
streams.at['Turb-Cool', 'Q'] = prop.h_p(streams.at['Turb-Cool', 'H'],streams.at['Turb-Cool', 'P'],streams.at['Turb-Cool', 'X'])["Q"]
streams.at['Turb-Cool', 'S'] = prop.h_p(streams.at['Turb-Cool', 'H'],streams.at['Turb-Cool', 'P'],streams.at['Turb-Cool', 'X'])["S"]




print(streams)