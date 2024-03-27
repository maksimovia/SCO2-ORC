import pandas as pd
import prop
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

streams = pd.read_excel('streams.xlsx', index_col=0, sheet_name='recomp')

def Calc(Input):
    H_percent = Input[0]
    dT_lt = Input[1]
    pi=4
    Tmax = 600
    Pk = 7.8
    KPDcomp = 0.9
    KPDturb = 0.9
    Tmin = 32
    x = 0.6
    Gas ='CO2'
    G0 = 1
    #dT_lt = 10
    #heat
    streams.at['Heat-Turb', 'T'] = Tmax
    streams.at['Heat-Turb', 'P'] = Pk*pi
    streams.at['Heat-Turb', 'X'] = Gas
    streams.at['Heat-Turb', 'H'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['H']
    streams.at['Heat-Turb', 'G'] = G0
    streams.at['Heat-Turb', 'Q'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['Q']
    streams.at['Heat-Turb', 'S'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['S']

    #turb
    streams.at['Turb-HTH', 'P'] = Pk
    streams.at['Turb-HTH', 'X'] = Gas
    streams.at['Turb-HTH', 'G'] = G0
    Ht = prop.p_s(streams.at['Turb-HTH', 'P'], streams.at['Heat-Turb', 'S'], streams.at['Turb-HTH', 'X'])['H']
    streams.at['Turb-HTH', 'H'] = streams.at['Heat-Turb', 'H'] - (streams.at['Heat-Turb', 'H']-Ht)*KPDturb
    streams.at['Turb-HTH', 'T'] = prop.h_p(streams.at['Turb-HTH', 'H'],streams.at['Turb-HTH', 'P'],streams.at['Turb-HTH', 'X'])["T"]
    streams.at['Turb-HTH', 'Q'] = prop.h_p(streams.at['Turb-HTH', 'H'],streams.at['Turb-HTH', 'P'],streams.at['Turb-HTH', 'X'])["Q"]
    streams.at['Turb-HTH', 'S'] = prop.h_p(streams.at['Turb-HTH', 'H'],streams.at['Turb-HTH', 'P'],streams.at['Turb-HTH', 'X'])["S"]

    #ht hot
    streams.at['HTH-LTH', 'P'] = Pk
    streams.at['HTH-LTH', 'X'] = Gas
    streams.at['HTH-LTH', 'G'] = G0


    #lt hot
    streams.at['LTH-SPLT', 'P'] = Pk
    streams.at['LTH-SPLT', 'X'] = Gas
    streams.at['LTH-SPLT', 'G'] = G0

    #s1
    streams.at['SPLT-Cool', 'P'] = Pk
    streams.at['SPLT-Cool', 'X'] = Gas
    streams.at['SPLT-Cool', 'G'] = G0*x

    #Cool
    streams.at['Cool-MC', 'P'] = Pk
    streams.at['Cool-MC', 'X'] = Gas
    streams.at['Cool-MC', 'G'] = streams.at['SPLT-Cool', 'G']
    streams.at['Cool-MC', 'T'] = Tmin
    streams.at['Cool-MC', 'H'] = prop.t_p(streams.at['Cool-MC', 'T'],streams.at['Cool-MC', 'P'],Gas)['H']
    streams.at['Cool-MC', 'Q'] = prop.t_p(streams.at['Cool-MC', 'T'],streams.at['Cool-MC', 'P'],Gas)['Q']
    streams.at['Cool-MC', 'S'] = prop.t_p(streams.at['Cool-MC', 'T'],streams.at['Cool-MC', 'P'],Gas)['S']

    #MC
    streams.at['MC-LTH', 'P'] = streams.at['SPLT-Cool', 'P']*pi
    streams.at['MC-LTH', 'X'] = Gas
    streams.at['MC-LTH', 'G'] = streams.at['Cool-MC', 'G']
    Ht = prop.p_s(streams.at['MC-LTH', 'P'], streams.at['Cool-MC', 'S'], Gas)['H']
    streams.at['MC-LTH', 'H'] = streams.at['Cool-MC', 'H'] + (Ht - streams.at['Cool-MC', 'H'])/KPDcomp
    streams.at['MC-LTH', 'T'] = prop.h_p(streams.at['MC-LTH', 'H'],streams.at['MC-LTH', 'P'],streams.at['MC-LTH', 'X'])["T"]
    streams.at['MC-LTH', 'Q'] = prop.h_p(streams.at['MC-LTH', 'H'],streams.at['MC-LTH', 'P'],streams.at['MC-LTH', 'X'])["Q"]
    streams.at['MC-LTH', 'S'] = prop.h_p(streams.at['MC-LTH', 'H'],streams.at['MC-LTH', 'P'],streams.at['MC-LTH', 'X'])["S"]

    #lt hot
    streams.at['LTH-SPLT', 'T'] = streams.at['MC-LTH', 'T'] + dT_lt
    streams.at['LTH-SPLT', 'H'] = prop.t_p(streams.at['LTH-SPLT', 'T'],streams.at['LTH-SPLT', 'P'],Gas)['H']
    streams.at['LTH-SPLT', 'Q'] = prop.t_p(streams.at['LTH-SPLT', 'T'],streams.at['LTH-SPLT', 'P'],Gas)['Q']
    streams.at['LTH-SPLT', 'S'] = prop.t_p(streams.at['LTH-SPLT', 'T'],streams.at['LTH-SPLT', 'P'],Gas)['S']

    #H_percent = 0.33
    #ht hot
    streams.at['HTH-LTH', 'H'] = (streams.at['Turb-HTH', 'H']-streams.at['LTH-SPLT', 'H'])*H_percent + streams.at['LTH-SPLT', 'H']
    streams.at['HTH-LTH', 'T'] = prop.h_p(streams.at['HTH-LTH', 'H'],streams.at['HTH-LTH', 'P'],Gas)['T']
    streams.at['HTH-LTH', 'Q'] = prop.h_p(streams.at['HTH-LTH', 'H'],streams.at['HTH-LTH', 'P'],Gas)['Q']
    streams.at['HTH-LTH', 'S'] = prop.h_p(streams.at['HTH-LTH', 'H'],streams.at['HTH-LTH', 'P'],Gas)['S']

    #lt cold
    streams.at['LTH-MIX', 'P'] = streams.at['MC-LTH', 'P']
    streams.at['LTH-MIX', 'X'] = Gas
    streams.at['LTH-MIX', 'G'] = streams.at['Cool-MC', 'G']
    streams.at['LTH-MIX', 'H'] = streams.at['MC-LTH', 'H'] + (G0/streams.at['LTH-MIX', 'G'])*(streams.at['HTH-LTH', 'H']-streams.at['LTH-SPLT', 'H'])
    streams.at['LTH-MIX', 'T'] = prop.h_p(streams.at['LTH-MIX', 'H'],streams.at['LTH-MIX', 'P'],Gas)['T']
    streams.at['LTH-MIX', 'Q'] = prop.h_p(streams.at['LTH-MIX', 'H'],streams.at['LTH-MIX', 'P'],Gas)['Q']
    streams.at['LTH-MIX', 'S'] = prop.h_p(streams.at['LTH-MIX', 'H'],streams.at['LTH-MIX', 'P'],Gas)['S']

    #s2
    streams.at['SPLT-RC', 'P'] = Pk
    streams.at['SPLT-RC', 'X'] = Gas
    streams.at['SPLT-RC', 'G'] = G0*(1-x)

    #RC
    streams.at['RC-MIX', 'P'] = streams.at['SPLT-RC', 'P']*pi
    streams.at['RC-MIX', 'X'] = Gas
    streams.at['RC-MIX', 'G'] = streams.at['SPLT-RC', 'G']
    Ht = prop.p_s(streams.at['RC-MIX', 'P'], streams.at['LTH-SPLT', 'S'], Gas)['H']
    streams.at['RC-MIX', 'H'] = streams.at['LTH-SPLT', 'H'] + (Ht - streams.at['LTH-SPLT', 'H'])/KPDcomp
    streams.at['RC-MIX', 'T'] = prop.h_p(streams.at['RC-MIX', 'H'],streams.at['RC-MIX', 'P'],streams.at['RC-MIX', 'X'])["T"]
    streams.at['RC-MIX', 'Q'] = prop.h_p(streams.at['RC-MIX', 'H'],streams.at['RC-MIX', 'P'],streams.at['RC-MIX', 'X'])["Q"]
    streams.at['RC-MIX', 'S'] = prop.h_p(streams.at['RC-MIX', 'H'],streams.at['RC-MIX', 'P'],streams.at['RC-MIX', 'X'])["S"]

    #MIX
    streams.at['MIX-HTH', 'P'] = Pk*pi
    streams.at['MIX-HTH', 'X'] = Gas
    streams.at['MIX-HTH', 'G'] = streams.at['RC-MIX', 'G'] + streams.at['MC-LTH', 'G']
    streams.at['MIX-HTH', 'H'] = (streams.at['RC-MIX', 'G']*streams.at['RC-MIX', 'H']+streams.at['LTH-MIX', 'G']*streams.at['LTH-MIX', 'H'])/streams.at['MIX-HTH', 'G']
    streams.at['MIX-HTH', 'T'] = prop.h_p(streams.at['MIX-HTH', 'H'],streams.at['MIX-HTH', 'P'],Gas)["T"]
    streams.at['MIX-HTH', 'Q'] = prop.h_p(streams.at['MIX-HTH', 'H'],streams.at['MIX-HTH', 'P'],Gas)["Q"]
    streams.at['MIX-HTH', 'S'] = prop.h_p(streams.at['MIX-HTH', 'H'],streams.at['MIX-HTH', 'P'],Gas)["S"]

    #ht cold
    streams.at['HTH-Heat', 'P'] = Pk*pi
    streams.at['HTH-Heat', 'X'] = Gas
    streams.at['HTH-Heat', 'G'] = streams.at['MIX-HTH', 'G']
    streams.at['HTH-Heat', 'H'] = streams.at['Turb-HTH', 'H'] - streams.at['HTH-LTH', 'H'] + streams.at['MIX-HTH', 'H']
    streams.at['HTH-Heat', 'T'] = prop.h_p(streams.at['HTH-Heat', 'H'],streams.at['HTH-Heat', 'P'],Gas)['T']
    streams.at['HTH-Heat', 'Q'] = prop.h_p(streams.at['HTH-Heat', 'H'],streams.at['HTH-Heat', 'P'],Gas)['Q']
    streams.at['HTH-Heat', 'S'] = prop.h_p(streams.at['HTH-Heat', 'H'],streams.at['HTH-Heat', 'P'],Gas)['S']

    # Разбиение и проверка темп. напора
    n = 10
    T_co2_ht_hot = np.zeros(n)
    H_co2_ht_hot = np.linspace(streams.at['Turb-HTH', 'H'], streams.at['HTH-LTH', 'H'], n)
    for i in range(0, n):
        T_co2_ht_hot[i] = prop.h_p(H_co2_ht_hot[i], streams.at['HTH-LTH', 'P'], Gas)["T"]
    Q_co2_ht_hot = (H_co2_ht_hot - streams.at['HTH-LTH', 'H']) * streams.at['HTH-LTH', 'G']

    T_co2_ht_cold = np.zeros(n)
    H_co2_ht_cold = np.linspace(streams.at['HTH-Heat', 'H'], streams.at['MIX-HTH', 'H'], n)
    for i in range(0, n):
        T_co2_ht_cold[i] = prop.h_p(H_co2_ht_cold[i], streams.at['MIX-HTH', 'P'], Gas)["T"]
    Q_co2_ht_cold = (H_co2_ht_cold - streams.at['MIX-HTH', 'H']) * streams.at['MIX-HTH', 'G']
    # plt.plot(Q_co2_ht_hot,T_co2_ht_hot)
    # plt.plot(Q_co2_ht_cold,T_co2_ht_cold)
    # plt.show()

    n = 10
    T_co2_lt_hot = np.zeros(n)
    H_co2_lt_hot = np.linspace(streams.at['HTH-LTH', 'H'], streams.at['LTH-SPLT', 'H'], n)
    for i in range(0, n):
        T_co2_lt_hot[i] = prop.h_p(H_co2_lt_hot[i], streams.at['HTH-LTH', 'P'], Gas)["T"]
    Q_co2_lt_hot = (H_co2_lt_hot - streams.at['LTH-SPLT', 'H']) * streams.at['HTH-LTH', 'G']

    T_co2_lt_cold = np.zeros(n)
    H_co2_lt_cold = np.linspace(streams.at['LTH-MIX', 'H'], streams.at['MC-LTH', 'H'], n)
    for i in range(0, n):
        T_co2_lt_cold[i] = prop.h_p(H_co2_lt_cold[i], streams.at['MC-LTH', 'P'], Gas)["T"]
    Q_co2_lt_cold = (H_co2_lt_cold - streams.at['MC-LTH', 'H']) * streams.at['MC-LTH', 'G']
    # plt.plot(Q_co2_lt_hot,T_co2_lt_hot)
    # plt.plot(Q_co2_lt_cold,T_co2_lt_cold)
    # plt.show()

    #KPD
    NCO2_MC = streams.at['Cool-MC', 'G']*(streams.at['MC-LTH', 'H']-streams.at['Cool-MC', 'H'])
    NCO2_RC = streams.at['RC-MIX', 'G'] * (streams.at['RC-MIX', 'H'] - streams.at['LTH-SPLT', 'H'])
    NCO2_Turb = G0*(streams.at['Heat-Turb', 'H']-streams.at['Turb-HTH', 'H'])
    NCO2_Heat = G0*(streams.at['Heat-Turb', 'H']-streams.at['HTH-Heat', 'H'])
    KPD = 0.99*0.99*(NCO2_Turb-(NCO2_MC+NCO2_RC))/NCO2_Heat*100

    # print(streams)
    print(min(T_co2_lt_hot-T_co2_lt_cold))
    print(min(T_co2_ht_hot-T_co2_ht_cold))
    print(KPD)
    return min(T_co2_ht_hot-T_co2_ht_cold)-10, min(T_co2_lt_hot-T_co2_lt_cold)-10
res = opt.root(Calc,(0.3,10),method='hybr',tol=0.001)
Calc(res.x)
