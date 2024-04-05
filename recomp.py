import pandas as pd
import prop
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

streams = pd.read_excel('streams.xlsx', index_col=0, sheet_name='recomp')
def Calc(pi,x,Pk,Tmax):
    def CalcRecomp(Input):
        H_percent = Input[0]
        dT_lt = Input[1]
        #pi=4
        #Tmax = 500
        #Pk = 7.8
        KPDcomp = 0.9
        KPDturb = 0.9
        Tmin = 32
        #x = 0.6
        Gas ='CO2'
        G0 = 1
        #dT_lt = 10
        dPpot = 0.05
        #heat
        streams.at['Heat-Turb', 'T'] = Tmax
        streams.at['Heat-Turb', 'P'] = Pk*pi - 3*dPpot
        streams.at['Heat-Turb', 'X'] = Gas
        streams.at['Heat-Turb', 'H'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['H']
        streams.at['Heat-Turb', 'G'] = G0
        streams.at['Heat-Turb', 'Q'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['Q']
        streams.at['Heat-Turb', 'S'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['S']

        #turb
        streams.at['Turb-HTH', 'P'] = Pk + 3*dPpot
        streams.at['Turb-HTH', 'X'] = Gas
        streams.at['Turb-HTH', 'G'] = G0
        Ht = prop.p_s(streams.at['Turb-HTH', 'P'], streams.at['Heat-Turb', 'S'], streams.at['Turb-HTH', 'X'])['H']
        streams.at['Turb-HTH', 'H'] = streams.at['Heat-Turb', 'H'] - (streams.at['Heat-Turb', 'H']-Ht)*KPDturb
        streams.at['Turb-HTH', 'T'] = prop.h_p(streams.at['Turb-HTH', 'H'],streams.at['Turb-HTH', 'P'],streams.at['Turb-HTH', 'X'])["T"]
        streams.at['Turb-HTH', 'Q'] = prop.h_p(streams.at['Turb-HTH', 'H'],streams.at['Turb-HTH', 'P'],streams.at['Turb-HTH', 'X'])["Q"]
        streams.at['Turb-HTH', 'S'] = prop.h_p(streams.at['Turb-HTH', 'H'],streams.at['Turb-HTH', 'P'],streams.at['Turb-HTH', 'X'])["S"]

        #ht hot
        streams.at['HTH-LTH', 'P'] = Pk + 2*dPpot
        streams.at['HTH-LTH', 'X'] = Gas
        streams.at['HTH-LTH', 'G'] = G0


        #lt hot
        streams.at['LTH-SPLT', 'P'] = Pk - dPpot
        streams.at['LTH-SPLT', 'X'] = Gas
        streams.at['LTH-SPLT', 'G'] = G0

        #s1
        streams.at['SPLT-Cool', 'P'] = Pk - dPpot
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
        streams.at['LTH-MIX', 'P'] = streams.at['MC-LTH', 'P'] - dPpot
        streams.at['LTH-MIX', 'X'] = Gas
        streams.at['LTH-MIX', 'G'] = streams.at['Cool-MC', 'G']
        streams.at['LTH-MIX', 'H'] = streams.at['MC-LTH', 'H'] + (G0/streams.at['LTH-MIX', 'G'])*(streams.at['HTH-LTH', 'H']-streams.at['LTH-SPLT', 'H'])
        streams.at['LTH-MIX', 'T'] = prop.h_p(streams.at['LTH-MIX', 'H'],streams.at['LTH-MIX', 'P'],Gas)['T']
        streams.at['LTH-MIX', 'Q'] = prop.h_p(streams.at['LTH-MIX', 'H'],streams.at['LTH-MIX', 'P'],Gas)['Q']
        streams.at['LTH-MIX', 'S'] = prop.h_p(streams.at['LTH-MIX', 'H'],streams.at['LTH-MIX', 'P'],Gas)['S']

        #s2
        streams.at['SPLT-RC', 'P'] = Pk + dPpot
        streams.at['SPLT-RC', 'X'] = Gas
        streams.at['SPLT-RC', 'G'] = G0*(1-x)

        #RC
        streams.at['RC-MIX', 'P'] = Pk*pi - dPpot
        streams.at['RC-MIX', 'X'] = Gas
        streams.at['RC-MIX', 'G'] = streams.at['SPLT-RC', 'G']
        Ht = prop.p_s(streams.at['RC-MIX', 'P'], streams.at['LTH-SPLT', 'S'], Gas)['H']
        streams.at['RC-MIX', 'H'] = streams.at['LTH-SPLT', 'H'] + (Ht - streams.at['LTH-SPLT', 'H'])/KPDcomp
        streams.at['RC-MIX', 'T'] = prop.h_p(streams.at['RC-MIX', 'H'],streams.at['RC-MIX', 'P'],streams.at['RC-MIX', 'X'])["T"]
        streams.at['RC-MIX', 'Q'] = prop.h_p(streams.at['RC-MIX', 'H'],streams.at['RC-MIX', 'P'],streams.at['RC-MIX', 'X'])["Q"]
        streams.at['RC-MIX', 'S'] = prop.h_p(streams.at['RC-MIX', 'H'],streams.at['RC-MIX', 'P'],streams.at['RC-MIX', 'X'])["S"]

        #MIX
        streams.at['MIX-HTH', 'P'] = Pk*pi - dPpot
        streams.at['MIX-HTH', 'X'] = Gas
        streams.at['MIX-HTH', 'G'] = streams.at['RC-MIX', 'G'] + streams.at['MC-LTH', 'G']
        streams.at['MIX-HTH', 'H'] = (streams.at['RC-MIX', 'G']*streams.at['RC-MIX', 'H']+streams.at['LTH-MIX', 'G']*streams.at['LTH-MIX', 'H'])/streams.at['MIX-HTH', 'G']
        streams.at['MIX-HTH', 'T'] = prop.h_p(streams.at['MIX-HTH', 'H'],streams.at['MIX-HTH', 'P'],Gas)["T"]
        streams.at['MIX-HTH', 'Q'] = prop.h_p(streams.at['MIX-HTH', 'H'],streams.at['MIX-HTH', 'P'],Gas)["Q"]
        streams.at['MIX-HTH', 'S'] = prop.h_p(streams.at['MIX-HTH', 'H'],streams.at['MIX-HTH', 'P'],Gas)["S"]

        #ht cold
        streams.at['HTH-Heat', 'P'] = Pk*pi - 2*dPpot
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
        # global KPD
        # global T
        T = streams.at['LTH-SPLT', 'T']
        KPD = 0.99*0.99*(NCO2_Turb-(NCO2_MC+NCO2_RC))/NCO2_Heat*100

        # print(streams)
        # print(min(T_co2_lt_hot-T_co2_lt_cold))
        # print(min(T_co2_ht_hot-T_co2_ht_cold))
        # print(KPD)







        return min(T_co2_ht_hot-T_co2_ht_cold)-10, min(T_co2_lt_hot-T_co2_lt_cold)-10
    res = opt.root(CalcRecomp,(0.3,10),method='hybr',tol=0.001)
    H_percent = res.x[0]
    dT_lt = res.x[1]
    # pi=4
    # Tmax = 500
    # Pk = 7.8
    KPDcomp = 0.9
    KPDturb = 0.9
    Tmin = 32
    # x = 0.6
    Gas = 'CO2'
    G0 = 1
    # dT_lt = 10
    dPpot = 0.05
    # heat
    streams.at['Heat-Turb', 'T'] = Tmax
    streams.at['Heat-Turb', 'P'] = Pk * pi - 3 * dPpot
    streams.at['Heat-Turb', 'X'] = Gas
    streams.at['Heat-Turb', 'H'] = \
    prop.t_p(streams.at['Heat-Turb', 'T'], streams.at['Heat-Turb', 'P'], streams.at['Heat-Turb', 'X'])['H']
    streams.at['Heat-Turb', 'G'] = G0
    streams.at['Heat-Turb', 'Q'] = \
    prop.t_p(streams.at['Heat-Turb', 'T'], streams.at['Heat-Turb', 'P'], streams.at['Heat-Turb', 'X'])['Q']
    streams.at['Heat-Turb', 'S'] = \
    prop.t_p(streams.at['Heat-Turb', 'T'], streams.at['Heat-Turb', 'P'], streams.at['Heat-Turb', 'X'])['S']

    # turb
    streams.at['Turb-HTH', 'P'] = Pk + 3 * dPpot
    streams.at['Turb-HTH', 'X'] = Gas
    streams.at['Turb-HTH', 'G'] = G0
    Ht = prop.p_s(streams.at['Turb-HTH', 'P'], streams.at['Heat-Turb', 'S'], streams.at['Turb-HTH', 'X'])['H']
    streams.at['Turb-HTH', 'H'] = streams.at['Heat-Turb', 'H'] - (streams.at['Heat-Turb', 'H'] - Ht) * KPDturb
    streams.at['Turb-HTH', 'T'] = \
    prop.h_p(streams.at['Turb-HTH', 'H'], streams.at['Turb-HTH', 'P'], streams.at['Turb-HTH', 'X'])["T"]
    streams.at['Turb-HTH', 'Q'] = \
    prop.h_p(streams.at['Turb-HTH', 'H'], streams.at['Turb-HTH', 'P'], streams.at['Turb-HTH', 'X'])["Q"]
    streams.at['Turb-HTH', 'S'] = \
    prop.h_p(streams.at['Turb-HTH', 'H'], streams.at['Turb-HTH', 'P'], streams.at['Turb-HTH', 'X'])["S"]

    # ht hot
    streams.at['HTH-LTH', 'P'] = Pk + 2 * dPpot
    streams.at['HTH-LTH', 'X'] = Gas
    streams.at['HTH-LTH', 'G'] = G0

    # lt hot
    streams.at['LTH-SPLT', 'P'] = Pk - dPpot
    streams.at['LTH-SPLT', 'X'] = Gas
    streams.at['LTH-SPLT', 'G'] = G0

    # s1
    streams.at['SPLT-Cool', 'P'] = Pk - dPpot
    streams.at['SPLT-Cool', 'X'] = Gas
    streams.at['SPLT-Cool', 'G'] = G0 * x

    # Cool
    streams.at['Cool-MC', 'P'] = Pk
    streams.at['Cool-MC', 'X'] = Gas
    streams.at['Cool-MC', 'G'] = streams.at['SPLT-Cool', 'G']
    streams.at['Cool-MC', 'T'] = Tmin
    streams.at['Cool-MC', 'H'] = prop.t_p(streams.at['Cool-MC', 'T'], streams.at['Cool-MC', 'P'], Gas)['H']
    streams.at['Cool-MC', 'Q'] = prop.t_p(streams.at['Cool-MC', 'T'], streams.at['Cool-MC', 'P'], Gas)['Q']
    streams.at['Cool-MC', 'S'] = prop.t_p(streams.at['Cool-MC', 'T'], streams.at['Cool-MC', 'P'], Gas)['S']

    # MC
    streams.at['MC-LTH', 'P'] = streams.at['SPLT-Cool', 'P'] * pi
    streams.at['MC-LTH', 'X'] = Gas
    streams.at['MC-LTH', 'G'] = streams.at['Cool-MC', 'G']
    Ht = prop.p_s(streams.at['MC-LTH', 'P'], streams.at['Cool-MC', 'S'], Gas)['H']
    streams.at['MC-LTH', 'H'] = streams.at['Cool-MC', 'H'] + (Ht - streams.at['Cool-MC', 'H']) / KPDcomp
    streams.at['MC-LTH', 'T'] = \
    prop.h_p(streams.at['MC-LTH', 'H'], streams.at['MC-LTH', 'P'], streams.at['MC-LTH', 'X'])["T"]
    streams.at['MC-LTH', 'Q'] = \
    prop.h_p(streams.at['MC-LTH', 'H'], streams.at['MC-LTH', 'P'], streams.at['MC-LTH', 'X'])["Q"]
    streams.at['MC-LTH', 'S'] = \
    prop.h_p(streams.at['MC-LTH', 'H'], streams.at['MC-LTH', 'P'], streams.at['MC-LTH', 'X'])["S"]

    # lt hot
    streams.at['LTH-SPLT', 'T'] = streams.at['MC-LTH', 'T'] + dT_lt
    streams.at['LTH-SPLT', 'H'] = prop.t_p(streams.at['LTH-SPLT', 'T'], streams.at['LTH-SPLT', 'P'], Gas)['H']
    streams.at['LTH-SPLT', 'Q'] = prop.t_p(streams.at['LTH-SPLT', 'T'], streams.at['LTH-SPLT', 'P'], Gas)['Q']
    streams.at['LTH-SPLT', 'S'] = prop.t_p(streams.at['LTH-SPLT', 'T'], streams.at['LTH-SPLT', 'P'], Gas)['S']

    # H_percent = 0.33
    # ht hot
    streams.at['HTH-LTH', 'H'] = (streams.at['Turb-HTH', 'H'] - streams.at['LTH-SPLT', 'H']) * H_percent + streams.at[
        'LTH-SPLT', 'H']
    streams.at['HTH-LTH', 'T'] = prop.h_p(streams.at['HTH-LTH', 'H'], streams.at['HTH-LTH', 'P'], Gas)['T']
    streams.at['HTH-LTH', 'Q'] = prop.h_p(streams.at['HTH-LTH', 'H'], streams.at['HTH-LTH', 'P'], Gas)['Q']
    streams.at['HTH-LTH', 'S'] = prop.h_p(streams.at['HTH-LTH', 'H'], streams.at['HTH-LTH', 'P'], Gas)['S']

    # lt cold
    streams.at['LTH-MIX', 'P'] = streams.at['MC-LTH', 'P'] - dPpot
    streams.at['LTH-MIX', 'X'] = Gas
    streams.at['LTH-MIX', 'G'] = streams.at['Cool-MC', 'G']
    streams.at['LTH-MIX', 'H'] = streams.at['MC-LTH', 'H'] + (G0 / streams.at['LTH-MIX', 'G']) * (
                streams.at['HTH-LTH', 'H'] - streams.at['LTH-SPLT', 'H'])
    streams.at['LTH-MIX', 'T'] = prop.h_p(streams.at['LTH-MIX', 'H'], streams.at['LTH-MIX', 'P'], Gas)['T']
    streams.at['LTH-MIX', 'Q'] = prop.h_p(streams.at['LTH-MIX', 'H'], streams.at['LTH-MIX', 'P'], Gas)['Q']
    streams.at['LTH-MIX', 'S'] = prop.h_p(streams.at['LTH-MIX', 'H'], streams.at['LTH-MIX', 'P'], Gas)['S']

    # s2
    streams.at['SPLT-RC', 'P'] = Pk + dPpot
    streams.at['SPLT-RC', 'X'] = Gas
    streams.at['SPLT-RC', 'G'] = G0 * (1 - x)

    # RC
    streams.at['RC-MIX', 'P'] = Pk * pi - dPpot
    streams.at['RC-MIX', 'X'] = Gas
    streams.at['RC-MIX', 'G'] = streams.at['SPLT-RC', 'G']
    Ht = prop.p_s(streams.at['RC-MIX', 'P'], streams.at['LTH-SPLT', 'S'], Gas)['H']
    streams.at['RC-MIX', 'H'] = streams.at['LTH-SPLT', 'H'] + (Ht - streams.at['LTH-SPLT', 'H']) / KPDcomp
    streams.at['RC-MIX', 'T'] = \
    prop.h_p(streams.at['RC-MIX', 'H'], streams.at['RC-MIX', 'P'], streams.at['RC-MIX', 'X'])["T"]
    streams.at['RC-MIX', 'Q'] = \
    prop.h_p(streams.at['RC-MIX', 'H'], streams.at['RC-MIX', 'P'], streams.at['RC-MIX', 'X'])["Q"]
    streams.at['RC-MIX', 'S'] = \
    prop.h_p(streams.at['RC-MIX', 'H'], streams.at['RC-MIX', 'P'], streams.at['RC-MIX', 'X'])["S"]

    # MIX
    streams.at['MIX-HTH', 'P'] = Pk * pi - dPpot
    streams.at['MIX-HTH', 'X'] = Gas
    streams.at['MIX-HTH', 'G'] = streams.at['RC-MIX', 'G'] + streams.at['MC-LTH', 'G']
    streams.at['MIX-HTH', 'H'] = (streams.at['RC-MIX', 'G'] * streams.at['RC-MIX', 'H'] + streams.at['LTH-MIX', 'G'] *
                                  streams.at['LTH-MIX', 'H']) / streams.at['MIX-HTH', 'G']
    streams.at['MIX-HTH', 'T'] = prop.h_p(streams.at['MIX-HTH', 'H'], streams.at['MIX-HTH', 'P'], Gas)["T"]
    streams.at['MIX-HTH', 'Q'] = prop.h_p(streams.at['MIX-HTH', 'H'], streams.at['MIX-HTH', 'P'], Gas)["Q"]
    streams.at['MIX-HTH', 'S'] = prop.h_p(streams.at['MIX-HTH', 'H'], streams.at['MIX-HTH', 'P'], Gas)["S"]

    # ht cold
    streams.at['HTH-Heat', 'P'] = Pk * pi - 2 * dPpot
    streams.at['HTH-Heat', 'X'] = Gas
    streams.at['HTH-Heat', 'G'] = streams.at['MIX-HTH', 'G']
    streams.at['HTH-Heat', 'H'] = streams.at['Turb-HTH', 'H'] - streams.at['HTH-LTH', 'H'] + streams.at['MIX-HTH', 'H']
    streams.at['HTH-Heat', 'T'] = prop.h_p(streams.at['HTH-Heat', 'H'], streams.at['HTH-Heat', 'P'], Gas)['T']
    streams.at['HTH-Heat', 'Q'] = prop.h_p(streams.at['HTH-Heat', 'H'], streams.at['HTH-Heat', 'P'], Gas)['Q']
    streams.at['HTH-Heat', 'S'] = prop.h_p(streams.at['HTH-Heat', 'H'], streams.at['HTH-Heat', 'P'], Gas)['S']

    NCO2_MC = streams.at['Cool-MC', 'G'] * (streams.at['MC-LTH', 'H'] - streams.at['Cool-MC', 'H'])
    NCO2_RC = streams.at['RC-MIX', 'G'] * (streams.at['RC-MIX', 'H'] - streams.at['LTH-SPLT', 'H'])
    NCO2_Turb = G0 * (streams.at['Heat-Turb', 'H'] - streams.at['Turb-HTH', 'H'])
    NCO2_Heat = G0 * (streams.at['Heat-Turb', 'H'] - streams.at['HTH-Heat', 'H'])
    KPD = 0.99*0.99*(NCO2_Turb-(NCO2_MC+NCO2_RC))/NCO2_Heat*100

    # TS
    m = 50
    # Расширение в турбине
    Hturb = np.zeros(m)
    Tturb = np.zeros(m)
    Sturb = np.zeros(m)
    Hturb[0] = streams.at['Heat-Turb', 'H']
    S0 = streams.at['Heat-Turb', 'S']
    Pturb = np.linspace(streams.at['Heat-Turb', 'P'], streams.at['Turb-HTH', 'P'], m)
    for i in range(0, m):
        Hturb[i] = Hturb[0] - (Hturb[0] - prop.p_s(Pturb[i], S0, 'CO2')['H']) * KPDturb
        Sturb[i] = prop.h_p(Hturb[i], Pturb[i], 'CO2')['S']
        Tturb[i] = prop.h_p(Hturb[i], Pturb[i], 'CO2')['T']
    # Сжатие в компрессоре MC
    Hcomp1 = np.zeros(m)
    Tcomp1 = np.zeros(m)
    Scomp1 = np.zeros(m)
    Hcomp1[0] = streams.at['Cool-MC', 'H']
    S0 = streams.at['Cool-MC', 'S']
    Pcomp1 = np.linspace(streams.at['Cool-MC', 'P'], streams.at['MC-LTH', 'P'], m)
    for i in range(0, m):
        Hcomp1[i] = Hcomp1[0] + (prop.p_s(Pcomp1[i], S0, 'CO2')['H'] - Hcomp1[0]) / KPDcomp
        Scomp1[i] = prop.h_p(Hcomp1[i], Pcomp1[i], 'CO2')['S']
        Tcomp1[i] = prop.h_p(Hcomp1[i], Pcomp1[i], 'CO2')['T']
    # Сжатие в компрессоре RC
    Hcomp2 = np.zeros(m)
    Tcomp2 = np.zeros(m)
    Scomp2 = np.zeros(m)
    Hcomp2[0] = streams.at['LTH-SPLT', 'H']
    S0 = streams.at['LTH-SPLT', 'S']
    Pcomp2 = np.linspace(streams.at['LTH-SPLT', 'P'], streams.at['RC-MIX', 'P'], m)
    for i in range(0, m):
        Hcomp2[i] = Hcomp2[0] + (prop.p_s(Pcomp2[i], S0, 'CO2')['H'] - Hcomp2[0]) / KPDcomp
        Scomp2[i] = prop.h_p(Hcomp2[i], Pcomp2[i], 'CO2')['S']
        Tcomp2[i] = prop.h_p(Hcomp2[i], Pcomp2[i], 'CO2')['T']
    # heat
    Pheat = np.linspace(streams.at['HTH-Heat', 'P'], streams.at['Heat-Turb', 'P'], m)
    Hheat = np.linspace(streams.at['HTH-Heat', 'H'], streams.at['Heat-Turb', 'H'], m)
    Theat = np.zeros(m)
    Sheat = np.zeros(m)
    for i in range(0, m):
        Theat[i] = prop.h_p(Hheat[i], Pheat[i], 'CO2')['T']
        Sheat[i] = prop.h_p(Hheat[i], Pheat[i], 'CO2')['S']
    Pheat2 = np.linspace(streams.at['MIX-HTH', 'P'], streams.at['HTH-Heat', 'P'], m)
    Hheat2 = np.linspace(streams.at['MIX-HTH', 'H'], streams.at['HTH-Heat', 'H'], m)
    Theat2 = np.zeros(m)
    Sheat2 = np.zeros(m)
    for i in range(0, m):
        Theat2[i] = prop.h_p(Hheat2[i], Pheat2[i], 'CO2')['T']
        Sheat2[i] = prop.h_p(Hheat2[i], Pheat2[i], 'CO2')['S']
    Pheat3 = np.linspace(streams.at['MC-LTH', 'P'], streams.at['LTH-MIX', 'P'], m)
    Hheat3 = np.linspace(streams.at['MC-LTH', 'H'], streams.at['LTH-MIX', 'H'], m)
    Theat3 = np.zeros(m)
    Sheat3 = np.zeros(m)
    for i in range(0, m):
        Theat3[i] = prop.h_p(Hheat3[i], Pheat3[i], 'CO2')['T']
        Sheat3[i] = prop.h_p(Hheat3[i], Pheat3[i], 'CO2')['S']
    # # Холодный источние
    Pcool = np.linspace(streams.at['Cool-MC', 'P'], streams.at['LTH-SPLT', 'P'], m)
    Hcool = np.linspace(streams.at['Cool-MC', 'H'], streams.at['LTH-SPLT', 'H'], m)
    Tcool = np.zeros(m)
    Scool = np.zeros(m)
    for i in range(0, m):
        Tcool[i] = prop.h_p(Hcool[i], Pcool[i], 'CO2')['T']
        Scool[i] = prop.h_p(Hcool[i], Pcool[i], 'CO2')['S']
    Pcool1 = np.linspace(streams.at['LTH-SPLT', 'P'], streams.at['HTH-LTH', 'P'], m)
    Hcool1 = np.linspace(streams.at['LTH-SPLT', 'H'], streams.at['HTH-LTH', 'H'], m)
    Tcool1 = np.zeros(m)
    Scool1 = np.zeros(m)
    for i in range(0, m):
        Tcool1[i] = prop.h_p(Hcool1[i], Pcool1[i], 'CO2')['T']
        Scool1[i] = prop.h_p(Hcool1[i], Pcool1[i], 'CO2')['S']
    Pcool2 = np.linspace(streams.at['HTH-LTH', 'P'], streams.at['Turb-HTH', 'P'], m)
    Hcool2 = np.linspace(streams.at['HTH-LTH', 'H'], streams.at['Turb-HTH', 'H'], m)
    Tcool2 = np.zeros(m)
    Scool2 = np.zeros(m)
    for i in range(0, m):
        Tcool2[i] = prop.h_p(Hcool2[i], Pcool2[i], 'CO2')['T']
        Scool2[i] = prop.h_p(Hcool2[i], Pcool2[i], 'CO2')['S']
    # линия насыщения
    Tsat = np.linspace(1, 30, m)
    Tsat = np.append(Tsat, [30.5, 30.6, 30.7, 30.8, 30.9, 30.95, 30.96, 30.97, 30.9782])
    Ssat0 = np.zeros(len(Tsat))
    Ssat1 = np.zeros(len(Tsat))
    Tsats = np.append(Tsat, np.flip(Tsat))
    for i in range(0, len(Tsat)):
        Ssat0[i] = prop.t_q(Tsat[i], 0, 'CO2')['S']
        Ssat1[i] = prop.t_q(Tsat[i], 1, 'CO2')['S']
    Ssat = np.append(Ssat0, np.flip(Ssat1))


    plt.plot(Sturb, Tturb, color='navy')
    plt.plot(Scomp1, Tcomp1, color='navy')
    plt.plot(Scomp2, Tcomp2, color='navy')
    plt.plot(Sheat, Theat, color='navy')
    plt.plot(Sheat2, Theat2, color='orange')
    plt.plot(Sheat3, Theat3, color='orange')
    plt.plot(Scool, Tcool, color='navy')
    plt.plot(Scool1, Tcool1, color='orange')
    plt.plot(Scool2, Tcool2, color='orange')
    plt.plot(Ssat, Tsats, color='red')
    plt.show()
    return KPD

print(Calc(4,0.6,7.8,500))






















# N = 3
# M = 3
# pi = np.linspace(1.1,10,N)
# x = np.linspace(0.5,0.9,M)
# T1 = np.array([490,590,690,790])
# KPD1 = np.zeros((N,M),dtype='float32')
# Tres = np.zeros((N,M),dtype='float32')
# for i in range(0,len(T1)):
#     for k in range(0,M):
#         for j in range(0,N):
#             res = Calc(pi[j],x[k],T1[i])
#             KPD1[j, k] = res[0]
#             Tres[j, k] = res[1]
#     df = pd.DataFrame(KPD1)
#     df2 = pd.DataFrame(Tres)
#     with pd.ExcelWriter("calc-recomp-res.xlsx",if_sheet_exists='replace',mode='a') as writer:
#         df.to_excel(writer,sheet_name=str(T1[i]))
#     with pd.ExcelWriter("calc-recomp-res.xlsx",if_sheet_exists='overlay',mode='a') as writer2:
#         df2.to_excel(writer2,sheet_name=str(T1[i]),startcol=M+1)
N = 2
M = 2
# global data
# data = open("res490.txt", "a")
# pi = np.linspace(1.1,7,N)
# x = np.linspace(0.5,0.9,N)
# pk = np.linspace(7.6,8,M)
# T1 = np.array([490,590,690,790])
# for k in range(0,N):
#     for j in range(0,N):
#         for i in range(0,M):
#             res = Calc(pi[j], x[k], pk[i], T1[0])
#             print(res)





