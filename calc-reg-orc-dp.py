import pandas as pd
import prop
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

streams = pd.read_excel('streams.xlsx', index_col=0, sheet_name='reg')


def Calc(pi, pi_orc, G_orc, H_percent):
    Tmax = 600
    Pk = 7.8
    # pi=12
    KPDcomp = 0.9
    KPDturb = 0.9
    Tmin = 32
    dTreg = 50
    dPpot = 0.01
    dPpot_orc = 0
    # cool
    streams.at['Cool-Comp', 'T'] = Tmin
    streams.at['Cool-Comp', 'P'] = Pk
    streams.at['Cool-Comp', 'X'] = "CO2"
    streams.at['Cool-Comp', 'H'] = \
        prop.t_p(streams.at['Cool-Comp', 'T'], streams.at['Cool-Comp', 'P'], streams.at['Cool-Comp', 'X'])['H']
    streams.at['Cool-Comp', 'G'] = 1
    streams.at['Cool-Comp', 'Q'] = 1
    streams.at['Cool-Comp', 'S'] = \
        prop.h_p(streams.at['Cool-Comp', 'H'], streams.at['Cool-Comp', 'P'], streams.at['Cool-Comp', 'X'])['S']

    # comp
    streams.at['Comp-Reg', 'P'] = streams.at['Cool-Comp', 'P'] * pi
    streams.at['Comp-Reg', 'X'] = streams.at['Cool-Comp', 'X']
    Ht = prop.p_s(streams.at['Comp-Reg', 'P'], streams.at['Cool-Comp', 'S'], streams.at['Comp-Reg', 'X'])['H']
    streams.at['Comp-Reg', 'H'] = streams.at['Cool-Comp', 'H'] + (Ht - streams.at['Cool-Comp', 'H']) / KPDcomp
    streams.at['Comp-Reg', 'G'] = streams.at['Cool-Comp', 'G']
    streams.at['Comp-Reg', 'T'] = \
        prop.h_p(streams.at['Comp-Reg', 'H'], streams.at['Comp-Reg', 'P'], streams.at['Comp-Reg', 'X'])["T"]
    streams.at['Comp-Reg', 'Q'] = \
        prop.h_p(streams.at['Comp-Reg', 'H'], streams.at['Comp-Reg', 'P'], streams.at['Comp-Reg', 'X'])["Q"]
    streams.at['Comp-Reg', 'S'] = \
        prop.h_p(streams.at['Comp-Reg', 'H'], streams.at['Comp-Reg', 'P'], streams.at['Comp-Reg', 'X'])["S"]

    dP_reg_cold = dPpot
    # reg_cold
    streams.at['Reg-Heat', 'P'] = streams.at['Comp-Reg', 'P'] - dP_reg_cold
    streams.at['Reg-Heat', 'X'] = streams.at['Comp-Reg', 'X']
    streams.at['Reg-Heat', 'G'] = streams.at['Comp-Reg', 'G']

    dP_heater = dPpot
    # #heat
    streams.at['Heat-Turb', 'T'] = Tmax
    streams.at['Heat-Turb', 'P'] = streams.at['Reg-Heat', 'P'] - dP_heater
    streams.at['Heat-Turb', 'X'] = streams.at['Reg-Heat', 'X']
    streams.at['Heat-Turb', 'G'] = streams.at['Reg-Heat', 'G']
    streams.at['Heat-Turb', 'H'] = \
        prop.t_p(streams.at['Heat-Turb', 'T'], streams.at['Heat-Turb', 'P'], streams.at['Heat-Turb', 'X'])['H']
    streams.at['Heat-Turb', 'Q'] = \
        prop.t_p(streams.at['Heat-Turb', 'T'], streams.at['Heat-Turb', 'P'], streams.at['Heat-Turb', 'X'])['Q']
    streams.at['Heat-Turb', 'S'] = \
        prop.t_p(streams.at['Heat-Turb', 'T'], streams.at['Heat-Turb', 'P'], streams.at['Heat-Turb', 'X'])['S']

    dP_reg_hot = dPpot
    dP_evap_hot = dPpot
    dP_cool = dPpot
    # turb
    streams.at['Turb-Reg', 'P'] = Pk + dP_reg_hot + dP_evap_hot + dP_cool
    streams.at['Turb-Reg', 'X'] = streams.at['Heat-Turb', 'X']
    streams.at['Turb-Reg', 'G'] = streams.at['Heat-Turb', 'G']
    Ht = prop.p_s(streams.at['Turb-Reg', 'P'], streams.at['Heat-Turb', 'S'], streams.at['Turb-Reg', 'X'])['H']
    streams.at['Turb-Reg', 'H'] = streams.at['Heat-Turb', 'H'] - (streams.at['Heat-Turb', 'H'] - Ht) * KPDturb
    streams.at['Turb-Reg', 'T'] = \
        prop.h_p(streams.at['Turb-Reg', 'H'], streams.at['Turb-Reg', 'P'], streams.at['Turb-Reg', 'X'])["T"]
    streams.at['Turb-Reg', 'Q'] = \
        prop.h_p(streams.at['Turb-Reg', 'H'], streams.at['Turb-Reg', 'P'], streams.at['Turb-Reg', 'X'])["Q"]
    streams.at['Turb-Reg', 'S'] = \
        prop.h_p(streams.at['Turb-Reg', 'H'], streams.at['Turb-Reg', 'P'], streams.at['Turb-Reg', 'X'])["S"]

    # reg_hot
    streams.at['Reg-Evap', 'T'] = streams.at['Comp-Reg', 'T'] + dTreg
    streams.at['Reg-Evap', 'P'] = streams.at['Turb-Reg', 'P'] - dP_reg_hot
    streams.at['Reg-Evap', 'X'] = streams.at['Turb-Reg', 'X']
    streams.at['Reg-Evap', 'G'] = streams.at['Turb-Reg', 'G']
    streams.at['Reg-Evap', 'H'] = \
        prop.t_p(streams.at['Reg-Evap', 'T'], streams.at['Reg-Evap', 'P'], streams.at['Reg-Evap', 'X'])['H']
    streams.at['Reg-Evap', 'Q'] = \
        prop.t_p(streams.at['Reg-Evap', 'T'], streams.at['Reg-Evap', 'P'], streams.at['Reg-Evap', 'X'])['Q']
    streams.at['Reg-Evap', 'S'] = \
        prop.t_p(streams.at['Reg-Evap', 'T'], streams.at['Reg-Evap', 'P'], streams.at['Reg-Evap', 'X'])['S']

    # reg_cold
    streams.at['Reg-Heat', 'H'] = streams.at['Comp-Reg', 'H'] + (
            streams.at['Turb-Reg', 'H'] - streams.at['Reg-Evap', 'H'])
    streams.at['Reg-Heat', 'T'] = \
        prop.h_p(streams.at['Reg-Heat', 'H'], streams.at['Reg-Heat', 'P'], streams.at['Reg-Heat', 'X'])['T']
    streams.at['Reg-Heat', 'Q'] = \
        prop.h_p(streams.at['Reg-Heat', 'H'], streams.at['Reg-Heat', 'P'], streams.at['Reg-Heat', 'X'])['Q']
    streams.at['Reg-Heat', 'S'] = \
        prop.h_p(streams.at['Reg-Heat', 'H'], streams.at['Reg-Heat', 'P'], streams.at['Reg-Heat', 'X'])['S']

    # H_percent = 0.6
    # evap
    streams.at['Evap-Cool', 'H'] = (streams.at['Reg-Evap', 'H'] - streams.at['Cool-Comp', 'H']) * H_percent + \
                                   streams.at['Cool-Comp', 'H']
    streams.at['Evap-Cool', 'P'] = Pk + dP_cool
    streams.at['Evap-Cool', 'X'] = streams.at['Reg-Evap', 'X']
    streams.at['Evap-Cool', 'G'] = streams.at['Reg-Evap', 'G']
    streams.at['Evap-Cool', 'T'] = \
        prop.h_p(streams.at['Evap-Cool', 'H'], streams.at['Evap-Cool', 'P'], streams.at['Evap-Cool', 'X'])["T"]
    streams.at['Evap-Cool', 'Q'] = \
        prop.h_p(streams.at['Evap-Cool', 'H'], streams.at['Evap-Cool', 'P'], streams.at['Evap-Cool', 'X'])["Q"]
    streams.at['Evap-Cool', 'S'] = \
        prop.h_p(streams.at['Evap-Cool', 'H'], streams.at['Evap-Cool', 'P'], streams.at['Evap-Cool', 'X'])["S"]

    Tmin_orc = 30
    Fluid_orc = "R236ea"
    Pk_orc = prop.t_q(Tmin_orc, 0, Fluid_orc)["P"]
    # G_orc = 0.66
    # ORC Cool
    streams.at['OCool-OPump', 'P'] = Pk_orc
    streams.at['OCool-OPump', 'T'] = Tmin_orc
    streams.at['OCool-OPump', 'X'] = Fluid_orc
    streams.at['OCool-OPump', 'G'] = G_orc
    streams.at['OCool-OPump', 'H'] = \
        prop.t_p(streams.at['OCool-OPump', 'T'], streams.at['OCool-OPump', 'P'], streams.at['OCool-OPump', 'X'])[
            "H"]
    streams.at['OCool-OPump', 'Q'] = \
        prop.t_p(streams.at['OCool-OPump', 'T'], streams.at['OCool-OPump', 'P'], streams.at['OCool-OPump', 'X'])[
            "Q"]
    streams.at['OCool-OPump', 'S'] = \
        prop.t_p(streams.at['OCool-OPump', 'T'], streams.at['OCool-OPump', 'P'], streams.at['OCool-OPump', 'X'])[
            "S"]


    KPDpump = KPDcomp
    # ORC Pump
    streams.at['OPump-OEvap', 'P'] = streams.at['OCool-OPump', 'P'] * pi_orc
    streams.at['OPump-OEvap', 'X'] = streams.at['OCool-OPump', 'X']
    Ht = prop.p_s(streams.at['OPump-OEvap', 'P'], streams.at['OCool-OPump', 'S'], streams.at['OPump-OEvap', 'X'])[
        'H']
    streams.at['OPump-OEvap', 'H'] = streams.at['OCool-OPump', 'H'] + (
            Ht - streams.at['OCool-OPump', 'H']) / KPDpump
    streams.at['OPump-OEvap', 'G'] = streams.at['OCool-OPump', 'G']
    streams.at['OPump-OEvap', 'T'] = \
        prop.h_p(streams.at['OPump-OEvap', 'H'], streams.at['OPump-OEvap', 'P'], streams.at['OPump-OEvap', 'X'])[
            "T"]
    streams.at['OPump-OEvap', 'Q'] = \
        prop.h_p(streams.at['OPump-OEvap', 'H'], streams.at['OPump-OEvap', 'P'], streams.at['OPump-OEvap', 'X'])[
            "Q"]
    streams.at['OPump-OEvap', 'S'] = \
        prop.h_p(streams.at['OPump-OEvap', 'H'], streams.at['OPump-OEvap', 'P'], streams.at['OPump-OEvap', 'X'])[
            "S"]

    # ORC Evap
    streams.at['OEvap-OTurb', 'P'] = streams.at['OPump-OEvap', 'P'] - dPpot_orc
    streams.at['OEvap-OTurb', 'H'] = (streams.at['Cool-Comp', 'G'] * (
            streams.at['Reg-Evap', 'H'] - streams.at['Evap-Cool', 'H'])) / G_orc + streams.at[
                                         'OPump-OEvap', 'H']
    streams.at['OEvap-OTurb', 'X'] = "R236ea"
    streams.at['OEvap-OTurb', 'G'] = G_orc
    streams.at['OEvap-OTurb', 'T'] = \
        prop.h_p(streams.at['OEvap-OTurb', 'H'], streams.at['OEvap-OTurb', 'P'], streams.at['OEvap-OTurb', 'X'])[
            "T"]
    streams.at['OEvap-OTurb', 'Q'] = \
        prop.h_p(streams.at['OEvap-OTurb', 'H'], streams.at['OEvap-OTurb', 'P'], streams.at['OEvap-OTurb', 'X'])[
            "Q"]
    streams.at['OEvap-OTurb', 'S'] = \
        prop.h_p(streams.at['OEvap-OTurb', 'H'], streams.at['OEvap-OTurb', 'P'], streams.at['OEvap-OTurb', 'X'])[
            "S"]

    KPDturb_orc = KPDturb
    # ORC Turb
    streams.at['OTurb-OCool', 'P'] = streams.at['OCool-OPump', 'P']
    streams.at['OTurb-OCool', 'X'] = streams.at['OCool-OPump', 'X']
    streams.at['OTurb-OCool', 'G'] = streams.at['OCool-OPump', 'G']
    Ht = prop.p_s(streams.at['OTurb-OCool', 'P'], streams.at['OEvap-OTurb', 'S'], streams.at['OTurb-OCool', 'X'])[
        'H']
    streams.at['OTurb-OCool', 'H'] = streams.at['OEvap-OTurb', 'H'] - (
            streams.at['OEvap-OTurb', 'H'] - Ht) * KPDturb_orc
    streams.at['OTurb-OCool', 'T'] = \
        prop.h_p(streams.at['OTurb-OCool', 'H'], streams.at['OTurb-OCool', 'P'], streams.at['OTurb-OCool', 'X'])[
            "T"]
    streams.at['OTurb-OCool', 'Q'] = \
        prop.h_p(streams.at['OTurb-OCool', 'H'], streams.at['OTurb-OCool', 'P'], streams.at['OTurb-OCool', 'X'])[
            "Q"]
    streams.at['OTurb-OCool', 'S'] = \
        prop.h_p(streams.at['OTurb-OCool', 'H'], streams.at['OTurb-OCool', 'P'], streams.at['OTurb-OCool', 'X'])[
            "S"]

    Q = streams.at['Reg-Heat', 'G'] * (streams.at['Heat-Turb', 'H'] - streams.at['Reg-Heat', 'H'])
    Nco2 = 0.99 * streams.at['Heat-Turb', 'G'] * (streams.at['Heat-Turb', 'H'] - streams.at['Turb-Reg', 'H']) - \
           streams.at['Cool-Comp', 'G'] * (streams.at['Comp-Reg', 'H'] - streams.at['Cool-Comp', 'H']) / 0.99
    Norc = 0.99 * streams.at['OTurb-OCool', 'G'] * (
            streams.at['OEvap-OTurb', 'H'] - streams.at['OTurb-OCool', 'H']) - streams.at[
               'OEvap-OTurb', 'G'] * (streams.at['OPump-OEvap', 'H'] - streams.at['OCool-OPump', 'H']) / 0.99
    KPD = (Nco2 + Norc) / Q * 100
    # Разбиение и проверка темп. напора
    n = 50
    T_co2_evap = np.zeros(n)
    H_co2_evap = np.linspace(streams.at['Reg-Evap', 'H'], streams.at['Evap-Cool', 'H'], n)
    for i in range(0, n):
        T_co2_evap[i] = prop.h_p(H_co2_evap[i], streams.at['Reg-Evap', 'P'], streams.at['Reg-Evap', 'X'])["T"]
    Q_co2_evap = -(H_co2_evap - streams.at['Reg-Evap', 'H']) * streams.at['Cool-Comp', 'G']

    T_orc_evap = np.zeros(n)
    H_orc_evap = np.linspace(streams.at['OEvap-OTurb', 'H'], streams.at['OPump-OEvap', 'H'], n)
    for i in range(0, n):
        T_orc_evap[i] = prop.h_p(H_orc_evap[i], streams.at['OPump-OEvap', 'P'], streams.at['OPump-OEvap', 'X'])["T"]
    Q_orc_evap = -(H_orc_evap - streams.at['OEvap-OTurb', 'H']) * streams.at['OEvap-OTurb', 'G']
    deltaT = T_co2_evap - T_orc_evap
    # print(streams)
    # print(min(deltaT))
    # print(Pk_orc)
    # plt.plot(Q_co2_evap,T_co2_evap)
    # plt.plot(Q_orc_evap,T_orc_evap)
    # plt.show()
    if min(deltaT) < 10:
        KPD = float(0)
    if streams.at['OEvap-OTurb', 'Q'] != 1:
        KPD = float(0)
    print(round(pi, 3),
          round(G_orc, 3),
          round(H_percent, 3),
          round(pi_orc, 3),
          round(KPD, 3))
    return -KPD

print(Calc(8, 5,0.6,0.6))




