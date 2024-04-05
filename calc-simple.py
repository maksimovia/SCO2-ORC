import pandas as pd
import prop
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
streams = pd.read_excel('streams.xlsx', index_col=0, sheet_name='simple')

def Calc(pi,Pk,Tmax):
    #Tmax = 500
    #Pk = 7.8
    KPDcomp = 0.9
    KPDturb = 0.9
    Tmin = 32
    dPpot = 0.05
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
    streams.at['Heat-Turb', 'P'] = streams.at['Comp-Heat', 'P'] - dPpot
    streams.at['Heat-Turb', 'X'] = streams.at['Comp-Heat', 'X']
    streams.at['Heat-Turb', 'H'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['H']
    streams.at['Heat-Turb', 'G'] = streams.at['Comp-Heat', 'G']
    streams.at['Heat-Turb', 'Q'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['Q']
    streams.at['Heat-Turb', 'S'] = prop.t_p(streams.at['Heat-Turb', 'T'],streams.at['Heat-Turb', 'P'],streams.at['Heat-Turb', 'X'])['S']

    #turb
    streams.at['Turb-Cool', 'P'] = Pk + dPpot
    streams.at['Turb-Cool', 'X'] = streams.at['Heat-Turb', 'X']
    streams.at['Turb-Cool', 'G'] = streams.at['Heat-Turb', 'G']
    Ht = prop.p_s(streams.at['Turb-Cool', 'P'], streams.at['Heat-Turb', 'S'], streams.at['Turb-Cool', 'X'])['H']
    streams.at['Turb-Cool', 'H'] = streams.at['Heat-Turb', 'H'] - (streams.at['Heat-Turb', 'H']-Ht)*KPDturb
    streams.at['Turb-Cool', 'T'] = prop.h_p(streams.at['Turb-Cool', 'H'],streams.at['Turb-Cool', 'P'],streams.at['Turb-Cool', 'X'])["T"]
    streams.at['Turb-Cool', 'Q'] = prop.h_p(streams.at['Turb-Cool', 'H'],streams.at['Turb-Cool', 'P'],streams.at['Turb-Cool', 'X'])["Q"]
    streams.at['Turb-Cool', 'S'] = prop.h_p(streams.at['Turb-Cool', 'H'],streams.at['Turb-Cool', 'P'],streams.at['Turb-Cool', 'X'])["S"]

    Q = streams.at['Comp-Heat', 'G']*(streams.at['Heat-Turb', 'H']-streams.at['Comp-Heat', 'H'])
    Nco2 =streams.at['Heat-Turb', 'G']*(streams.at['Heat-Turb', 'H']-streams.at['Turb-Cool', 'H']) - streams.at['Heat-Turb', 'G']*(streams.at['Comp-Heat', 'H']-streams.at['Cool-Comp', 'H'])
    KPD = 0.99*0.99*Nco2/Q*100
    T = streams.at['Turb-Cool', 'T']

    #TS
    m = 50
    #Расширение в турбине
    Hturb = np.zeros(m)
    Tturb = np.zeros(m)
    Sturb = np.zeros(m)
    Hturb[0] = streams.at['Heat-Turb', 'H']
    S0 = streams.at['Heat-Turb', 'S']
    Pturb = np.linspace(streams.at['Heat-Turb', 'P'],streams.at['Turb-Cool', 'P'],m)
    for i in range(0,m):
        Hturb[i] = Hturb[0] - (Hturb[0] - prop.p_s(Pturb[i],S0,'CO2')['H'])*KPDturb
        Sturb[i] = prop.h_p(Hturb[i],Pturb[i],'CO2')['S']
        Tturb[i] = prop.h_p(Hturb[i], Pturb[i], 'CO2')['T']
    #Сжатие в компрессоре
    Hcomp = np.zeros(m)
    Tcomp = np.zeros(m)
    Scomp = np.zeros(m)
    Hcomp[0] = streams.at['Cool-Comp', 'H']
    S0 = streams.at['Cool-Comp', 'S']
    Pcomp = np.linspace(streams.at['Cool-Comp', 'P'], streams.at['Comp-Heat', 'P'], m)
    for i in range(0, m):
        Hcomp[i] = Hcomp[0] + (prop.p_s(Pcomp[i], S0, 'CO2')['H'] - Hcomp[0]) / KPDcomp
        Scomp[i] = prop.h_p(Hcomp[i], Pcomp[i], 'CO2')['S']
        Tcomp[i] = prop.h_p(Hcomp[i], Pcomp[i], 'CO2')['T']
    #Прогрев гоев
    Pheat = np.linspace(streams.at['Comp-Heat', 'P'],streams.at['Heat-Turb', 'P'],m)
    Hheat = np.linspace(streams.at['Comp-Heat', 'H'],streams.at['Heat-Turb', 'H'],m)
    Theat = np.zeros(m)
    Sheat = np.zeros(m)
    for i in range(0,m):
        Theat[i] = prop.h_p(Hheat[i],Pheat[i],'CO2')['T']
        Sheat[i] = prop.h_p(Hheat[i], Pheat[i], 'CO2')['S']

    # Холодный источние
    Pcool = np.linspace(streams.at['Cool-Comp', 'P'], streams.at['Turb-Cool', 'P'], m)
    Hcool = np.linspace(streams.at['Cool-Comp', 'H'], streams.at['Turb-Cool', 'H'], m)
    Tcool = np.zeros(m)
    Scool= np.zeros(m)
    for i in range(0, m):
        Tcool[i] = prop.h_p(Hcool[i], Pcool[i], 'CO2')['T']
        Scool[i] = prop.h_p(Hcool[i], Pcool[i], 'CO2')['S']
    #линия насыщения
    Tsat = np.linspace(1, 30, m)
    Tsat=np.append(Tsat,[30.5, 30.6,30.7,30.8,30.9,30.95,30.96,30.97,30.9782])
    Ssat0 = np.zeros(len(Tsat))
    Ssat1 = np.zeros(len(Tsat))
    Tsats = np.append(Tsat,np.flip(Tsat))
    for i in range(0,len(Tsat)):
        Ssat0[i] = prop.t_q(Tsat[i],0,'CO2')['S']
        Ssat1[i] = prop.t_q(Tsat[i],1,'CO2')['S']
    Ssat = np.append(Ssat0,np.flip(Ssat1))

    print(Ssat0)
    plt.plot(Sturb, Tturb,color='navy')
    plt.plot(Scomp, Tcomp,color='navy')
    plt.plot(Sheat, Theat,color='navy')
    plt.plot(Scool, Tcool,color='navy')
    plt.plot(Ssat, Tsats,color='red')
    plt.show()

    return KPD, T
Calc(4,7.8,500)
# pi = np.linspace(1,30,N)
# pk = np.linspace(7.5,8.5,M)
# T = np.array([490,590,690,790])
# KPD = np.zeros((N,M),dtype='float32')
# Tres = np.zeros((N,M),dtype='float32')
# for i in range(0,len(T)):
#     for k in range(0,M):
#         for j in range(0,N):
#             res = Calc(pi[j],pk[k],T[i])
#             KPD[j, k] = res[0]
#             Tres[j, k] = res[1]
#     df = pd.DataFrame(KPD)
#     df2 = pd.DataFrame(Tres)
#     with pd.ExcelWriter("calc-simple-res.xlsx",if_sheet_exists='replace',mode='a') as writer:
#         df.to_excel(writer,sheet_name=str(T[i]))
#     with pd.ExcelWriter("calc-simple-res.xlsx",if_sheet_exists='overlay',mode='a') as writer2:
#         df2.to_excel(writer2,sheet_name=str(T[i]),startcol=M+1)
# N = 200
# M = 20
# data = open("simple790.txt", "a")
# pi = np.linspace(1.1,7,N)
# pk = np.linspace(7.5,8.5,M)
# T1 = np.array([490,590,690,790])
# for j in range(0,N):
#     for i in range(0,M):
#         res = Calc(pi[j], pk[i], T1[3])
#         print(res)
#         data.write('\n' + str(round(pi[j], 5))  + '\t' + str(
#             round(pk[i], 5)) + '\t' + str(round(res[1], 5)) + '\t' + str(round(res[0], 5)))





