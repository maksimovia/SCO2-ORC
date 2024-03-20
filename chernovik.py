import pandas as pd
import prop
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
streams = pd.read_excel('streams.xlsx', index_col=0, sheet_name='simple')





pi_orc=5
def CalcPP(pi,pi_orc):
    #pi = 2
    #pi_orc = 10
    def CalcGH(Input):
        G_orc = Input[0]
        H_percent = Input[1]
        Tmax = 600
        Pk = 7.8

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
        streams.at['Turb-Evap', 'P'] = Pk
        streams.at['Turb-Evap', 'X'] = streams.at['Heat-Turb', 'X']
        streams.at['Turb-Evap', 'G'] = streams.at['Heat-Turb', 'G']
        Ht = prop.p_s(streams.at['Turb-Evap', 'P'], streams.at['Heat-Turb', 'S'], streams.at['Turb-Evap', 'X'])['H']
        streams.at['Turb-Evap', 'H'] = streams.at['Heat-Turb', 'H'] - (streams.at['Heat-Turb', 'H']-Ht)*KPDturb
        streams.at['Turb-Evap', 'T'] = prop.h_p(streams.at['Turb-Evap', 'H'],streams.at['Turb-Evap', 'P'],streams.at['Turb-Evap', 'X'])["T"]
        streams.at['Turb-Evap', 'Q'] = prop.h_p(streams.at['Turb-Evap', 'H'],streams.at['Turb-Evap', 'P'],streams.at['Turb-Evap', 'X'])["Q"]
        streams.at['Turb-Evap', 'S'] = prop.h_p(streams.at['Turb-Evap', 'H'],streams.at['Turb-Evap', 'P'],streams.at['Turb-Evap', 'X'])["S"]

        #H_percent = 0.6
        #evap
        streams.at['Evap-Cool', 'H'] = (streams.at['Turb-Evap', 'H']-streams.at['Cool-Comp', 'H'])*H_percent+streams.at['Cool-Comp', 'H']
        streams.at['Evap-Cool', 'P'] = Pk
        streams.at['Evap-Cool', 'X'] = streams.at['Turb-Evap', 'X']
        streams.at['Evap-Cool', 'G'] = streams.at['Turb-Evap', 'G']
        streams.at['Evap-Cool', 'T'] = prop.h_p(streams.at['Evap-Cool', 'H'],streams.at['Evap-Cool', 'P'],streams.at['Evap-Cool', 'X'])["T"]
        streams.at['Evap-Cool', 'Q'] = prop.h_p(streams.at['Evap-Cool', 'H'],streams.at['Evap-Cool', 'P'],streams.at['Evap-Cool', 'X'])["Q"]
        streams.at['Evap-Cool', 'S'] = prop.h_p(streams.at['Evap-Cool', 'H'],streams.at['Evap-Cool', 'P'],streams.at['Evap-Cool', 'X'])["S"]


        Tmin_orc = 30
        Fluid_orc = "R236ea"
        Pk_orc = prop.t_q(Tmin_orc,0,Fluid_orc)["P"]
        #G_orc = 0.5
        #ORC Cool
        streams.at['OCool-OPump', 'P'] = Pk_orc
        streams.at['OCool-OPump', 'T'] = Tmin_orc
        streams.at['OCool-OPump', 'X'] = Fluid_orc
        streams.at['OCool-OPump', 'G'] = G_orc
        streams.at['OCool-OPump', 'H'] = prop.t_p(streams.at['OCool-OPump', 'T'],streams.at['OCool-OPump', 'P'],streams.at['OCool-OPump', 'X'])["H"]
        streams.at['OCool-OPump', 'Q'] = prop.t_p(streams.at['OCool-OPump', 'T'],streams.at['OCool-OPump', 'P'],streams.at['OCool-OPump', 'X'])["Q"]
        streams.at['OCool-OPump', 'S'] = prop.t_p(streams.at['OCool-OPump', 'T'],streams.at['OCool-OPump', 'P'],streams.at['OCool-OPump', 'X'])["S"]


        KPDpump = KPDcomp
        #ORC Pump
        streams.at['OPump-OEvap', 'P'] = streams.at['OCool-OPump', 'P']*pi_orc
        streams.at['OPump-OEvap', 'X'] = streams.at['OCool-OPump', 'X']
        Ht = prop.p_s(streams.at['OPump-OEvap', 'P'], streams.at['OCool-OPump', 'S'], streams.at['OPump-OEvap', 'X'])['H']
        streams.at['OPump-OEvap', 'H'] = streams.at['OCool-OPump', 'H'] + (Ht - streams.at['OCool-OPump', 'H'])/KPDpump
        streams.at['OPump-OEvap', 'G'] = streams.at['OCool-OPump', 'G']
        streams.at['OPump-OEvap', 'T'] = prop.h_p(streams.at['OPump-OEvap', 'H'],streams.at['OPump-OEvap', 'P'],streams.at['OPump-OEvap', 'X'])["T"]
        streams.at['OPump-OEvap', 'Q'] = prop.h_p(streams.at['OPump-OEvap', 'H'],streams.at['OPump-OEvap', 'P'],streams.at['OPump-OEvap', 'X'])["Q"]
        streams.at['OPump-OEvap', 'S'] = prop.h_p(streams.at['OPump-OEvap', 'H'],streams.at['OPump-OEvap', 'P'],streams.at['OPump-OEvap', 'X'])["S"]

        #ORC Evap
        streams.at['OEvap-OTurb', 'P'] = streams.at['OPump-OEvap', 'P']
        streams.at['OEvap-OTurb', 'H'] = (streams.at['Cool-Comp', 'G']*(streams.at['Turb-Evap', 'H']-streams.at['Evap-Cool', 'H']))/G_orc + streams.at['OPump-OEvap', 'H']
        streams.at['OEvap-OTurb', 'X'] = "R236ea"
        streams.at['OEvap-OTurb', 'G'] = G_orc
        streams.at['OEvap-OTurb', 'T'] = prop.h_p(streams.at['OEvap-OTurb', 'H'],streams.at['OEvap-OTurb', 'P'],streams.at['OEvap-OTurb', 'X'])["T"]
        streams.at['OEvap-OTurb', 'Q'] = prop.h_p(streams.at['OEvap-OTurb', 'H'],streams.at['OEvap-OTurb', 'P'],streams.at['OEvap-OTurb', 'X'])["Q"]
        streams.at['OEvap-OTurb', 'S'] = prop.h_p(streams.at['OEvap-OTurb', 'H'],streams.at['OEvap-OTurb', 'P'],streams.at['OEvap-OTurb', 'X'])["S"]

        KPDturb_orc = KPDturb
        #ORC Turb
        streams.at['OTurb-OCool', 'P'] = streams.at['OCool-OPump', 'P']
        streams.at['OTurb-OCool', 'X'] = streams.at['OCool-OPump', 'X']
        streams.at['OTurb-OCool', 'G'] = streams.at['OCool-OPump', 'G']
        Ht = prop.p_s(streams.at['OTurb-OCool', 'P'], streams.at['OEvap-OTurb', 'S'], streams.at['OTurb-OCool', 'X'])['H']
        streams.at['OTurb-OCool', 'H'] = streams.at['OEvap-OTurb', 'H'] - (streams.at['OEvap-OTurb', 'H']-Ht)*KPDturb_orc
        streams.at['OTurb-OCool', 'T'] = prop.h_p(streams.at['OTurb-OCool', 'H'],streams.at['OTurb-OCool', 'P'],streams.at['OTurb-OCool', 'X'])["T"]
        streams.at['OTurb-OCool', 'Q'] = prop.h_p(streams.at['OTurb-OCool', 'H'],streams.at['OTurb-OCool', 'P'],streams.at['OTurb-OCool', 'X'])["Q"]
        streams.at['OTurb-OCool', 'S'] = prop.h_p(streams.at['OTurb-OCool', 'H'],streams.at['OTurb-OCool', 'P'],streams.at['OTurb-OCool', 'X'])["S"]

        Q = streams.at['Comp-Heat', 'G']*(streams.at['Heat-Turb', 'H']-streams.at['Comp-Heat', 'H'])
        Nco2 =0.99*streams.at['Heat-Turb', 'G']*(streams.at['Heat-Turb', 'H']-streams.at['Turb-Evap', 'H']) - streams.at['Heat-Turb', 'G']*(streams.at['Comp-Heat', 'H']-streams.at['Cool-Comp', 'H'])/0.99
        Norc =0.99*streams.at['OTurb-OCool', 'G']*(streams.at['OEvap-OTurb', 'H']-streams.at['OTurb-OCool', 'H']) - streams.at['OEvap-OTurb', 'G']*(streams.at['OPump-OEvap', 'H']-streams.at['OCool-OPump', 'H'])/0.99
        KPD = (Nco2+Norc)/Q*100
        #Разбиение и проверка темп. напора
        n = 50
        T_co2_evap = np.zeros(n)
        H_co2_evap = np.linspace(streams.at['Turb-Evap', 'H'], streams.at['Evap-Cool', 'H'], n)
        for i in range(0,n):
            T_co2_evap[i] = prop.h_p(H_co2_evap[i],streams.at['Turb-Evap', 'P'],streams.at['Turb-Evap', 'X'])["T"]
        Q_co2_evap = -(H_co2_evap-streams.at['Turb-Evap', 'H'])*streams.at['Cool-Comp', 'G']

        T_orc_evap = np.zeros(n)
        H_orc_evap = np.linspace(streams.at['OEvap-OTurb', 'H'], streams.at['OPump-OEvap', 'H'], n)
        for i in range(0,n):
            T_orc_evap[i] = prop.h_p(H_orc_evap[i],streams.at['OPump-OEvap', 'P'],streams.at['OPump-OEvap', 'X'])["T"]
        Q_orc_evap = -(H_orc_evap-streams.at['OEvap-OTurb', 'H'])*streams.at['OEvap-OTurb', 'G']
        deltaT = T_co2_evap - T_orc_evap
        #print(streams)
        #print(min(deltaT))
        #print(Pk_orc)
        if min(deltaT)<10:
            KPD = 0
        if streams.at['OEvap-OTurb', 'Q'] != 1:
            KPD = 0
        print(round(pi,3),
              round(G_orc,3),
              round(H_percent,3),
              round(pi_orc,3),
              round(KPD,3))
        # plt.plot(Q_co2_evap,T_co2_evap)
        # plt.plot(Q_orc_evap,T_orc_evap)
        # plt.show()

        return -KPD

    res=opt.minimize(CalcGH,(0.66,0.66),method='Nelder-Mead', bounds=((0.2,5),(0.05,0.95)), tol=10**-2)
    KPD = -CalcGH(res.x)
    return res.x[0],res.x[1],KPD
n = 9
pi = np.linspace(2,10,n)
pi_orc = np.linspace(5,29,n)
kpd = np.zeros((n,n),dtype="float32")
Gopt = np.zeros((n,n),dtype="float32")
Hopt = np.zeros((n,n),dtype="float32")
for i in range(0,n):
    for j in range(0,n):
        res = CalcPP(pi[j], pi_orc[i])

        kpd[i, j] = res[2]
        Gopt[i, j] = res[0]
        Hopt[i, j] = res[1]
print(kpd)

plt.plot(pi,kpd[:,0],label=str(pi_orc[0]))
plt.plot(pi,kpd[:,1],label=str(pi_orc[1]))
plt.plot(pi,kpd[:,2],label=str(pi_orc[2]))
plt.plot(pi,kpd[:,3],label=str(pi_orc[3]))
# plt.plot(pi,kpd[:,4],label=str(pi_orc[4]))
# plt.plot(pi,kpd[:,5],label=str(pi_orc[5]))
# plt.plot(pi,kpd[:,6],label=str(pi_orc[6]))
# plt.plot(pi,kpd[:,7],label=str(pi_orc[7]))
# plt.plot(pi,kpd[:,8],label=str(pi_orc[8]))
plt.legend()
plt.show()

plt.plot(pi,Gopt[:,0],label=str(pi_orc[0]))
plt.plot(pi,Gopt[:,1],label=str(pi_orc[1]))
plt.plot(pi,Gopt[:,2],label=str(pi_orc[2]))
plt.plot(pi,Gopt[:,3],label=str(pi_orc[3]))
# plt.plot(pi,Gopt[:,4],label=str(pi_orc[4]))
# plt.plot(pi,Gopt[:,5],label=str(pi_orc[5]))
# plt.plot(pi,Gopt[:,6],label=str(pi_orc[6]))
# plt.plot(pi,Gopt[:,7],label=str(pi_orc[7]))
#plt.plot(pi,Gopt[:,8],label=str(pi_orc[8]))
plt.legend()
plt.show()

plt.plot(pi,Hopt[:,0],label=str(pi_orc[0]))
plt.plot(pi,Hopt[:,1],label=str(pi_orc[1]))
plt.plot(pi,Hopt[:,2],label=str(pi_orc[2]))
plt.plot(pi,Hopt[:,3],label=str(pi_orc[3]))
# plt.plot(pi,Hopt[:,4],label=str(pi_orc[4]))
# plt.plot(pi,Hopt[:,5],label=str(pi_orc[5]))
# plt.plot(pi,Hopt[:,6],label=str(pi_orc[6]))
# plt.plot(pi,Hopt[:,7],label=str(pi_orc[7]))
# plt.plot(pi,Hopt[:,8],label=str(pi_orc[8]))
plt.legend()
plt.show()