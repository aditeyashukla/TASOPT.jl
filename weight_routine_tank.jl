include("tanksize.jl")
include("tankWmech.jl")
include("tankWthermal.jl")


v_f = 67.3 #specific volume as saturated liquid
v_g = 2.5 #specific volume as saturated vapor
x = 0.03 #quality (suggested by energies)

rhoFuel = 0.96 * (1 / ((1 / v_f) + x * ((1 / v_g) - (1 / v_f)))) #96% due to pvent = 2atm in energies literature

gee = 9.81

deltap = 150000#150000#101325  #based on 0.5 bar cabin pressure, 2 or 1.5 bar LH2 pressure (Verstraete / Brewer)

Rfuse = 2.07 #radius of fuselage of A320

dRfuse = 0 #extended portion of fuselage at bottom (take as zero due to cylindrical tank)

hconvgas = 0 #no purged gas in foam insulated tank

h_LH2 = 750 #based on elsevier pub film boiling heat transfer properties of liq hydrogen

Tfuel = 20 #20 K

Tair = 293 #cabin heated temperature

h_v = 447000 #Engineering toolbox, Hydrogen - Thermophysical properties

t_cond = [0.1, 1.524e-5, 0.1, 1.524e-5, 1.57e-2] #assumed from energies

k = [20e-3, 20e-3, 20e-3, 20e-3, 20e-3] #foam conductivities

#http://www.hysafe.org/download/1196/BRHS_Chap1_V1p2.pdf
hconvair = 25 #from sciencedirect.com https://www.sciencedirect.com/topics/engineering/convection-heat-transfer-coefficient

time_flight = 18000 #seconds, 5 hr flight

fstring = 0.34

ffadd = 0.2

wfb = 0

nfweb = 0

sigskin = 172.4e6 #AL 2219 Brewer / energies stress for operating conditions (290e6 ultimate operatoin)

rho_insul = [35.24, 14764, 35.24, 14764, 83] #energies

rhoskin = 2825 #Al 2219 / energies
AR = 2 #ellipsoidal head aspect ratio
#Wfuel = 8000
Wfuel = 800 * 9.81 #8000 * 9.81

threshold_percent = 0.1

#r_tank = 1.722

mode = 2

clearance_fuse = 0.13 #m , Brewer suggests 5 in clearance

Wtank_total, thickness_insul, lshell, mdot_boiloff, Vfuel, Wfuel_tot, m_boiloff, tskin, t_head, Rtank, Whead, Wcyl, Winsul_sum, Winsul = tanksize(gee, rhoFuel, deltap,
                      Rfuse, dRfuse, hconvgas, h_LH2, Tfuel, Tair,
                      h_v, t_cond, k, hconvair, time_flight, fstring,ffadd,
                      wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, mode, clearance_fuse, AR)
#comments
#gee is gravitational acceleration, 9.81 m/s^2
#rhoFuel is calculated based on 3% quality (energies paper)
#deltap is based on pressure difference across tank wall,



#remove when running (only for test)
#m_boiloff = 0.001 * Wfuel
#thickness_insul = sum(t_cond)
#lshell = 2.38
#r_tank = 1.722
