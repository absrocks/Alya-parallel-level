#
# Ganser model, different coefficients used by Alya to be plotted by Gnuplot
# Cd
# Cd*Re
# d(Cd*Re)/dRe
# phi = 1/Cd * d(Cd*Re)/dRe
#
reset
psi                = 1.0
k1                 = 3.0 / ( 1.0 + 2.0/sqrt(psi) )
k2                 = 10.0**( 1.84148 * (-log10(psi))**0.5743 )
Re(x)              = x>10**5 ? 10**5 : x
CdRe_ganser(x)     = 24.0/k1 * ( 1.0 + 0.1118*(Re(x)*k1*k2)**0.6567 ) + 0.4305*Re(x)*Re(x)*k2/ ( Re(x) + 3305.0/(k1*k2) )
Cd_ganser(x)       = 1.0/Re(x)*(24.0/k1 * ( 1.0 + 0.1118*(Re(x)*k1*k2)**0.6567 ) + 0.4305*Re(x)*Re(x)*k2/ ( Re(x) + 3305.0/(k1*k2) ))
dCdRedRe_ganser(x) = 24.0/k1 * ( 0.1118*0.6567*(k1*k2)**0.6567*Re(x)**(1.0-0.6567) ) \
                     + 0.4305 * Re(x) * 2.0 * k2 / ( Re(x) + 3305.0/(k1*k2) ) \
                     - 0.4305 * Re(x) * Re(x) * k2 / ( Re(x) + 3305.0/(k1*k2) ) ** 2
phi_ganser(x)      = dCdRedRe_ganser(x)/Cd_ganser(x)                     

Re2(x)              = x
CdRe_ganser2(x)     = 24.0/k1 * ( 1.0 + 0.1118*(Re2(x)*k1*k2)**0.6567 ) + 0.4305*Re2(x)*Re2(x)*k2/ ( Re2(x) + 3305.0/(k1*k2) )
Cd_ganser2(x)       = 1.0/Re2(x)*(24.0/k1 * ( 1.0 + 0.1118*(Re2(x)*k1*k2)**0.6567 ) + 0.4305*Re2(x)*Re2(x)*k2/ ( Re2(x) + 3305.0/(k1*k2) ))
dCdRedRe_ganser2(x) = 24.0/k1 * ( 0.1118*0.6567*(k1*k2)**0.6567*Re2(x)**(1.0-0.6567) ) \
                     + 0.4305 * Re2(x) * 2.0 * k2 / ( Re2(x) + 3305.0/(k1*k2) ) \
                     - 0.4305 * Re2(x) * Re2(x) * k2 / ( Re2(x) + 3305.0/(k1*k2) ) ** 2
phi_ganser2(x)      = dCdRedRe_ganser2(x)/Cd_ganser2(x)                     


set format x "10^{%2T}
set format y "10^{%2T}
set xlabel 'Re'
set xrange[1.0e-3:1e6]
set log x
set log y
plot phi_ganser(x)
