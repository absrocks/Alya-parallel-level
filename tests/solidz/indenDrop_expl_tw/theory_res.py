from math import pi, sqrt
# Units: SI (mm)

r_i    = 8.              # Indenter radius (mm)
m_i    = 0.00341         # Mass (Tonne)
h_i    = r_i*3.          # Indenter total height (mm)
energy = 35*1000.        # Energy (N*mm)
g      = 9810.           # Gravity (mm/s^2)
h_ini  = energy/m_i/g    # Total height at weight drop (mm)
h_mid  = 4.53            # User-defined height (mm)

# Indenter drop characteristics (a=ctt)
# t_ini: v_ini = 0
# t_mid: Velocity imposed as initial condition (v_mid)
# t_fin: v_fin
# t_sim: Simulation time

t_ini = 0.0
v_ini = 0.0
rho   = m_i/(pi*r_i**2*(h_i-r_i) + 2/3.*pi*r_i**3)
v_fin = sqrt(energy*2/m_i)
t_mid = sqrt((h_ini - h_mid)*2./g)
t_fin = (v_fin - v_ini)/g
h_fin = h_ini - 0.5*g*t_fin**2

# Time and velocity used in the simulation (inputs)
t_sim = t_fin - t_mid
v_mid = v_fin - g*t_sim

print('Density:')
print('rho:',rho)

print('Times:')
print('t_ini:',t_ini)
print('t_mid:',t_mid)
print('t_fin:',t_fin)
print('t_sim:',t_sim)

print('Velocities:')
print('v_ini:',v_ini)
print('v_mid:',v_mid)
print('v_fin:',v_fin)

print('Heights:')
print('h_ini:',h_ini)
print('h_mid:',h_mid)
print('h_fin:',h_fin)
