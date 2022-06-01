from math import pi, sqrt
# Units: SI (mm)

r_i    = 8.              # Indenter radius (mm)
m_i    = 2.0e-3          # Mass (Tonne)
h_i    = r_i*1.          # Indenter total height (mm)
energy = 10*1000.        # Energy (N*mm)
g      = 9810.           # Gravity (mm/s^2)
h_o    = energy/m_i/g    # Total height at weight drop (mm)
h_user = 0.01            # User-defined height (mm)

# Indenter drop characteristics (a=ctt)
# t_o:    v_o = 0 (free fall)
# t_user: Total time imposedVelocity imposed as initial condition (v_mid)
# t_f:    v_fin
# t_simu: Simulation time before contact

t_o = 0.0
v_o = 0.0
rho     = m_i/(pi*r_i**2*(h_o-r_i) + 2/3.*pi*r_i**3)
v_f     = sqrt(energy*2/m_i)
t_mid   = sqrt((h_o - h_user)*2./g)
t_f     = (v_f - v_o)/g
t_sim   = t_f - t_mid
v_user  = v_f - g*t_sim

print('')
print('USER INPUTS')
print('Impactor density, rho (Tonne/mm^2)',rho)
print('Initial Velocity, v_user (mm/s):', v_user)
print('Impactor height, h_user: (mm):',h_user)
print('')
print('OTHER CALCULATIONS')
print('Time before contact (assuming free fall), t_f (s):',t_f)
print('Time before contact (assuming an initial velocity, v_user), t_sim (s):',t_sim)
print('Final Velocity, v_f (mm/s):',v_f)
print('Impactor height (assuming free fall), h_o (mm):',h_o)
print('')
