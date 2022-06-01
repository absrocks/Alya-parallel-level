#!/usr/bin/python

# pzt.py
# Implemented by Gerard Guillamet
#
# Last revision: 22 March 2018
#

from math import pi, cos, sin, exp, sqrt
from numpy import arange, array, transpose, linalg, dot, zeros, concatenate
import matplotlib.pyplot as plt

def laminate_bounds(sqt):
    """
    Laminate bounds at laminate mid-plane
    """
    z = 0.
    z_bounds = [-sum(sqt)/2.]
    for t_ply in sqt:
        z += t_ply
        z_bounds.append(z - sum(sqt)/2.)

    return z_bounds

def compliance_matrix(e11, e22, nu12, g12):
    """
    Compliance matrix S
    Material or local coordinates
    """
    s = array([[1./e11, -nu12/e11, 0.],
               [-nu12/e11, 1./e22, 0.],
               [0., 0., 1./g12]])

    return s

def transformed_matrix(sq):
    """
    Compute transformed matrix for each layer
    """
    t = []
    tg = []
    for theta_deg in sq:
        theta = theta_deg*pi/180.0
        t_ = array([[cos(theta)**2, sin(theta)**2,  2*cos(theta)*sin(theta)],
                    [sin(theta)**2, cos(theta)**2, -2*cos(theta)*sin(theta)],
                    [-cos(theta)*sin(theta), cos(theta)*sin(theta), cos(theta)**2-sin(theta)**2]])
        t.append(t_)
        tg.append(transpose(linalg.inv(t_)))

    return t, tg

def transformed_compliance_matrix(sq, e11, e22, nu12, g12):
    """
    Transformed Compliance matrix for each lamina [s_bar]
    Compliance matrix in global coordinates
    """
    # Compute Compliance matrix
    s = compliance_matrix(e11, e22, nu12, g12)

    # Compute transformation matrices t and t_gamma
    s_bar = []
    for i in range(len(sq)):
        t = transformed_matrix(sq)[0][i]
        tg = transformed_matrix(sq)[1][i]

        # Compute compliance transformed matrices
        s_bar.append(dot(dot(linalg.inv(tg), s), t))

    return s_bar

def transformed_stiffness_matrix(sq, e11, e22, nu12, g12):
    """
    Transformed Stiffness matrix for each lamina [Qb]
    Stiffness matrix in global coordinates
    """
    # Compute compliance transformed matrices
    s_bar = transformed_compliance_matrix(sq, e11, e22, nu12, g12)

    # Compute stiffness transformed matrices
    q_bar = []
    for i in range(len(sq)):
        q_bar.append(linalg.inv(s_bar[i]))

    return q_bar

def mat_abd(sq, sqt, e11, e22, nu12, g12):
    """
    Compute matrix [[A,B],[B,D]]
    """
    # Laminate Bounds
    z = laminate_bounds(sqt)

    # Compute transformed stiffness matrix
    q_bar = transformed_stiffness_matrix(sq, e11, e22, nu12, g12)

    a = zeros((3, 3))
    b = zeros((3, 3))
    d = zeros((3, 3))
    for i in range(len(sq)):
        for j in range(3):
            for k in range(3):
                a[j][k] += q_bar[i][j][k]*(z[i+1] - z[i])
                b[j][k] += q_bar[i][j][k]*(z[i+1]**2 - z[i]**2)*1./2
                d[j][k] += q_bar[i][j][k]*(z[i+1]**3 - z[i]**3)*1./3

    # Laminate matrix
    abd = concatenate((concatenate((a, b), axis=0),
                       concatenate((b, d), axis=0)), axis=1)

    return abd

def eng_constants_in_plane(sq, sqt, e11, e22, nu12, g12):
    """
    In-plane laminate engineering constants
    """
    # Total laminate thickness
    h = sum(sqt)

    # Compliance matrices
    abd = linalg.inv(mat_abd(sq, sqt, e11, e22, nu12, g12))

    # ENGINEERING CONSTANTS
    # Moduli
    ex = 1/(h*abd[0][0])
    ey = 1/(h*abd[1][1])
    gxy = 1/(h*abd[2][2])

    # Poisson ratios
    nuxy = -abd[0][1]/abd[0][0]
    nuyx = -abd[0][1]/abd[1][1]

    return ex, ey, gxy, nuxy, nuyx

def heaviside(x):
    """
    Step function Heaviside
    """
    if (x >= 0):
        result = 1
    elif (x == 0):
        result = 0.5
    else:
        result = 0

    return result

def radial_deformation(h_LMT, E_LMT, v_LMT, h_PZT, E_PZT, v_PZT, R, d31, Vin):
    """
    Equivalent radial deformation along the PZT top circumference
    """
    # Laminate parameters with the PZT
    A_tilde = (1. - v_LMT)/(1. - v_PZT)
    B_tilde = 2.*E_LMT*E_PZT*h_LMT*h_PZT*(1. - v_LMT)*(2.*h_LMT**2 + 3.*h_LMT*h_PZT + 2.*h_PZT**2)
    C_tilde = (E_LMT**2*h_LMT**4*(1. - v_PZT)**2 + E_PZT**2*h_PZT**4*(1. - v_LMT)**2)/(1. - v_PZT)
    D_tilde = (E_LMT*E_PZT*h_PZT*(4.*E_LMT*h_LMT**3 + 3.*E_LMT*h_LMT**2*h_PZT + A_tilde*E_PZT*h_PZT**3))/(B_tilde + C_tilde)
    # Equivalent Young Modulus
    E_tilde = A_tilde*D_tilde*E_PZT/E_LMT - E_PZT/(1. - v_PZT)

    # Equivalent radial deformation
    dr = R * d31/h_PZT*(E_tilde*(1. - v_PZT)/E_PZT + 1.)*Vin

    return dr

def exitation_signal(P, Np, fc, t):
    """
    Transient potetinal signal (input signal to the actuator)
    described by five cycle
    sinosuidal toneburst, modulated by Hanning
    """
    H_value1 = heaviside(t)
    H_value2 = heaviside(t - Np/fc)
    voltage_in = P*(H_value1 - H_value2)*(1. - cos(2.*pi*fc*t/Np))*sin(2.*pi*fc*t)

    return voltage_in

if __name__ == "__main__":
    """
    Write
    """
    print("###########################################")
    print(" ")
    print("-------------------------------------------")

    # Outfile parameters
    fileName = "d-t_Amp-50kHz.bc"

    # Laminate properties
    E_LMT = 100.5e9   # Equivalent laminate Modulus
    v_LMT = 0.36      # Equivalent laminate Poisson ratio
    h_LMT = 2.208e-3  # Laminate total thickness

    # Piezo-electric material properties
    E_PZT = 7.25e+10  # Young Modulus
    v_PZT = 0.31      # Poisson ratio
    h_PZT = 1.0e-3    # Total thickness
    R     = 5.0e-3    # Radius
    d31   = -170.0e-12  # Induced strain in direction 1 per unit electric field applied in direction 3

    # Excitation parameters
    Np = 5.           # Peak number
    P  = 55.          # Signal intensity (ICL has 20)
    fc = 50e+3        # Central frequency (Hz)
    T  = 1.1e-4       # Total time (s) 1.5e-4

    freq = 1e-6
    time_range = list(arange(0, T+freq, freq))
    f0 = open(fileName, 'w')
    f0.write('%8.0d\n' % (len(time_range)))
    t_list, v_in_list, d_r_list, h_list = [], [], [], []
    for time in time_range:
        v_in = exitation_signal(P, Np, fc, time)
        d_r  = radial_deformation(h_LMT, E_LMT, v_LMT, h_PZT, E_PZT, v_PZT, R, d31, v_in)
        t_list.append(time)
        v_in_list.append(v_in)
        d_r_list.append(-d_r)
        # Write file
        f0.write('%8.8e %8.8e %8.8e %8.8e\n' % (time, -d_r, 0.0, 0.0))
    f0.close()

    # Plot figure
    fig = plt.figure(num=fileName, figsize=(18, 7.0))
    a1 = fig.add_subplot(111)
    a1.plot(t_list, d_r_list, 'b-', label='d')
    a1.grid(True)
    a1.set_xlim([0,T])
    #a1.set_ylim([-1,1])
    a1.set_xlabel('Time(s)')
    a1.set_ylabel('Radial displacement (m)')

    #fig.savefig(job_name+'.png', bbox_inches='tight', dpi=300)

    plt.show()

    print("###########################################")

