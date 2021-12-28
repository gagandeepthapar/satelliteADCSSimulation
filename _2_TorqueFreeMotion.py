from functools import update_wrapper
import numpy as np
import matplotlib.pyplot as plt
from _helperFuncs import *
import _1_MassProperties as p1

def main():
    # givens
    quatInitial = quaternion(np.array([0, 0, 0]), 1)
    angVel_ECI = np.array([0.001, -0.001, 0.002]).reshape(3,1)
    scIbar = np.array([812.0396, 545.3729, 627.7083])

    h = 53335.2
    ecc = 0
    inc = 98.43
    raan = 0
    arg = 0
    ta = 0
    orbitCOES = np.array([h, ecc, inc, raan, arg, ta])
    scOrbit = Orbit(coesList= orbitCOES)
    rInitial = np.array([scOrbit.rPath[0][0], scOrbit.rPath[1][0], scOrbit.rPath[2][0]])
    vInitial = np.array([scOrbit.vPath[0][0], scOrbit.vPath[1][0], scOrbit.vPath[2][0]])
    

    # computations:
        # QUATERNION
            # body to eci frame
    lvlh_eci_frame = lvlh2eci(rInitial, vInitial)
    body_lvlh_frame = quat2matr(quatInitial)
    body_eci_frame = np.matmul(body_lvlh_frame, lvlh_eci_frame)

            # initial quaternion in ECI frame
    initialQuat_ECI = matr2quat(body_eci_frame)

            # progression of quaternion across the period
    quatSol = quatProgressionSol(initialQuat_ECI, angVel_ECI, scIbar, scOrbit.period)

        # EULER ANGLES
            # initial euler angles
    phi = np.arctan2(body_eci_frame[1][2],body_eci_frame[2][2])
    theta = -np.arcsin(body_eci_frame[0][2])
    psi = np.arctan2(body_eci_frame[0][1],body_eci_frame[0][0])

            # progresion of euler angles across the period
    initialEuler = np.array([phi, theta, psi])
    eulerSol = eulerAngleSol(initialEuler, angVel_ECI, scIbar, scOrbit.period)

        # ANGULAR VELOCITY
    angVelSol = angVelProgressionSol(angVel_ECI, scIbar, scOrbit.period)

        # extract data to plot            
            # angle data
    angVel_time_prog = eulerSol.t
    phi_prog = eulerSol.y[0]
    theta_prog = eulerSol.y[1]
    psi_prog = eulerSol.y[2]

            # quaternion data
    quat_time_prog = quatSol.t
    e1_prog = quatSol.y[0]
    e2_prog = quatSol.y[1]
    e3_prog = quatSol.y[2]
    n_prog = quatSol.y[3]

            # angular velocity data
    angVel_time_prog = angVelSol.t
    wx_prog = angVelSol.y[0]
    wy_prog = angVelSol.y[1]
    wz_prog = angVelSol.y[2]

        # plot figures
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharey = False, sharex = True)

            # angle
    ax1.plot(angVel_time_prog, phi_prog, label = 'Phi')
    ax1.plot(angVel_time_prog, theta_prog, label = 'Theta')
    ax1.plot(angVel_time_prog, psi_prog, label = 'Psi')
    ax1.legend(loc = 'upper right')
    # ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Angle [rad]")
    ax1.grid(True)
    ax1.set_title("Euler Angle")

            # quaternions
    ax2.plot(quat_time_prog, e1_prog, label = 'q1')
    ax2.plot(quat_time_prog, e2_prog, label = "q2")
    ax2.plot(quat_time_prog, e3_prog, label = "q3")
    ax2.plot(quat_time_prog, n_prog, label = 'n')
    # ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Quaternion Parameter")
    ax2.legend(loc = 'upper right')
    ax2.grid(True)
    ax2.set_title("Quaternion Parameters")

            # angular velocity
    ax3.plot(angVel_time_prog, wx_prog, label = 'wX')
    ax3.plot(angVel_time_prog, wy_prog, label = 'wY')
    ax3.plot(angVel_time_prog, wz_prog, label = 'wZ')
    ax3.legend(loc = 'upper right')
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("Angular Velocity [rad/s]")
    ax3.grid(True)
    ax3.set_title("Angular Velocity")

    # ALTERNATIVE METHOD: SINGLE ODE
    torqueFreeProg = torqueFreeSol(initialQuat_ECI, initialEuler, angVel_ECI, scIbar, scOrbit.period)

    time = torqueFreeProg.t
    e1Dot = torqueFreeProg.y[0]
    e2Dot = torqueFreeProg.y[1]
    e3Dot = torqueFreeProg.y[2]
    nDot = torqueFreeProg.y[3]
    phiDot = torqueFreeProg.y[4]
    thetaDot = torqueFreeProg.y[5]
    psiDot = torqueFreeProg.y[6]
    wxDot = torqueFreeProg.y[7]
    wyDot = torqueFreeProg.y[8]
    wzDot = torqueFreeProg.y[9]

    f, (ax4, ax5, ax6) = plt.subplots(3,1, sharey = False, sharex = True)

    ax5.plot(time, e1Dot, label = 'q1')
    ax5.plot(time, e2Dot, label = 'q2')
    ax5.plot(time, e3Dot, label = 'q3')
    ax5.plot(time, nDot, label = 'n')
    ax5.legend(loc = 'upper right')
    ax5.grid(True)
    ax5.set_xlabel("Time [sec]")
    ax5.set_ylabel("Quaternion Parameters")
    ax5.set_title("Quaternion Parameters")

    ax4.plot(time, phiDot, label = 'Phi')
    ax4.plot(time, thetaDot, label = 'Theta')
    ax4.plot(time, psiDot, label = 'Psi')
    ax4.legend(loc = 'upper right')
    ax4.grid(True)
    ax4.set_xlabel("Time [sec]")
    ax4.set_ylabel("Angle [rad]")
    ax4.set_title("Euler Angle")
    
    ax6.plot(time, wxDot, label = 'wX')
    ax6.plot(time, wyDot, label = 'wY')
    ax6.plot(time, wzDot, label = 'wZ')
    ax6.legend(loc = 'upper right')
    ax6.grid(True)
    ax6.set_xlabel("Time [sec]")
    ax6.set_ylabel("Angular Velocity [rad/s]")
    ax6.set_title("Angular Velocity")

    plt.show()

if __name__ == '__main__':
    main()