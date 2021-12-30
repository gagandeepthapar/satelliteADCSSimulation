import numpy as np
import matplotlib.pyplot as plt
from _helperFuncs import *
from scipy.integrate import solve_ivp

# PART SPECIFIC FUNCTIONS
# ODE to solve for progression of quaternion over a single period
def quatProgressionSol(quaternion, angularVelocity, inertiaVector, period):
    def quatProgression(t, state):
        ix = inertiaVector[0]
        iy = inertiaVector[1]
        iz = inertiaVector[2]

        e1 = state[0]
        e2 = state[1]
        e3 = state[2]
        n = state[3]
        wx = state[4]
        wy = state[5]
        wz = state[6]

        e1Dot = ((n*wx) - (e3*wy) + (e2*wz))/2
        e2Dot = ((e3*wx) + (n*wy) - (e1*wz))/2
        e3Dot = ((-e2*wx) + (e1*wy) + (n*wz))/2

        E = np.array([e1, e2, e3])
        W = np.array([wx, wy, wz])
        nDot = np.matmul(E, np.transpose(W)) * -0.5

        wxDot = (iy-iz)*wy*wz/ix
        wyDot = (iz-ix)*wx*wz/iy
        wzDot = (ix-iy)*wx*wy/iz

        dState = np.array([e1Dot, e2Dot, e3Dot, nDot, wxDot, wyDot, wzDot])

        return dState

    e1 = quaternion.E[0]
    e2 = quaternion.E[1]
    e3 = quaternion.E[2]
    n = quaternion.n
    wX = angularVelocity[0][0]
    wY = angularVelocity[1][0]
    wZ = angularVelocity[2][0]

    state = np.array([e1, e2, e3, n, wX, wY, wZ])

    return solve_ivp(quatProgression, (0, np.ceil(period)), state, method = 'RK23', t_eval = np.linspace(0, np.ceil(period), 1001))

# ODE to solve for progression of Euler Angles over a single period
def eulerAngleSol(initialEuler, angularVelocity, inertiaVector, period):
    def eulerAngle(t, state):
        ix = inertiaVector[0]
        iy = inertiaVector[1]
        iz = inertiaVector[2]

        phi = state[0]
        theta = state[1]
        psi = state[2]
        wx = state[3]
        wy = state[4]
        wz = state[5]

        cos = np.cos
        sin = np.sin

        angVelVec = np.transpose(np.array([wx, wy, wz]))

        phiMatr = np.array([cos(theta), sin(phi) * sin(theta), cos(phi) * sin(theta)])
        phiDot = (1/cos(theta)) * np.matmul(phiMatr, angVelVec)

        thetaMatr = np.array([0, cos(phi)*cos(theta), -sin(phi)*cos(theta)])
        thetaDot = (1/cos(theta)) * np.matmul(thetaMatr, angVelVec)

        psiMatr = np.array([0, sin(phi), cos(phi)])
        psiDot = (1/cos(theta)) * np.matmul(psiMatr, angVelVec)

        wxDot = (iy-iz)*wy*wz/ix
        wyDot = (iz-ix)*wx*wz/iy
        wzDot = (ix-iy)*wx*wy/iz

        dState = np.array([phiDot, thetaDot, psiDot, wxDot, wyDot, wzDot])
        
        return dState

    phi = initialEuler[0]
    theta = initialEuler[1]
    psi = initialEuler[2]
    wx = angularVelocity[0][0]
    wy = angularVelocity[1][0]
    wz = angularVelocity[2][0]

    state = np.array([phi, theta, psi, wx, wy, wz])

    return solve_ivp(eulerAngle, (0, np.ceil(period)), state, method = 'RK23', t_eval = np.linspace(0, np.ceil(period), 1001))

# ODE to solve for progression of angular velocity components over single period
def angVelProgressionSol(initialVelocity, inertiaVector, period):
    def angVelProgression(t, state):
        ix = inertiaVector[0]
        iy = inertiaVector[1]
        iz = inertiaVector[2]

        wx = state[0]
        wy = state[1]
        wz = state[2]

        wxDot = (iy-iz)*wy*wz/ix
        wyDot = (iz-ix)*wx*wz/iy
        wzDot = (ix-iy)*wx*wy/iz

        dState = np.array([wxDot, wyDot, wzDot])

        return dState

    wx = initialVelocity[0][0]
    wy = initialVelocity[1][0]
    wz = initialVelocity[2][0]

    state = np.array([wx, wy, wz])

    return solve_ivp(angVelProgression, (0, np.ceil(period)), state, method = 'RK23', t_eval = np.linspace(0, np.ceil(period), 1001))

# ODE to solve for quaternion, euler angles, AND angular velocity components over single period
def torqueFreeSol(initial_quaternion, initial_euler_angles, initial_angular_velocity, inertiaVector, period):
    def torqueFree(t, state):
        ix = inertiaVector[0]
        iy = inertiaVector[1]
        iz = inertiaVector[2]

        e1 = state[0]
        e2 = state[1]
        e3 = state[2]
        n = state[3]

        phi = state[4]
        theta = state[5]
        psi = state[6]

        wx = state[7]
        wy = state[8]
        wz = state[9]

        # QUATERNION
        e1Dot = ((n*wx) - (e3*wy) + (e2*wz))/2
        e2Dot = ((e3*wx) + (n*wy) - (e1*wz))/2
        e3Dot = ((-e2*wx) + (e1*wy) + (n*wz))/2

        E = np.array([e1, e2, e3])
        W = np.array([wx, wy, wz])
        nDot = np.matmul(E, np.transpose(W)) * -0.5

        # EULER ANGLES
        cos = np.cos
        sin = np.sin

        angVelVec = np.transpose(np.array([wx, wy, wz]))

        phiMatr = np.array([cos(theta), sin(phi) * sin(theta), cos(phi) * sin(theta)])
        phiDot = (1/cos(theta)) * np.matmul(phiMatr, angVelVec)

        thetaMatr = np.array([0, cos(phi)*cos(theta), -sin(phi)*cos(theta)])
        thetaDot = (1/cos(theta)) * np.matmul(thetaMatr, angVelVec)

        psiMatr = np.array([0, sin(phi), cos(phi)])
        psiDot = (1/cos(theta)) * np.matmul(psiMatr, angVelVec)

        # ANGULAR VELOCITY
        wxDot = (iy-iz)*wy*wz/ix
        wyDot = (iz-ix)*wx*wz/iy
        wzDot = (ix-iy)*wx*wy/iz

        dState = np.array([e1Dot, e2Dot, e3Dot, nDot, phiDot, thetaDot, psiDot, wxDot, wyDot, wzDot])

        return dState

    e1 = initial_quaternion.E[0]
    e2 = initial_quaternion.E[1]
    e3 = initial_quaternion.E[2]
    n = initial_quaternion.n

    phi = initial_euler_angles[0]
    theta = initial_euler_angles[1]
    psi = initial_euler_angles[2]

    wx = initial_angular_velocity[0][0]
    wy = initial_angular_velocity[1][0]
    wz = initial_angular_velocity[2][0]

    state = np.array([e1, e2, e3, n, phi, theta, psi, wx, wy, wz])

    return solve_ivp(torqueFree, (0, np.ceil(period)), state, method = 'RK23', t_eval = np.linspace(0, np.ceil(period), 1001))

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