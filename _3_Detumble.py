import numpy as np
from matplotlib import pyplot as plt
from _helperFuncs import *

# PART SPECIFIC FUNCTIONS
# ODE for angular velocity, quaternion, euler angle, and thruster torque progression
def detumbleProgSol(initial_angVel, initial_quat, initial_eulerAng, initial_Thruster, gain, intertiaVector, period):
    
    def detumbleProg(t, state):
        # unpacking variables
        ix = intertiaVector[0]
        iy = intertiaVector[1]
        iz = intertiaVector[2]

        k = gain
        
        wx = state[0]
        wy = state[1]
        wz = state[2]

        e1 = state[3]
        e2 = state[4]
        e3 = state[5]
        n = state[6]

        phi = state[7]
        theta = state[8]
        psi = state[9]

        tx = state[10]
        ty = state[11]
        tz = state[12]

        # defining time derivative forms...
        # ANGULAR VELOCITY
        wxDot = (tx + (iy-iz)*wy*wz)/ix
        wyDot = (ty + (iz-ix)*wx*wz)/iy
        wzDot = (tz + (ix-iy)*wx*wy)/iz

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

        # THRUSTER TORQUES
        txDot = k * wxDot
        tyDot = k * wyDot
        tzDot = k * wzDot

        dState = np.array([wxDot, wyDot, wzDot, e1Dot, e2Dot, e3Dot, nDot, phiDot, thetaDot, psiDot, txDot, tyDot, tzDot])

        return dState

    # unpacking parameters
    wx = initial_angVel[0]
    wy = initial_angVel[1]
    wz = initial_angVel[2]

    e1 = initial_quat.E[0]
    e2 = initial_quat.E[1]
    e3 = initial_quat.E[2]
    n = initial_quat.n

    phi = initial_eulerAng[0]
    theta = initial_eulerAng[1]
    psi = initial_eulerAng[2]

    tx = initial_Thruster[0]
    ty = initial_Thruster[1]
    tz = initial_Thruster[2]

    state = np.array([wx, wy, wz, e1, e2, e3, n, phi, theta, psi, tx, ty, tz])

    return solve_ivp(detumbleProg, (0, 5*np.ceil(period)), state, method = 'RK23', t_eval = np.linspace(0, 5*np.ceil(period), 5*1000 + 1))

def main():
    # givens
    quatInitial = quaternion(np.array([0, 0, 0]), 1)
    scIbar = np.array([426.6667, 426.6667, 426.6667])
    detum_angVel = np.array([-0.05, 0.03, 0.2])

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
        # EULER ANGLES
            # initial euler angles
    phi = np.arctan2(body_eci_frame[1][2],body_eci_frame[2][2])
    theta = -np.arcsin(body_eci_frame[0][2])
    psi = np.arctan2(body_eci_frame[0][1],body_eci_frame[0][0])

            # progresion of euler angles across the period
    initialEuler = np.array([phi, theta, psi])
    gainVector = np.arange(-0.8, 0.1, 0.1)
    thrusterCombos = [gain * detum_angVel for gain in gainVector]
    
    solSet = [None] * np.size(gainVector)

    for i in range(np.size(gainVector)):
        tempSol = detumbleProgSol(detum_angVel, initialQuat_ECI, initialEuler, thrusterCombos[i], gainVector[i], scIbar, scOrbit.period)
        solSet[i] = tempSol.y
        print(f"TEMP SOL CALCULATED: {i+1}/{np.size(gainVector)}")

    time = tempSol.t

    print(f"Solution Size: {np.shape(solSet)}")

    figa, angVels = plt.subplots(3,3, sharey = True,sharex = True)
    figb, quats = plt.subplots(3,3, sharey = True, sharex = True)
    figc, euler = plt.subplots(3,3, sharey = True, sharex = True)
    figd, torq = plt.subplots(3,3, sharey = True, sharex = True)

    figa.canvas.manager.set_window_title("Angular Velocity")
    figb.canvas.manager.set_window_title("Quaternions")
    figc.canvas.manager.set_window_title("Euler Angles")
    figd.canvas.manager.set_window_title("Thruster Torques")

    for i in range(np.size(gainVector)):
        row, col = np.divmod(i, 3)

        # unpack variables
        wx, wy, wz, e1, e2, e3, n, phi, theta, psi, tx, ty, tz = solSet[i]

        # plot graphs
            # angular velocity
        angVels[row,col].plot(time, wx, label = 'wX')
        angVels[row,col].plot(time, wy, label = 'wY')
        angVels[row,col].plot(time, wz, label = 'wZ')

            # quaternions
        quats[row,col].plot(time, e1, label = 'q1')
        quats[row,col].plot(time, e2, label = 'q2')
        quats[row,col].plot(time, e3, label = 'q3')
        quats[row,col].plot(time, n, label = 'n')

            # euler angles
        euler[row,col].plot(time, phi, label = 'Phi')
        euler[row,col].plot(time, theta, label = 'Theta')
        euler[row,col].plot(time, psi, label = 'Psi')

            # thruster torques
        torq[row,col].plot(time, tx, label = 'tX')
        torq[row,col].plot(time, ty, label = 'tY')
        torq[row,col].plot(time, tz, label = 'tZ')
        

        # graph commentary
        plotTitle = f"GAIN: {np.round(gainVector[i], 1)}"
            # angular velocity
        angVels[row,col].set_title(plotTitle)
        angVels[row,col].grid(True)
        angVels[row,col].legend(loc = 'upper right')

            # quaternion parameters
        quats[row,col].set_title(plotTitle)
        quats[row,col].grid(True)
        quats[row,col].legend(loc = 'upper right')

            # euler angles
        euler[row,col].set_title(plotTitle)
        euler[row,col].grid(True)
        euler[row,col].legend(loc = 'upper right')

            # thruster torques
        torq[row,col].set_title(plotTitle)
        torq[row,col].grid(True)
        torq[row,col].legend(loc = 'upper right')

        if(row > 1):    # only set bottom row x-axis
            angVels[row,col].set_xlabel("Time [sec]")
            quats[row,col].set_xlabel("Time [sec]")
            euler[row,col].set_xlabel("Time [sec]")
            torq[row,col].set_xlabel("Time [sec]")

        if(col == 0):   # only set first col y-axis
            angVels[row,col].set_ylabel("Angular Velocity [rad/s]")
            quats[row,col].set_ylabel("Quaternion Parameters")
            euler[row,col].set_ylabel("Euler Angle [rad]")
            torq[row,col].set_ylabel("Thruster Torque [N-m]")

    plt.show()

    return
    
if __name__ == '__main__':
    main()