import numpy as np
from numpy.linalg import solve
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

# GENERAL
# class creation for arbitrary central body
class CentralBody:
    def __init__(self, name, radius, Mu = None, mass = None):
        self.name = name
        self.rad = radius
        if(Mu != None):
            self.mu = Mu
        if(mass != None):
            G = 6.6743015*(10**-11)
            self.mu = G*mass

    def __repr__(self):
        return (f"{self.name}: Radius {self.rad} km; Mu {self.mu} km^3*s^-2")

# class creation for orbit around arbitrary body
class Orbit:
    def __init__(self, stateVector = None, coesList = None, centralBody = CentralBody("Earth", 6378, Mu = 398600)):
        if(coesList.all() == None):
            self.R = stateVector[0]
            self.V = stateVector[1]

            self.h, self.ecc, self.inc, self.raan, self.arg, self.ta, self.a = RV2coes(self.R, self.V, centralBody)

        else:
            self.h, self.ecc, self.inc, self.raan, self.arg, self.ta = coesList

            self.R, self.V = coes2RV(self.h, self.ecc, self.inc, self.raan, self.arg, self.ta)
        
        a = -centralBody.mu/(np.linalg.norm([self.V])**2 - (2*centralBody.mu/np.linalg.norm([self.R]))) # equation for semi-major axis
        self.period = 2*np.pi*(a**1.5)/(np.sqrt(centralBody.mu)) # equation for period of orbit

        self.body = centralBody
        self.createPath()
    
    def __repr__(self):
        ra = self.h ** 2 / self.body.mu * (1/(1 + self.ecc))
        rp = self.h ** 2 / self.body.mu * (1/(1 - self.ecc))
        
        a = ra + rp

        beta = np.arccos(-1*self.ecc)
        b = a * ((1-self.ecc**2)/(1 + self.ecc * np.cos(beta)))

        return f"{round(a, 3)} km x {round(b,3)} km orbit with inclination of {self.inc} deg"

    def createPath(self):
        sol = twoBodySol(self.R, self.V, self.body)
        self.tFrames = sol.t
        rX = sol.y[0]
        rY = sol.y[1]
        rZ = sol.y[2]
        vX = sol.y[3]
        vY = sol.y[4]
        vZ = sol.y[5]

        self.rPath = np.array([rX, rY, rZ])
        self.vPath = np.array([vX, vY, vZ])

        return np.array([self.rPath, self.vPath])

# class creation for arbitrary quaternion
class quaternion:
    def __init__(self, vector, scalar):
        self.E = vector
        self.n = scalar

    def __repr__(self):
        return f"\nE = [{self.E[0]}, {self.E[1]}, {self.E[2]}]; n = {self.n}\n"

# ODE solver for Kepler's Equations of Motion
def twoBodySol(R, V, centralBody = CentralBody("Earth", 6378, Mu = 398600)):

    def twoBody(t, state):
        R = np.array([state[0], state[1], state[2]])
        V = np.array([state[3], state[4], state[5]])

        rad = np.linalg.norm(R)
        
        ax = -1*centralBody.mu*R[0] / rad**3
        ay = -1*centralBody.mu*R[1] / rad**3
        az = -1*centralBody.mu*R[2] / rad**3

        # dState = np.array([V[0], V[1], V[2], ax, ay, az])
        # print(f"DSTATE SHAPE: {dState.shape}")

        dState = [V[0], V[1], V[2], ax, ay, az]

        return dState
    
    # state = np.array([R[0], R[1], R[2],
    #                     V[0], V[1], V[2]]).reshape(1,6)

    state = np.array([R[0], R[1], R[2], V[0], V[1], V[2]])

    # print(f"STATE: {state.shape}")

    a = -centralBody.mu/(np.linalg.norm([V[0], V[1], V[2]])**2 - (2*centralBody.mu/np.linalg.norm([R[0], R[1], R[2]]))) # equation for semi-major axis
    period = 2*np.pi*(a**1.5)/(np.sqrt(centralBody.mu)) # equation for period of orbit

                        
    return solve_ivp(twoBody, (0, np.ceil(period)), state, method = 'RK23', t_eval = np.linspace(0, np.ceil(period), 1001))

# returning state vector from classical orbital elements
def coes2RV(angularMomentum, eccentricity, inclination, raan, argumentOfPerigee, trueAnomaly, centralBody = CentralBody("Earth", 6378, Mu = 398600)):
    # Explanation of parameters
        # inc - [deg] inclination of orbit
        # raan - [deg] right ascenscion of ascending node
        # ecc - [~] eccentricity of orbit
        # arg - [deg] argument of perigee
        # h - [km3s2] specific angular momentum of orbit
        # ta - [deg] true anomaly of spacecraft on orbit
        # a - [km] semi major axis of orbit
        # M - [deg] mean anomaly (different from true anomaly)
    # renaming parameters for use
    h = angularMomentum
    ecc = eccentricity
    raan = np.deg2rad(raan)
    inc = np.deg2rad(inclination)
    arg = np.deg2rad(argumentOfPerigee)
    ta = np.deg2rad(trueAnomaly)

    cos = np.cos
    sin = np.sin

    periRConst = (h**2/centralBody.mu) * (1/(1 + ecc * cos(ta)))
    periVConst = centralBody.mu/h

    periRBar = periRConst * np.array([cos(ta), sin(ta), 0])
    periVBar = periVConst * np.array([-1*sin(ta), ecc + cos(ta), 0])

    periRBar = np.reshape(periRBar, (3,1))
    periVBar = np.reshape(periVBar, (3,1))

    Qa = rotZ(arg)
    
    Qb = rotX(inc)
                
    Qc = rotZ(raan)
    
    Q = np.matmul(np.matmul(Qa, Qb), Qc)
    Q = np.transpose(Q)

    R = np.matmul(Q, periRBar)
    V = np.matmul(Q, periVBar)

    R = R.reshape(1,3)
    V = V.reshape(1,3)

    return R[0], V[0]

# returning classical orbital elements from state vector
def RV2coes(R, V, centralBody = CentralBody("Earth", 6378, Mu = 398600)):
    # Explanation of parameters
        # R - [km, km, km] position vector of spacecraft at some time, t
        # V - [km/s, km/s, km/s] velocity vector of spacecraft at time t
    # Explanation of variables
        # h - angularMomentum [km3s-2]
        # inc - inclination [deg]
        # ecc - eccentricity
        # raan - right ascenscion of ascending node [deg]
        # arg - argumentOfPerigee [deg]
        # ta - trueAnomaly [deg]
        # a - semiMajorAxis [km]]

    rad = np.linalg.norm(R)
    vel = np.linalg.norm(V)
    
    radVel = np.dot(R, V)/rad
    
    hBar = np.cross(R, V)
    h = np.linalg.norm(hBar)

    inc = np.rad2deg(np.arccos(hBar[2]/h))

    nodeBar = np.cross(np.array([0, 0, 1]), hBar)
    node = np.linalg.norm(nodeBar)

    if node == 0:
        raan = 0
    else:
        raan = np.rad2deg(np.arccos(nodeBar[0]/node))
        if nodeBar[1] < 0:
            raan = 360-raan
    
    eccBar = 1/centralBody.mu * ((vel**2 - centralBody.mu/rad) * R - rad*radVel*V)
    ecc = np.linalg.norm(eccBar)

    if node == 0:
        arg = 0
    else:
        arg = np.rad2deg(np.arccos(np.dot(nodeBar/node , eccBar/ecc)))
        if eccBar[2] < 0:
            arg = 360-arg
    
    ta = np.rad2deg(np.arccos(np.dot(eccBar/ecc , R/rad)))
    if radVel < 0:
        ta = 360-ta
    
    rP = (h**2 / centralBody.mu) * (1/(1 + ecc))
    rA = (h**2 / centralBody.mu) * (1/(1 - ecc))
    
    a = 0.5*(rP + rA)

    # return [angularMomentum [km3s-2], eccentricity, inclination [deg], raan [deg], argumentOfPerigee [deg], trueAnomaly [deg], semiMajorAxis [km]]
    return np.array([h, ecc, inc, raan, arg, ta, a])

# 3x3 representation of cross of vector
def vectorCross(vector):
    rowA = np.array([0, -vector[2], vector[1]])
    rowB = np.array([vector[2], 0, -vector[0]])
    rowC = np.array([-vector[1], vector[0], 0])

    return np.array([rowA, rowB, rowC]).reshape(3,3)

# rotation matrix about x-axis for phi rad
def rotX(phi):
    cos = np.cos
    sin = np.sin

    rowA = np.array([1, 0, 0])
    rowB = np.array([0, cos(phi), sin(phi)])
    rowC = np.array([0, -sin(phi), cos(phi)])

    return np.array([rowA, rowB, rowC]).reshape(3,3)

# rotation matrix about y-axis for theta rad
def rotY(theta):
    cos = np.cos
    sin = np.sin

    rowA = np.array([cos(theta), 0, -sin(theta)])
    rowB = np.array([0, 1, 0])
    rowC = np.array([sin(theta), 0, cos(theta)])

    return np.array([rowA, rowB, rowC]).reshape(3,3)

# rotation matrix about z-axis for psi rad
def rotZ(psi):

    cos = np.cos
    sin = np.sin

    rowA = np.array([cos(psi), sin(psi), 0])
    rowB = np.array([-sin(psi), cos(psi), 0])
    rowC = np.array([0, 0, 1])

    return np.array([rowA, rowB, rowC]).reshape(3,3)

# rotation matrix between (non-standard) LVLH frame and ECI frame
def lvlh2eci(R, V):
    rad = np.linalg.norm(R)
    zz = -R/rad

    zzCross = vectorCross(R)
    ww = np.matmul(zzCross, V)

    yy = -ww/np.sqrt(np.matmul(ww, np.transpose(ww)))

    xx = np.cross(yy, zz)

    return np.transpose(np.array([xx[0], yy[0], zz[0],
                        xx[1], yy[1], zz[1],
                        xx[2], yy[2], zz[2]]).reshape(3,3))

# quaternion between two radius vectors
def rads2quat(r1, r2):
    Vc = np.cross(r1, r2)
    Vnc = Vc/np.linalg.norm(Vc)

    c = (np.dot(r1, r2))/(np.linalg.norm(r1) * np.linalg.norm(r2))
    s = -1*np.sqrt((1-c)/2)

    return quaternion(np.array([Vnc[0], Vnc[1], Vnc[2]]), -s)

# rotation matrix form of quaternion
def quat2matr(quat):
    q1 = quat.E[0]
    q2 = quat.E[1]
    q3 = quat.E[2]
    n = quat.n

    matrA = (2*(n**2)-1) * np.array([1, 0, 0, 0, 1, 0, 0, 0, 1]).reshape(3,3)

    tempB = 2*np.array([q1, q2, q3]).reshape(1,3)
    matrB = np.multiply(quat.E, tempB)

    matrC = 2*n*n*vectorCross(quat.E)

    return matrA + matrB - matrC

# quaternion form of rotation matrix
def matr2quat(matr):
    n = (np.sqrt(np.trace(matr) + 1))/2

    e1 = (matr[1][2] - matr[2][1])/(4*n)
    e2 = (matr[2][0] - matr[0][2])/(4*n)
    e3 = (matr[0][1] - matr[1][0])/(4*n)

    return quaternion(np.array([e1, e2, e3]), n)

if __name__ == '__main__':
    R = np.array([7136.6, 0, 0])
    V = np.array([0, -1.0956, 7.3927])

    matr = lvlh2eci(R, V)

    print(matr2quat(matr))