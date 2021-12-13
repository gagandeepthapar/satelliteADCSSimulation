import numpy as np
from matplotlib import pyplot as plt

# Motivation/Thought process:
# Center of the spacecraft body will act as coordinate axes origin
# Will use the spacecraft schematic (see README) for axis reference

class Component:
    # name - name of component
    # mass - [kg] mass of component
    # xDim/yDim/zDim - [m] length of given dimension
    # delX/delY/delZ - [m] component of vector pointing from bus CM to local CM

    def __init__(self, name, mass, localCM, vect2CM):
        self.name = name
        self.mass = mass
        self.CM = localCM
        self.vect2CM = vect2CM

        self.delX = self.vect2CM[0]
        self.delY = self.vect2CM[1]
        self.delZ = self.vect2CM[2]
    
    def __repr__(self):
        return f"{self.name}: {self.mass} kg; CM: {self.localCM[0]} x {self.localCM[1]} x {self.localCM[2]} m"
    
class Spacecraft:
    # mass - [kg] total mass of spacecraft
    # CM - [m] coordinates of center of mass of spacecraft relative to bus CM
    # Ibar - [kg*m^2] Inertia matrix of spacecraft relative to spacecraft CM

    # bus - spacecraft main bus
    # sensor - spacecraft sensor
    # lPanel - spacecraft left solar panel
    # rPanel - spacecraft right solar panel

    mass = None;
    CM = np.array([None, None, None])
    Ibar = np.array([[None, None, None], [None, None, None], [None, None, None]])

    def __init__(self, componentList):

        self.components = componentList
        self.setMass() # compute mass at init time
        self.setCM()  # compute CM at init time
        self.setIbar() # compute Inertia matrix at init time

    def __repr__(self):
        return f"Mass: {self.mass()} kg; CM: {self.CMs()} m"
    
    def setMass(self):
        mass = 0
        for comp in self.components:
            mass += comp.mass
        self.mass = mass
        return self.mass

    def setCM(self):
        # recall Center of Mass = sum(mass * distance)/total mass
        CMx = sum([comp.mass*comp.delX for comp in self.components])/self.mass
        CMy = sum([comp.mass*comp.delY for comp in self.components])/self.mass
        CMz = sum([comp.mass*comp.delZ for comp in self.components])/self.mass

        self.CM = np.array([CMx, CMy, CMz])
        return self.CM
    
    def setIbar(self):
        return
    
    def plotSat(self):
        fig = plt.figure()
        ax = plt.axes(projection = '3d')
        ax.set_xlabel("X [m]")
        ax.set_ylabel("Y [m]")
        ax.set_zlabel("Z [m]")
        ax.set_title("Vectors Pointing from Spacecraft CM to Component CM (Normal Operations)")

        compNames = []
        for comp in self.components:
            vect = cmVect(self.CM, comp.vect2CM)
            plt.plot(vect[0], vect[1], vect[2])
            compNames.append(comp.name)

        ax.legend(compNames)
        plt.show()

        return

class detumbleSat:
    # CM - [m] coordinates of center of mass of spacecraft relative to bus CM
    # Ibar - [kg*m^2] Inertia matrix of spacecraft relative to spacecraft CM

    CM = np.array([0,0,0]) # setting geo center to CM
    Ibar = np.array([[None, None, None], [None, None, None], [None, None, None]])


    def __init__(self, mass, dimensions, vect2CM, angularVelocity):
        # mass - [kg] total mass of spacecraft
        # xDim/yDim/zDim - [m] length of given dimension
        # wX/wY/wZ - [rad/s] initial rotation of spacecraft during detumble

        self.mass = mass
        self.dims = dimensions
        self.CM = vect2CM
        self.angVel = angularVelocity

        self.wX = self.angVel[0]
        self.wY = self.angVel[1]
        self.wZ = self.angVel[2]
        self.xDim = self.dims[0]
        self.yDim = self.dims[1]
        self.zDim = self.dims[2]

        self.setIbar()
    
    def __repr__(self):
        return f"Mass: {self.mass} kg; {self.xDim} x {self.yDim} x {self.zDim} m"
    
    def setIbar(self):
        const = self.mass/12;
        # in relation to the CM...
        Ix = (self.yDim ** 2) + (self.zDim ** 2)
        Iy = (self.xDim ** 2) + (self.zDim ** 2)
        Iz = (self.xDim ** 2) + (self.yDim ** 2)

        Matr = np.array([Ix, 0, 0, 0, Iy, 0, 0, 0, Iz])
        self.Ibar = const * (np.reshape(Matr, (3,3)))

        return self.Ibar

def cmVect(CM, localCM):
    x = np.array([CM[0],localCM[0] - CM[0]])
    y = np.array([CM[1],localCM[1] - CM[1]])
    z = np.array([CM[2],localCM[2] - CM[2]])
    return [x, y, z]

def main():
    # initialization

    # components
    bus = Component("BUS", 500, np.array([2, 2, 2]), np.array([0, 0, 0]))
    sensor = Component("SENSOR", 100, np.array([0.25, 0.25, 1]), np.array([0, 0, 1.5]))
    lPanel = Component("LEFT PANEL", 20, np.array([2, 3, 0.05]), np.array([0, -2.5, 0]))
    rPanel = Component("RIGHT PANEL", 20, np.array([2, 3, 0.05]), np.array([0, 2.5, 0]))
    componentList = [bus, sensor, lPanel, rPanel]

    # sat during detumble phase
    detSAT = detumbleSat(640, np.array([2, 2, 2]), np.array([0, 0, 0]), np.array([-0.05, 0.03, 0.2]))

    # sat during normal operations
    normalSAT = Spacecraft(componentList)

    # Outputs
    outFile = open("OutputFiles/1_MassPropertiesOutput.txt", "a")
    outFile.truncate(0)
    print(f"MASS OF SPACECRAFT:", file = outFile)
    print(f"Total mass of spacecraft (Detumble): {detSAT.mass} kg", file = outFile)
    print(f"Total mass of spacecraft (Normal): {normalSAT.mass} kg\n", file = outFile)
    
    print(f"CENTER OF MASS OF SPACECRAFT:", file = outFile)
    print(f"Satellite center of mass relative to bus center of mass (Detumble): [{detSAT.CM[0]}, {detSAT.CM[1]}, {detSAT.CM[2]}] m", file = outFile)
    print(f"Satellite center of mass relative to bus center of mass (Normal): [{normalSAT.CM[0]}, {normalSAT.CM[1]}, {normalSAT.CM[2]:.3}] m\n", file = outFile)

    print(f"INERTIA MATRIX OF SPACECRAFT:", file = outFile)
    print(f"Inertia Matrix relative to center of mass (Detumble):\n{detSAT.Ibar} kg*m^2", file = outFile)
    print(f"Inertia Matrix relative to center of mass (Normal):\n{normalSAT.Ibar} kg*m^2", file = outFile)

    print(f"\nNotice: the Z-axis for the spacecraft is positive in the direction of the sensor", file = outFile)

    normalSAT.plotSat()

    return

if __name__ == "__main__":
    # main()
    # plotSat()
    main()
