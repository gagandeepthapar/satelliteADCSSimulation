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

    def __init__(self, name, mass, xDim, yDim, zDim, delX, delY, delZ):
        self.name = name
        self.mass = mass
        self.xDim = xDim
        self.yDim = yDim
        self.zDim = zDim
        self.delX = delX
        self.delY = delY
        self.delZ = delZ
    
    def __repr__(self):
        return f"{self.name}: {self.mass} kg; {self.xDim} x {self.yDim} x {self.zDim} m"
    
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

    def __init__(self, bus, sensor, leftPanel, rightPanel):
        self.bus = bus
        self.sensor = sensor
        self.lPanel = leftPanel
        self.rPanel = rightPanel

        self.components = [self.bus, self.sensor, self.lPanel, self.rPanel]
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

class detumbleSat:
    # CM - [m] coordinates of center of mass of spacecraft relative to bus CM
    # Ibar - [kg*m^2] Inertia matrix of spacecraft relative to spacecraft CM

    CM = np.array([0,0,0]) # setting geo center to CM
    Ibar = np.array([[None, None, None], [None, None, None], [None, None, None]])


    def __init__(self, mass, xDim, yDim, zDim, wX, wY, wZ):
        # mass - [kg] total mass of spacecraft
        # xDim/yDim/zDim - [m] length of given dimension
        # wX/wY/wZ - [rad/s] initial rotation of spacecraft during detumble

        self.mass = mass
        self.xDim = xDim
        self.yDim = yDim
        self.zDim = zDim
        self.wX = wX
        self.wY = wY
        self.wZ = wZ

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

def main():
    # initialization

    # components
    bus = Component("BUS", 500, 2, 2, 2, 0, 0, 0)
    sensor = Component("SENSOR", 100, 0.25, 0.25, 1, 0, 0, 1.5)
    lPanel = Component("LEFTPANEL", 20, 2, 3, 0.05, 0, -2.5, 0)
    rPanel = Component("RIGHTPANEL", 20, 2, 3, 0.05, 0, 2.5, 0)

    # sat during detumble phase
    detSAT = detumbleSat(640, 2, 2, 2, -0.05, 0.03, 0.2)

    # sat during normal operations
    normalSAT = Spacecraft(bus, sensor, lPanel, rPanel)

    # Outputs
    print(f"Total mass of spacecraft (Detumble): {detSAT.mass} kg")
    print(f"Total mass of spacecraft (Normal): {normalSAT.mass} kg\n")
    
    print(f"Satellite center of mass relative to bus center of mass (Detumble): [{detSAT.CM[0]}, {detSAT.CM[1]}, {detSAT.CM[2]}] m")
    print(f"Satellite center of mass relative to bus center of mass (Normal): [{normalSAT.CM[0]}, {normalSAT.CM[1]}, {normalSAT.CM[2]:.3}] m\n")

    print(f"Inertia Matrix relative to center of mass (Detumble):\n{detSAT.Ibar} kg*m^2")
    print(f"Inertia Matrix relative to center of mass (Normal):\n{normalSAT.Ibar} kg*m^2")
    return

if __name__ == "__main__":
    main()
