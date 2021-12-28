from functools import update_wrapper
import numpy as np
import matplotlib.pyplot as plt
from _helperFuncs import *

def main():
    # given
    quatInitial_ECI = quaternion(np.array([0, 0, 0]), 1)
    angVel_ECI = np.array([0.001, -0.001, 0.002]).reshape(3,1)
    h = 53335.2
    ecc = 0
    inc = 98.43
    raan = 0
    arg = 0
    ta = 0
    orbitCOES = np.array([h, ecc, inc, raan, arg, ta])
    scOrbit = Orbit(coesList= orbitCOES)

    print(f"{scOrbit.vPath[0][0]}, {scOrbit.vPath[1][0]}, {scOrbit.vPath[2][0]}")

    


if __name__ == '__main__':
    main()
    # A = np.array([0,1,2]).reshape(3,1)
    # B = np.array([3,4,5]).reshape(3,1)
    # C = np.array([6,7,8]).reshape(3,1)
    # a = np.array([A,B,C]).reshape(3,3)
    # print(A)
    