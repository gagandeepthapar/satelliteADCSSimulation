# **AERO - 421 Final Project**
**Spacecraft Attitude, Dynamics, and Control:** Simulation to determine and control a satellite's attitude in LEO.

## **Background**
AERO-421, or Spacecraft Attitude, Dynamics, and Controls, is a class taught at the California Polytechnic State University at San Luis Obispo (Cal Poly SLO) which serves as an...

> "Introduction to spacecraft attitude dynamics and control... [and] fundamentals of guidance and navigation systems... [with emphasis in] analysis and design of control systems for aerospace vehicles." - Cal Poly Aerospace Engineering Course Catalog

The final project in the course was to develop a simulation to determine a satellite's attitude in Low Earth Orbit (LEO), consider and model detumbling from a launch vehicle, consider and model disturbances due to external forces i.e., Solar Radiation Pressure (SRP), and to consider and model control via onboard reaction wheels. 

More information regarding the 6 different parts can be found below.

The initial project was developed in MATLAB, however, the project will be completely redeveloped in Python to showcase controls and software development skillsets.

## **Part 0: Context and Given Data**
NA

## **Part 1: Mass Properties**
Determine the mass and inertial properties of the spacecraft for both the detumble and the normal operations phases.

**Outputs:**
* Total mass of the spacecraft
* Center of mass relative to the spacecraft bus center of mass. The body frame will be located at the center of mass of the whole spacecraft
* Intertia matrix of the whole spacecraft about the center of mass of the spacecraft

## **Part 2: Torque Free Motion**
Model the torque free orbital and attitude motion of the spacecraft

**Outputs:**
Plots for...
* Euler angles and quaternions relating body to ECI reference frames
* Angular velocity of the spacecraft in body components for one orbit of the normal operations phase

## **Part 3: Detumble**
Simulate the motion of the satellite during the detumble phase. Assume fully modulated thrusters and use direct velocity feedback

**Outputs:**
Plots for... 
* Euler angles and quaternions relating body to ECI reference frames
* Angular velocity of the spacecraft in body components for the detumble phase
* Torque components in the body frame

## **Part 4: Disturbance Simulation**
Add the four disturbance models to the simulation:
* Atmospheric Drag
* Solar Pressure
* Gravity Gradient
* Earth Magnetic Field

Use the following model for the atmospheric density. Notice that **h** is the height above the Earth's surface in **kilometers** where **R_Earth** equals **6378km**

![Disturbance Model](./OutputFiles/0_DisturbanceModel.png)

Consider the simulation epoch to be March 20, 2021. Disregard any variations of the ECI representation of the sunlight direction during the simulation.

**Outputs:** Plots for...
* Euler angles and quaternions relating the body to the ECI reference frame
* Euler angles and quaternions relating the body to the LVLH reference frame
* Angular velocity of the spacecraft relative to the ECI frame expressed in body components
* Angular velocity of the spacecraft relative to the LVLH frame expressed in body components
* Torque components for atmospheric drag, solar radiation pressure, gravity gradient, and earth magnetic field

## **Part 5: Reaction Wheel Control**
Determine the control gains for a full state feedback 3-axis reaction wheel control system. Use the requirements of **Zeta = 0.65** and **t_s = 30 sec**

The positions of the 3 reaction wheels are **[1 0 0]**, **[0 1 0]**, and **[0 0 1]**. Each reaction wheel can be modeled as a simple cylinder with **radius** of **0.3 m** and a **height** of **0.02 m**

**Outputs:** Plots for...
* Euler angles and quaternions relating the body to ECI reference frame
* Euler angles and quaternions relating the body to LVLH reference frame
* Angular velocity of the spacecraft relative to the ECI reference frame expressed in body components
* Angular velocity of the spacecraft relative to LVLH frame expressed in body components
* Commanded moment from the determined control law
* Wheel speed of each reaction wheel

## **Part 6: Visualization**
Determine and animate the quanterions of the spacecraft, from ECI to body frame, for the duration of 1-5 revolutions. 

**Output**:
* Table of quaternion and time data
* Video or other animation file to show the revolution of the spacecraft