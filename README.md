# Orbit-Determination
Astrophysics Orbit Determination for 1988DZ4, including orbital elements calculator and ephemeris generator.
The purpose of this study is to determine the orbit of Mars-Crosser Asteroid 10737 (1988DZ4). Orbital elements are calculated using the Method of Gauss with observational data taken from 25 June 2020 to 22 July 2020 at the Central Washington University Observatory and the New Mexico Skies iTelescope Observatory. Using image processing softwares, the specific Right Ascension and Declination values for 1988DZ4 are determined, which are used as observational data for subsequent calculations. Specifically, the vectors for position and velocity are calculated from three pieces of observational data via the Method of Gauss. These vectors are then used to determine the orbital elements. After predicting preliminary values for the orbital elements, a Monte Carlo simulation is run in order to refine these results and provide an error on each of the orbital elements. It is these orbital elements which can be used to predict the position of the asteroid at any given time during its orbit. The results of the study contribute to a larger body of Near-Earth Asteroid data and provide the scientific community with a more accurate picture of the Earth’s celestial environment.

**Date: 07/15/2020**

## Contributors
**Main Contributor: Christopher Luey**


## Installation
Begin by downloading. Run LueyOD.py. Files will automatically unzip and directories will assemble. 

The input and output directories, respectively:

    /io/Observation
    /io/OrbitalElements
    
The following files will be used as tests:

    /io/Observation/2020test.txt
    /io/OrbitalElements/2020test.txt
    
Running the above test file should return the following orbital elements:
    
    Orbital Elements:
		a = 2.2061809479545387 AU
		e = 0.5709859541391656
		i = 18.773114736565933 deg
		ω = 280.0166781049967 deg
		Ω = 270.6062119686555 deg
		E = 2832.850394874371 deg
		M = 24.58713082789974 deg
		n = 0.005326177956890254 rad/day
		Last Perihelion Passage = 2457852.082865381 JD
		P = 3.229753615833997 yrs
    
The following libraries are required:
   
    numpy, astropy, pyephem

## Project Structure
Project structure as of 07/21:

    |-- io
        |-- Ephemeris
            |-- 2002QF15.txt
        |-- Images
            |-- 06-25
                |--corr_set3.fits
            |-- 07-03
                |--corr_set1.fits
                |--corr_set2.fits
                |--corr_set3.fits
            |-- 07-12
                |--corr_set1.fits
                |--corr_set2.fits
            |-- 07-22
                |--corr_set1.fits
                |--corr_set2.fits
        |-- MonteCarlo
            |-- MonteSimulation.txt
        |-- Observation
            |-- 1988DZ4.txt
            |-- 2020test.txt
        |-- OrbitalElements
            |-- 2020test.txt
    |-- Models
        |-- Distribution.jpeg
        |-- Model.py
        |-- Plot.py
    |-- odlib
        |-- ImageProcessing
            |-- __init__.py
            |-- AperaturePhotometry.py
            |-- Centroid.py
            |-- LSPR.py
        |-- odmath
            |-- NewtonRaphson.py
            |-- ScalarEquationLagrange.py
        |-- __init__.py
        |-- Asteroid.py
        |-- Constants.py
        |-- Conversion.py
        |-- Ephemeris.py
        |-- GaussMethod.py
        |-- OrbitalElements.py
    |-- LueyOD.py
    |-- MonteCarlo.py
    |-- README.md
    
## Usage
Run LueyOD.py
The code occasionally requests for user inputs. User input requests can be changed by inspecting LueyOD.py

When it asks for the input file, use this as a test file:

    io/Observation/LueyInput.txt
    
When it asks requests for the rows, use this for the 1st test case:
    
    0, 1, 3
    
However, you may choose to run all the combinations of test cases and find the one you're looking for.

## Project Status
Orbital determination of 1988DZ4 complete.