# Orbit-Determination
Astrophysics Orbit Determination for 1988DZ4, including orbital elements calculator and ephemeris generator.

**Date: 07/15/2020**

## Contributors
**Main Contributor: Christopher Luey**


## Installation
Begin by downloading and unzip the files. 

The input and output directories, respectively:

    /io/Observation
    /io/OrbitalElements
    
The following files will be used as tests:

    /io/Observation/2020test.txt
    /io/OrbitalElements/2020test.txt
    
The following libraries are required:
   
    numpy, astropy, pyephem

## Project Structure
Project structure as of 07/21:

    |-- io
        |-- Ephemeris
            |-- 2002QF15.txt
        |-- Observation
            |-- 1988DZ4.txt
            |-- 2020test.txt
        |-- OrbitalElements
            |-- 2020test.txt
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
    |-- README.md
    
## Usage
Run LueyOD.py
The code occasionally requests for user inputs. User input requests can be changed by inspecting LueyOD.py


## Project Status
1. OD component of project is complete. 
2. Implementation of Monte Carlo and differential correction are expected to be completed by 07/22. 
3. Testing OD with 1988DZ4 observation data is ready. Waiting for Monte Carlo.