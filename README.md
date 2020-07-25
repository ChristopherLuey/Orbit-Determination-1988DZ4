# Orbit-Determination
Astrophysics Orbit Determination for 1988DZ4, including orbital elements calculator and ephemeris generator.

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
            |-- monte.txt
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
Orbital determination of 1988DZ4 complete with Monte Carlo simulations