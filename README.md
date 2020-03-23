# Illumination conditions in pump-probe experiments

These tools were developed for evaluating the extent of excitation of light-absorbing molecules in a sample by pulsed laser illumination.

Both the python and Excel tool share the same basic functionality: given sample and pump illumination properties, the absorbed number of photons per molecule is calculated as a function of depth within the sample. Equal probabilities for single-and multi-photon absorption are assumed as a strong over-simplification. The python script additionally evaluates the fraction of molecules in single- and multiphoton regime.

Incident pump laser fluence is assumed to be the same for all area illuminated. This corresponds to the situation where the probe beam is very small, probing a cross sectional area of the excited sample over which the incident laser power did not change significantly. For example, this is the case in typical SFX experiments where the optical pump pulse is focused to ~100 µm diameter and the X-ray probe beam is focused to ~1 µm diameter.

## photon_calculation.py

### System requirements

This script is a standalone python script. It requires the following python packages:
- numpy
- matplotlib

It was tested with the following versions:
- python version 3.7.3
- numpy version 1.16.4
- matplotlib version 3.1.0 

### Instructions for use

Instructions how to run this script and what output to expect are listed at the beginning of the script.

Adapt the input parameters (sample properties, pump laser parameters) according to your system and in line with the provided instructions. Save as XXX.py file and run in python calling XXX.py. The script will output to the command shell / python environment, generate two figures which are both saved alongside an output .txt file summarising inputs and results.


## photon_calculation.xlsx

### System requirements

This is a standalone Excel file. It was developed and tested in Excel 2016 on a Windows machine (Windows 10)

### Instructions for use

Instructions how to use this file are listed at the beginning of the file.

Adapt the input parameters (sample properties, pump laser parameters) according to your system. The result is automatically updated and plotted.
