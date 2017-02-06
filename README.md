# FCS-calc-ImageJ-plugin
An ImageJ plugin for calculating correlation functions for FCS (Fluorescence Correlation Spectroscopy)

Installation:
1. Make sure that ImageJ is installed on your computer.
2. Download FCS_calc.java and Gmn.java source files and copy them to the ImageJ plugins folder.
3. Compile Gmn by running it from Plugins->Compile and Run... in ImageJ.
4. Compile FCS_calc by running it from Plugins->Compile and Run... in ImageJ.
5. The correlator is ready for use. To launch your calculation, select “FCS calc” under Plugins.

FCS_calc.java:
FCS_calc is a plugin for ImageJ that provides user interface for a software correlator realized in Gmn.java. “Time base” is given in units of time resolution of your experiment and must always be an integer; the value of 1 results in the highest time resolution of the correlation function and the slowest calculation. “N of cascades” and “N of points per cascade” are parameters used for constructing the semi-logarithmic X axis. The calculation can be run in either 'auto' (AxA) or 'cross' (AxA, BxB, AxB, and BxA) mode. Multiple files can be selected for sequential calculation; for data types that store channels A and B separately, select channel A, then channel B files.

Gmn.java:
The Gmn class has an initialization method, several utility methods, methods for converting input data files into a universal format with inter-photon arrival times and photon weights, and two methods for calculating auto and cross correlations. The correlation algorithm realized by Gmn is based on

Yang et al. (2009) Real-time data acquisition incorporating high-speed software correlator for single molecule spectroscopy. J of Microscopy 234, 302-310.

Symmetric normalization and the multi-tau scheme of time coarsening are described in

Schatzel et al. (1988) Photon correlation measurement at large lag times: improving statistical accuracy. J of Modern Optics 35, 711-718.
Wahl et al. (2003) Fast calculation of fluorescence correlation data with asynchronous time-correlated single-photon counting. Optics Express 11, 3583-3591.


Supported data types:
Presently three data types are supported: Flex 8-bit, Confocor2 and Confocor3 raw data files. If your data type is not supported, please email me with a description of the data format and an example file, and I will incorporate your data type into FCS_calc.

Future additions:
1. Display of intensity traces.
2. Calculation of higher order auto-correlations.
3. Support for other data formats.

Edited on 01/28/16.
