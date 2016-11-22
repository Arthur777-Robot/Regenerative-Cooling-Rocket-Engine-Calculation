# Regenerative-Cooling-Rocket-Engine-Calculation
This program will calculate regenerative cooling rocket engine specification.

work in unix environment.

If Mingw and gnuplot are installed to your Windows and path are set, it will work also.

offline CEA from NASA needs to be install in order to use this program.

https://www.grc.nasa.gov/WWW/CEAWeb/ceaguiDownload-unix.htm

FCEA2 command must be pathed.

If you are Japanese, this website will describe how to install CEA to MACOS.

It is still under developement.


# How to use
1: Change the valuables in "valuables.h" to your desired parameter

2: compile the project using "make"

3: run the program by typing ./regene

#To Do List
able to spit out DXF file for the chamber geometry

able to select cooling channel type. now it is only for helical coil seen like in LR-101

able to select conical and bell nozzle

if there are any request, I will think about it.


#done list
able to calculate gas properties from CEA

able to calculate chamber geometry from CEA result

able to calculate regenerative cooling with model of helical cooling channel 

able to plot the properties to txt file

able to graph the properties using gnuplot


#License
The MIT License (MIT)

Copyright (c) 2016 Seiji Arthur Murakami

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

bove copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
