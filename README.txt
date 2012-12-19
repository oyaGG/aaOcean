aaOcean v2.5
Author: Amaan Akram 
www.amaanakram.com
aaOcean is free software and can be redistributed and modified under the terms of the 
GNU General Public License (Version 3) as provided by the Free Software Foundation.
GNU General Public License http://www.gnu.org/licenses/gpl.html
Please note that this aaOcean core code is not the one used to generate the binaries that are downloadable from my website

Some key features
-Multi-threaded via openmp
-Has gone through a few levels of code optimization via profiling
-The shape of the ocean does not change between resolution changes. You can change resolution at any time without changing the shape/position of waves

Example of work done with aaOcean:
https://vimeo.com/42087457

This repository contains the following

-aaOcean core class
-aaOcean Mental Ray shaders
-Shader Definitions for Mental Ray shaders
-aaOcean Arnold shaders
-aaOcean Softimage ICE deformer
-aaOcean Maya Deformer
-Houdini FFT node, which I used to port aaOcean as a deformer to Houdini
-several helper functions that I often use

LINUX:
All plugins/shaders come with makefiles
You will need single-precision, threaded fftw library installed. See 'externals/fftw' for a ./configure script that 
I use to generate a static library to link with. Download source from fftw.org
You will also need OpenEXR and Z-lib installed

MAYA & LINUX:
To compile the Maya plugin, you will need gcc-4.2.4. 
This version is ABI-compatible with gcc-4.1.x and supports OpenMP which is used
by aaOcean

WINDOWS:
Look at the linux makefiles to make your Visual Studio projects
You will need single-precision fftw library to link with. Download from fftw.org

