aaOcean v2.6
Author: Amaan Akram 
www.amaanakram.com
aaOcean is free software and can be redistributed and modified under the terms of the 
GNU General Public License (Version 3) as provided by the Free Software Foundation.
GNU General Public License http://www.gnu.org/licenses/gpl.html


Some key features
-Multi-threaded via OpenMP
-Has gone through a few levels of code optimization via profiling
-FFT spectrum sampling ensures that the resulting ocean shape is only 
enhanced by increasing ocean resolution, and not completely changed
-Experimental support for shaders to output OpenEXR image sequences 
for vector displacement
-Comes with deformers and shaders for various 3D apps and renderers

Example of work done with aaOcean:
https://vimeo.com/42087457

This repository contains the following

-aaOcean core class
-aaOcean Mental Ray shaders
-Softimage Shader Definitions for Mental Ray shaders
-aaOcean Arnold shaders
-aaOcean Softimage ICE deformer
-aaOcean Maya Deformer
-aaOcean Houdini SOP
-several helper functions that I often use

LINUX:
All plugins/shaders come with makefiles
You will need single-precision, threaded fftw library installed. See 
'externals/fftw' for a ./configure script that I use to generate 
a static library to link with. Download source from fftw.org

If you enable OpenEXR support, you will also need:
IlmBase
OpenEXR
zlib

MAYA & LINUX:
To compile the Maya plugin, you will need gcc-4.2.4. 
This version is ABI-compatible with gcc-4.1.x and supports OpenMP which is used
by aaOcean

WINDOWS:
Look at the linux makefiles to make your Visual Studio projects
You will need single-precision fftw library to link with. Download from fftw.org

