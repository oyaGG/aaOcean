aaOcean v2.6
Author: Amaan Akram 
www.amaanakram.com
aaOcean is an implementation of Jerry Tessendorf's 2004 paper on Simulating Ocean Waves.
aaOcean is free software and can be redistributed and modified under the terms of the 
GNU General Public License (Version 3) as provided by the Free Software Foundation.
GNU General Public License http://www.gnu.org/licenses/gpl.html

Documentation:
https://bitbucket.org/amaanakram/aaocean/wiki/browse/

**FEATURES**
* Multi-threaded via OpenMP
* FFT spectrum sampling ensures that the resulting ocean shape is only enhanced by increasing ocean resolution, 
and not completely changed
* OpenEXR output for object-space vector displacement

Example of work done with aaOcean:
https://vimeo.com/42087457

This repository contains the following

* aaOcean core class
* aaOcean Mental Ray shaders
* Softimage Shader Definitions for Mental Ray shaders
* aaOcean Arnold shaders
* aaOcean Prman 19 RIS displacement shader
* aaOcean Softimage ICE deformer
* aaOcean Maya Deformer
* aaOcean Houdini SOP
* aaOcean standalong terminal/shell application
* everal helper functions that I often use

*****LINUX BUILD INSTRUCTIONS*******
Please see the wiki page at 
https://bitbucket.org/amaanakram/aaocean/wiki/Build%20instructions

MAYA & LINUX:
To compile the Maya plugin, you will need gcc-4.2.4 which is ABI-compatible with gcc-4.1.x and
supports OpenMP which is used by aaOcean

Acknowledgements for help and bug fixes: Frederic Servant, Fabrice Macagno, Phil Stopford, Andrew Helmer