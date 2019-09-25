Algebraic Lens Distortion Model Estimation (Basic Version, only bmp images)

Authors : Luis Alvarez, Luis Gomez, Rafael Sendra



INFORMATION ABOUT PROGRAM COMPILATION AND EXECUTION 

Step 1 : 
The content of this zip archive is : (This code is distributed under the terms of the GPLv3 license.)

        a) lens_distortion.h, lens_distortion.cpp  : lens distortion model estimation 
           library.
        b) ami_pol.h, ami_pol.cpp : polynomial real root estimation library.
        c) image : basic library to read/write BMP images.
        d) lens_distortion_estimation.cpp : main program with lens distortion 
           estimation function calls.
        e) readme.txt : this file.
        f) calibration_pattern.jpg : sample input image to test the program.
        g) calibration_pattern_line_primitives.dat : sample line primitives to test 
           the program.
        h) point2d.h : class to manage 2D points. 

Step 2 : To compile, once the zip archive files are saved in a directory we 
execute the "make" command : 

>make -e CXX=g++

Step 3 : To execute the sample test included we do :

>./lens_distortion_estimation calibration_pattern.bmp output.bmp calibration_pattern_line_primitives.dat output_results.dat


Function parameter explanation : 

(1) calibration_pattern.bmp is the input image that has to be BMP format.

(2) output.bmp is the output lens distortion corrected image, that has to be BMP format.

(3) calibration_pattern_line_primitives.dat is the input collection of points 
    which belong to distorted lines presented in the image. 
    This file is an ASCII file organized in the following way : First, in the file, 
        we write the number of lines, next, for each line we write the number of points 
        and the point coordinates. 

(4) output_results.dat is the output file where the coefficients of the estimated 
    distortion models and energy values associated to the problem are stored. 

Optional input parameters:
(5,6) [x_center y_center]: if this parameter is not indicated, the center of the 
loaded image is the center of the distortion (it is internally calculated from the 
width and height image sizes). If the x_center and y_center are indicated, 
these values are considered the center of distortion during the calculation carried 
out by the program. These values, x_center and y_center may be integer or double and 
are supposed to be in pixels units.

(7) [optimize_center]: if this parameter is not indicated, the center of distortion 
is not optimized by the patch search (pixel level) and gradient method (subpixel level) 
during the calculation carried out by the program. If so, it is optimized. This value 
is an integer number, which a value of 1 (optimize center) or 0 (do not optimize center).

Note that there are some valid combinations for the input parameters:

Basic configuration (first option):
(1,2,3,4): these input parameters must be always indicated (the center of the loaded 
image is the center of the distortion and it is not optimized)


Second option:
(1,2,3,4, optimize_center): i.e.:(1,2,3,4, 1)
the center of the loaded image is the center of the  distortion and it is optimized 
during calculations

Third option:
(1,2,3,4, x_center y_center): i.e.:(1,2,3,4, 100 250)
the point (100,200) is the center of the  distortion and it is not optimized 
during calculations

Fourth option:
(1,2,3,4, x_center y_center 1): i.e.:(1,2,3,4, 100 250 1)
the point (100,200) is the center of the  distortion and it is optimized 
during calculations


--------------------------------------------------------------------------------
Understanding the output_results.dat file:

This output file contains the solution for the following cases considered: 

  1. algebraic method and its improvement through the gradient 
    (without taking into account the zoom: see paper published in IPOL),
  2. algebraic method and its improvement through the gradient 
    (taking into account the zoom: see paper published in IPOL),


From these solutions (both valid), the one selected in the IPOL online demo 
as "solution" (best solution) is the solution  2, that is,  the one obtained 
by  first applying the algebraic method from the trivial solution,(with zoom) 
and then improving it by means of the simple gradient method  (with zoom). 

Only even distortion coefficients are used (k2, k4) and both are optimized.


When a parameter (distortion coefficient or the center of distortion) is 
optimized (active during the calculations) it is clearly indicated in the output file. 
For example,

MODIFIED VARIABLES THROUGH OPTIMIZATION:        k[2] k[4], (x0,y0)

indicates that k2 and k4 coefficients and the center of distortion are optimized 
(active during the calculations).

At the beginning of the output file it is indicated if the center of distortion 
is going to be optimized:

                        The center of distortion is going to be optimized.

or if not,

                        The center of distortion is not going to be optimized.



The set of solutions -for the case of not optimizing the center of distortion 
is presented as follows (see aclaration below):

--------------------------------------------------------------------------------
                                The center of distortion is not going to be optimized.

/* FIRST WE COMPUTE THE TRIVIAL SOLUTION FOR THE POLYNOMIAL GRADE 
Na=4 (note that this case includes the Na=2) */
Na=4; x0=(1,0,0,0,0)
        
        The center of distortion is shown. 
        The initial values for energies related to the problem are shown: 
            Emin, Vmin, D (see paper published in IPOL).
        The value of the initial (trivial) distortion parameters are shown.
        The CPU time is shown.
        
--------------------------------------------------------------------------------
-------------------------------------------------------------------------------- 

SOLUTION 1: ALGEBRAIC METHOD (FROM TRIVIAL SOLUTION) + GRADIENT METHOD. 
NO ZOOM APPLIED (degree of LDM polynom: 4)
2 parameters, 1 iteration

MODIFIED VARIABLES THROUGH GRADIENT OPTIMIZATION:       k[2]  k[4]


        The optimized values for energies related to the problem are shown: 
        Emin, Vmin, D (see paper published in IPOL).
        The value of the final distortion parameters are shown.
        The center of distortion is shown.
        The number of gradient iterations and the number of Distance (D) 
                function evaluations are shown.
        The CPU time is shown.


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
SOLUTION 2: ALGEBRAIC METHOD (FROM TRIVIAL SOLUTION) + GRADIENT METHOD. 
ZOOM APPLIED (degree of LDM polynom: 4)
2 parameters, 1 iteration

MODIFIED VARIABLES THROUGH GRADIENT OPTIMIZATION:       k[2]  k[4]

        The optimized values for energies related to the problem are shown: 
        Emin, Vmin, D (see paper published in IPOL).
        The value of the final distortion parameters are shown.
        The center of distortion is shown. 
        The number of gradient iterations and the number of Distance (D) 
                function evaluations are shown.
        The CPU time is shown.
--------------------------------------------------------------------------------

The set of solutions for the case of optimizing the center of distortion
is presented as follows (see aclaration below):

--------------------------------------------------------------------------------
                                The center of distortion is going to be optimized.

/* FIRST WE COMPUTE THE TRIVIAL SOLUTION FOR THE POLYNOMIAL GRADE Na=4 
(note that this case includes the Na=2) */
Na=4; x0=(1,0,0,0,0)
        
        The center of distortion is shown. 
        The initial values for energies related to the problem are shown: 
        Emin, Vmin, D (see paper published in IPOL).
        The value of the initial (trivial) distortion parameters are shown.
    The CPU time is shown.
        
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

SOLUTION 1: ALGEBRAIC METHOD (FROM TRIVIAL SOLUTION) + GRADIENT METHOD. 
NO ZOOM APPLIED (degree of LDM polynom: 4)
2 parameters, 1 iteration

MODIFIED VARIABLES THROUGH GRADIENT OPTIMIZATION:       k[2]  k[4], (x0,y0)


        The optimized values for energies related to the problem are shown: 
        Emin, Vmin, D (see paper published in IPOL).
        The value of the final distortion parameters are shown.
        The new center of distortion is shown.
    The number of gradient iterations and the number of Distance (D) 
    function evaluations are shown.
    The CPU time is shown.


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
SOLUTION 2: ALGEBRAIC METHOD (FROM TRIVIAL SOLUTION) + GRADIENT METHOD. 
ZOOM APPLIED (degree of LDM polynom: 4)
2 parameters, 1 iteration

MODIFIED VARIABLES THROUGH GRADIENT OPTIMIZATION:       k[2]  k[4], (x0,y0)

        The optimized values for energies related to the problem are shown: 
        Emin, Vmin, D (see paper published in IPOL).
        The value of the final distortion parameters are shown.
        The new center of distortion is shown. 
    The number of gradient iterations and the number of Distance (D) 
    function evaluations are shown.
    The CPU time is shown.
--------------------------------------------------------------------------------

