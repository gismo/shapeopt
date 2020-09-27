# shapeopt

This repository contains code for performing shape optimization with IGA
Written by Asger Limkilde, in collaboration with Angelos Mantzaflaris

The code was written as part of Asger Limkildes PhD at the technical university of denmark (DTU).

Two main apporaches are implemented:

    - A boundary driven approach where a parametrization technique to generate inner control points from the boundary control points are use, implemented by deriving from 'gsParamMethod'. This class is the combined with a state equation, typically a PDE, derived from 'gsStateEquation', and constraints on theJacobian determinant of the geometry map in 'gsDetJacConstraints'. It is combined in the class 'gsShapeOptProblem'. 

    - A regularization based approach where all control points (both the inner control points and the boundary control points) are treated as design variables. This approach is implemented in the class 'gsShapeOptWithReg'

A few different specific problems are implemented:
    
    - The 2D shape optimization problem of designing electromagnetic reflectors by maximizing the electrical field at a specific point. 'gsOptAntenna'

    - The 3D shape optimization problem of designing reflectors for free surface water waves. 'gsOptPotWaves'

    - The problem of finding parametrizations, formulated as a shape optimization problem with the final shape know. 'gsOptParam'.


Note: Only the files starting with "gs" (and main.cpp) are used. The other files are old version that will be updated or deleted at some point.


=============================================
 Guide for compilation of G+Smo and shapeopt
=============================================
NOTE: it is important that everything is done in the appropiate folder, such that you dont mess up the gismo folder with build file.. If you ar unlucky to build in the gismo folder you can clean up by running a command (DONT REMEMBER NAME ATM)

Go to your git directory, e.g.

	$ mkdir git
	$ cd git

Clone G+SMO from GitHub

	$ git clone https://github.com/gismo/gismo.git

Switch to the ie-shapeopt branch and fetch code

	$ git checkout ie-shapeopt
	$ git fetch (note: perhaps this step is not necessary)


Make a build folder that is seperate and go into this folder

	$ mkdir build
	$ cd build

Run cmake from the build folder linking to the source code

	$ cmake ../gismo

Run ccmake to choose the necessary options

	$ ccmake .

Choose options:

	CMAKE_BUILD_TYPE : RelWithDebInfo (To allow debugging, only use Release mode for fast execution when program is correct)

	GISMO_WITH_IPOPT : ON		(Gets and compiles Ipopt to do optimization)
	GISMO_BUILD_EXAMPLES : ON   (Optional, builds examples)

Press [c] and then [g] to generate

Use the new settings by running cmake

	$ cmake .

Perhaps with gcc as compiler and c++14 run cmake

	$ cmake . -DCMAKE_CXX_STANDARD=14 -DCMAKE_CXX_COMPILER=gcc

Build with make 

	$ make -j 3 (To use 3 cores)

------ Note ------

If there is an issue with the automatic download of the IpOpt code, you can download the file 

  	https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.6.tgz

manually and put it in the folder

	extensions/gsIpopt/IpOpt-prefix/src/

then in the file 

	../gismo/extensions/gsIpopt/CMakeLists.txt

comment out the line

   	#URL http://www.coin-or.org/download/source/Ipopt/Ipopt-${IPOPT_VER}.tgz 

and replace with the

	URL ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/src/IpOpt-${IPOPT_VER}.tgz

to tell gismo to look for the tgz file here and not on coin-or.org/

------------------


=================================== 
Using shapeopt
=================================== 

Get the code from git

	$ git clone https://github.com/gismo/shapeopt.git

Now we need to run cmake.

------ Note ------

You need c++ > 11 to compile the code. If cmake version is >3.1 you can write 

	set(CMAKE_CXX_STANDARD, 14) 

in the CMakeList.txt file, if cmake version is <3.1 you can write
	
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14")

------------------

again we have to set 

	CMAKE_BUILD_TYPE = RelWithDebInfo
	CMAKE_CXX_STANDARD = 14


