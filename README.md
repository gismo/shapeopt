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


# Guide for compilation of G+Smo and shapeopt

This guide is only tested on Ubuntu (>18.XX).

Note: it is important that everything is done in the appropiate folder, such that you dont mess up the gismo folder with build file.. If you ar unlucky to build in the gismo folder you can clean up by running a command (see G+Smo directions)

Go to your git directory, e.g.

	$ mkdir git
	$ cd git

Clone G+SMO from GitHub

	$ git clone https://github.com/gismo/gismo.git

Switch to the ie-shapeopt branch and fetch code

   $ cd git
	$ git checkout ie-shapeopt
	$ git fetch (note: perhaps this step is not necessary)


Make a build folder that is seperate and go into this folder

   $ cd ..
	$ mkdir build
	$ cd build

Run cmake from the build folder linking to the source code

	$ cmake ../gismo

Run ccmake to choose the necessary options

	$ ccmake .

Choose options:

	CMAKE_BUILD_TYPE : RelWithDebInfo (To allow debugging, only use Release mode for fast execution when program is correct)
   ( Obs: I had issues running the code in pure release mode, so be aware if going from RelWithDebInfo to release )

	GISMO_WITH_IPOPT : ON		(Gets and compiles Ipopt to do optimization)
	GISMO_BUILD_EXAMPLES : OFF   (Optional, builds examples)

Press [c] and then [g] to generate

Use the new settings by running cmake. Perhaps with gcc as compiler and c++14 run cmake 

	$ cmake . -DCMAKE_CXX_STANDARD=14 -DCMAKE_CXX_COMPILER=gcc

Build with make 

	$ make -j 3 gismo (To use 3 cores)

-------------
Note 

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

# Using shapeopt

Get the code from git

	$ git clone https://github.com/gismo/shapeopt.git

Now we need to run cmake.

-------------
Note

You need c++ > 11 to compile the code. If cmake version is >3.1 you can write 

	set(CMAKE_CXX_STANDARD, 14) 

in the CMakeList.txt file, if cmake version is <3.1 you can write
	
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14")

------------------

again we have to set 

	CMAKE_BUILD_TYPE = RelWithDebInfo
	CMAKE_CXX_STANDARD = 14


# Overview of the classes and the code

OptParam: Doing parametrization by means of shapeopt


class gsParamMethod

   This class is the abstract interface of a parametrization method, ie. a map from the boundary controlpoints to the interior controlpoints.
   Most importantly this method is the main 'bookkeeper' in my code. 

   One of the main members is 'm_mappers' which holds a vector with a gsDofMapper for each dimension. For each dimension this Dofmapper is a 'Degrees of Freedom Mapper' that maps local indicies to global indicies of the controlpoints. This also takes care of the coupling of interfaces and possible elimination/ fixation of the boundary controlpoints. Finally by tagging controlpoints we can flag the ones that has a special role later, e.g. the design variables of a shape optimixation step.
   For performing OptParam the tagged points denotes the boundary controlpoints that should end up in the configuration of the 'goal' (the shape we want to parametrize)
   When constructing the class the default behaivour is to glue the interfaces, fix and tag the controlpoints on the boundary. 

   In the following methods when handling controlpoints we use the format 

      c = [ c_x c_y c_z ]^T
   
  The important methods are
 
     * getTagged()                 // Returns a vector of the tagged controlpoints (cps)
     * updateTagged( gsVector<> x) // Updates the geometry from a vector of tagged cps
     * getFree()                   // Returns a vector of free controlpoints
     * updateFree( gsVector<> x)   // Updates the geometry from a vector of free cps
     * getFixed()                  // Returns a vector of fixed controlpoints
     * updateFixed( gsVector<> x)  // Updates the geometry from a vector of free cps
     * getControlPoints()          // Get a vector of all cps
     * updateControlPoints( gsVector<> x)  // Updates the geometry from a vector of all cps
     * getFlat()                   // Get a vector of all cps using local indicies (one patch at a time). This format is independent of the mappings in m_mappers 
   

class gsShapeOptProblem

   The base class for shape optimization problems. It inherits from the gsOptProblem<> from G+Smo which interfaces IpOpt.

   The methods 

     * evalObj()
     * gradObj()          // gradient with respect to the tagged cps
     * setupMappers()     // should construct the mapper for the specific problem

   needs to be overloaded for the specific problem. Other important methods are:

     * gradAll()          // gradient with respect to all cps
     * evalCon()          // evaluation of constraints
     * mapper_grad        // should return a gsDofMapper from the indicies used in gradAll and the indicies used in the problem

# Parametrizations by means of shape optimization

class gsOptParam

   This class is an interpretation of the parametrization challenge as a shape optimization problem
 
   Key members are: 
     * gsMultiPatch<> m_mp              // Holds the current design (inherited from gsShapeOptProblem)
     * gsSpringMethod m_mp_goal         // Holds the goal design as a gsParamMethod. The inner cps do not matter, it is only used to hold the boundary cps, and for bookkeeping.
                                        // In practice it is just there to return the mappers of the cps of the goal domain (in this case the boundary cps)
     *  

