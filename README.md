# CryptCellSeparation

## Code for "Interkinetic nuclear migration and basal tethering facilitates post-mitotic daughter separation in intestinal organoids"

This project contains the code necesary for running the simulations presented in Caroll et al "Interkinetic nuclear migration and basal tethering facilitates post-mitotic daughter separation in intestinal organoids" in the Journal of Cell Science [https://doi.org/10.1242%2Fjcs.211656](https://doi.org/10.1242%2Fjcs.211656) 

This model is an extension on a model originally presented in Dunn et al. (2016) doi: 10.1091/mbc.E15-12-0854 and makes use of the cell cycle and geometry model given here https://chaste.cs.ox.ac.uk/trac/wiki/PaperTutorials/CryptProliferationDistribution. 

Before looking at this, you may wish to look at some of the basic user tutorials for Chaste https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials.

## Getting the code and installing dependencies 

Before running these examples you will need to install Chaste's dependencies (https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides) and the source code for version 3.4 (http://www.cs.ox.ac.uk/chaste/download.html).
The easiest way to do this is using an Ubuntu machine (or an Ubuntu virtual machine) as discussed on https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/UbuntuPackage. 
Note that Chaste is only fully supported on Linux/Unix systems, so users of Windows or Mac OS X may need to follow the virtual machine route.

You will also need the source for the CryptCellSeparation project.  This can be done by checking out the version from this github repository in to the projects folder of the Chaste directory.

Now the project should be installed, and everything should compile and run correctly. 
You can now run the tests or simulations, or create your own test suites.

## Documentation
There are two folders - `src` and `test`.
 1. The `src` folder contains the classes necesary to run the simulation. These define the additional forces and boundary conditions not in the core chaste code.
 2. The `test` folder contains:
  * TestCryptSingleRun - this file can be run to generate the results in Figure .
  * run_script.sh - Script to run multiple simultions

## Running tests
You can then run tests and simulations with,
```
cd <Chaste3.4 path>
scons b=GccOpt cl=0 co=1 ts=projects/CryptCellSeparation/test/TestCryptSingleRun.hpp
```

Note that this will only compile the test. The following commands will run the parammeter sweep detailed in the paper:
```
cd projects/CryptCellSeparation/test/
sh run_script.sh
```

**NB**: the paper was developed with release version 3.4. It will not work with with release version 3.3 or under.

For further information on using Chaste, see the extensive guide material (https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides).
You may also wish to look at some of the basic user tutorials (https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials).
