# Markov Decision Process Simulation

This code simulates different models presented by the University of Washington's Foster School
of Business. The models use dynamic programming to simulate return and production of servers. 

The basic model is in C++ and the other models are written in D.

## Installing Requirements

1. Install D
    1. Download Windows installer from [here](http://downloads.dlang.org/releases/2.x/2.075.0/dmd-2.075.0.exe)
	2. Double click the downloaded executible to install
2. Download source
	1. Download it [here](https://github.com/pgluss/inventory-mdp/archive/master.zip)
	2. Unzip the downloaded file

## Running the Code

1. Open a Windows terminal (or Powershell)
2. Navigate to the unzipped source code directory
3. Run code using one of the following commands:
	- dub --config=ha --build=release-nobounds
	- dub --config=basic --build=release-nobounds
	- dub --config=extension --build=release-nobounds
	- dub --config=heuristic --build=release-nobounds

## Modifying Parameters

The parameters of the model are declared at the beginning of each file. To modify them, open up 
the file in a text editor and directly modify the values. Rerunning the code as described above 
will automatically compile and use the new parameters.

## Model Outputs

To do!
