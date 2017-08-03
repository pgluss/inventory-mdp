# Markov Decision Process Simulation

This code simulates different models presented by the University of Washington's Foster School
of Business. The models use dynamic programming to simulate return and production of servers. 

The basic model is in C++ and the other models are written in D.

## Installing Requirements

1. Install D
    1. Download Windows installer from [here](http://downloads.dlang.org/releases/2.x/2.075.0/dmd-2.075.0.exe)
	2. Double click the downloaded executible to install
    3. Choose "no" and do not install 64-bit support (You might need 64-bit support for very large state spaces)
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

## To Do

- [ ] Update boundary conditions
- [X] Implement new lower boundary on N in the extended model
- [X] Implement CSV file outputs for data analysis in Excel and R
    - [X] Initial model
    - [X] Extended model
    - [X] Heuristic
- [X] Move basic model from C++ to D
- [X] Comment and upload R code for heatmaps
    - [X] Commented R code
    - [X] Make heatmap boundaries flexible
- [ ] Edit code to take in batches of inputs
- [X] Fully comment code for non-tech users
    - [X] Initial model
    - [X] Extended model
    - [X] Heuristic
- [ ] Implement heuristic extension comparison
- [ ] Add printed checkpoints to code
