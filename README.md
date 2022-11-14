## eM2dRNAs: Enhanced Multiobjective Metaheuristic to Design RNA Sequences
### Overview
eM2dRNAs is a Multiobjective Evolutionary Algorithm (MOEA) to solve the RNA inverse folding problem.

![eM2dRNAs-image](https://user-images.githubusercontent.com/118007536/201639695-5b13b959-b435-4cbc-b50a-ba12f1866006.png)

### Installation
To compile and execute eM2dRNAs in your local machine the following are needed: 

* [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/) 
* RNALib library

#### Instructions

(For ubuntu 18.04, change whatever is needed in your system).

Install required packages:

```
sudo apt-get update
sudo apt-get install build-essential libgsl23 libgslcblas0
```
Download and install ViennaRNA package (not required to compile):

```
wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_18_04/viennarna_2.4.17-1_amd64.deb
sudo dpkg -i viennarna_2.4.17-1_amd64.deb
```

Download and install RNALib 2.4.17 libraries:

```
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.17.tar.gz
tar -zxvf ViennaRNA-2.4.17.tar.gz
cd ViennaRNA-2.4.17
./configure
make
sudo make install
```

Inside eM2dRNAs-source folder, compile:

```
make clean
make
```

In order to test the installation, you can either type "./e_m2dRNAs -h" or only "./e_m2dRNAs" to get a help message.

### Use guide

To execute eM2dRNAs just type from its folder:

```
./e_m2dRNAs <random-seed> <target-structure> <population-size> <stopping-criteria> <model> <solving-flag>
```

Being the parameters:

```
<random-seed>: number between 1-31 to generate a random-seed used.
<target-structure>: dot-bracket notation structure to solve.
<population-size>: number of individuals, the algorithm will output this number of RNA sequences at most.
<stopping-criteria>: positive or negative number. If positive, it indicates a number of evaluations. If negative, it indicates a number of seconds to evolve.
<model>: the following two models are available: TURNER1999 and TURNER2004.
<solving-flag>: If it is 1, it indicates that the algorithm stops as soon as it founds a valid sequence.
<decomposition-flag>: Set to 1 to decompose. Otherwise, it indicates that the algorithm will not decompose the target RNA structure (classical m2dRNAs).
```

All the parameters are mandatory.


Examples:

```
./e_m2dRNAs 1 "..((((....))))((((....))))((((...))))" 50 2700 TURNER2004 0 1
./e_m2dRNAs 1 $(cat ./data/input/Eterna100-V1/eterna62.ss) 50 -60 turner1999 1 1
```
