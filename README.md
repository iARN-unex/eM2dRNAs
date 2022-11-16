## eM2dRNAs: Enhanced Multiobjective Metaheuristic to Design RNA Sequences
### Overview
eM2dRNAs is a Multiobjective Evolutionary Algorithm (MOEA) to solve the RNA inverse folding problem. 

eM2dRNAs is an enhanced version of m2dRNAs, a previously developed MOEA approach to that problem. The main improvement is the recursive decomposition of the target structure into substructures, smaller and consequently easier to solve than the original one. These substructures will be solved independently of each other and subsequently combined to obtain a complete nucleotide sequence solution. The algorithm responsible for solving each substructure is m2dRNAs.

### Algorithm
Briefly, the process followed by eM2dRNAs is:

1. Decompose the target structure into substructures (RNAinv problems) in order to create a directed acyclic graph.
2. Deletion of certain RNAinv problems:
   - If the result of multiplying substructure size by its number of referenced times is lesser than 20. (e.g. ID4)
   - Every substructure containing only one substructure will be deleted, independently of its size. (e.g. ID3)
   - Whether a substructure is the only child of the target structure. (e.g. ID1)
3. The *participation* of each problem is calculated to proportionally distribute the given execution time (global stopping criterion) among the different RNAinv problems to be solved.
4. Topological sorting of the RNAinv problems to be solved. This is a lineal ordering of all the nodes of the directed acyclic graph previously created.
5. Solve each RNAinv problem independently using a modified version of the m2dRNAs algorithm, which is able to manage the subproblems contained by the input substructure. m2dRNAs algorithm receives as input parameters the substructure to be solved and the stopping criterion. The output will be a set of all the found RNA sequences, which will be saved. After the resolution of each problem, the participation of every remaining problem must be recalculated, since it is possible that the assigned time is not exhausted. This way, easier problems do not waste time that could be needed to solve more complex problems. If m2dRNAs fails to solve a subproblem, eM2dRNAs will eliminate that subproblem and restructure the graph.
6. The last problem to be solved will be the one containing the target structure so, after solving it, all the RNA sequences will be reported.


![eM2dRNAs-image](https://user-images.githubusercontent.com/118007536/201639695-5b13b959-b435-4cbc-b50a-ba12f1866006.png)

### Comparative Study
The performance of eM2dRNAs was evaluated with the widely used Eterna100 benchmark in its V1 and V2 version, as well as both Turner1999 and Turner2004 energy parameters sets. Input target structures can be found in [data/input](data/input) folder.

The parameter configuration of eM2dRNAs was: population size of 52 individuals (eM2dRNAs requires a population size multiple of 4 to perform
the binary tournament) and stopping criterion based on time (24 hours). All structures were attempted ten times by each combination of Eterna100-version/Turner-version. Resulting sequences can be found in [data/output](data/output) folder.

Results were compared against other RNA inverse folding methods.


### Instructions
To compile and execute eM2dRNAs in your local machine the following are needed: 

* [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/) 
* RNALib library

#### Installation guide

(For ubuntu 18.04, change whatever fits for your system).

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

Inside [src](src) folder, compile:

```
make clean
make
```

In order to test the installation, you can either type "./e_m2dRNAs -h" or just "./e_m2dRNAs" to get a help message.

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
./e_m2dRNAs 1 $(cat ../data/input/Eterna100-V1/eterna62.ss) 50 -60 turner1999 1 1
```
