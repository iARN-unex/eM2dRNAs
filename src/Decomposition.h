#ifndef DECOMPOSITION_H_
#define DECOMPOSITION_H_

#include "Global.h"
#include "Problem.h"

class Decomposition
{
  private:
	double index = 1.8;//1.4; // O(n^index) folding
	Problem* global; // global: a not decomposed problem
        vector<Problem*> problems;
	int totalParticipation;
	string targetRNA;
	int randomSeed;
	int popSize;
	double stopCrit;
	bool fastMode;
	int model;
	bool solvingFlag;
	bool decompositionFlag;

  public:
        Decomposition(string targetRNA_, int randomSeed_, int popSize_, double stopCrit_, bool fastMode_, int model_, bool solvingFlag_, bool decompositionFlag_);
	void updateParticipation();
	void decompose();
	Problem* recursive(string s);
	void topologicalSort(bool* visited, int numberOfSubproblems);
        void topologicalSortUtil(int v, bool visited[], stack<Problem*>& Stack);
	Problem* checkIfExists(string s);
        vector<Problem*> getProblems();
        int getProblemSize();
        Problem* getProblemAt(int i);
        Problem* getProblemByID(int id);
        void addProblem(Problem* p);
        void deleteProblem(int i, bool update);
        void deleteProblemByID(int id, bool update);
        void print();
        void printGraphviz();
	void solve();
	void printSolutions();
	~Decomposition();
};

#endif /* DECOMPOSITION_H_ */
