#ifndef PROBLEM_H_
#define PROBLEM_H_

#include "Global.h"

class Problem
{
  private:
	int ID;
        string targetRNA;
	string mask;
        vector< pair <int, Problem*> > children; // A problem may contain several problems
	vector< pair <int, Problem*> > parents; // A problem may be refered by several problems
        vector<string> sequences;
	double occupation;	// percentage of occupation within the Global Structure
	double participation;   // percentage of participations in the Evolution
	double stopCrit;
	bool attempted; // indicates if this problem has been attempted to solve
  public:
        Problem();
	int getID();
	void setID(int id);
        bool isAttempted();
	double getStopCrit();
	void setStopCrit(int stop);
        string getTargetRNA();
        string getTargetRNASubseq();
        string getCompressedTargetRNASubseq();
        string getMask();
	void setMask(string m);
        void setTargetRNA(string s, double occupation_);
	double getOccupation();
	double getParticipation();
	void setParticipation(double participation_);
	vector< pair <int, Problem*> > getParents();
        Problem* getParentAt(int i);
        pair<int, Problem*> getParentPairAt(int i);
	int getParentsSize();
	void addParent(int pos, Problem* p);
	void deleteParent(int id);
        vector< pair<int, Problem*> > getChildren();
        int getChildrenSize();
        Problem* getChild(int id);
        pair<int, Problem*> getChildrenPairAt(int i);
        void addChild(int pos, Problem* p);
	void deleteChild(int id);
        vector<string>* getSequences();
        int getSequencesSize();
        string getSequenceAt(int i);
        void addSequence(string s);
        void deleteSequence(int i);
        void print();
	void printGraphviz();
	void printSequences();
	void solve(int randomSeed, int popSize, double stopCrit, bool solvingFlag);
        double findStringSimilarity(std::string first, std::string second);
        int getEditDistance(std::string first, std::string second);
	void limitSequences(int n);
	void remove();
	void sortAndUnique();
	void sortChildrenParents();
	~Problem();
};

#endif /* PROBLEM_H_ */
