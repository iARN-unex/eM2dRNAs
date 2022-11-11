#ifndef NSGA2_H_
#define NSGA2_H_

#include "Rand.h"
#include "Global.h"
#include "Problem.h"

# define RnaInvFold

# define INF 1.0e14
# define EPS 1.0e-14
# define E  2.71828182845905
# define PI 3.14159265358979
# define GNUPLOT_COMMAND "gnuplot -persist"

typedef struct
{
    int rank;
    double constr_violation;
    double *xreal;
    int **gene;
    double *xbin;
    double *obj;
    double EuclDistance;
    double *constr;
    double crowd_dist;
    int *rnaMatches;
    char* mfe;
    char* sequence;
}
individual;

typedef struct
{
    individual *ind;
}
population;

typedef struct lists
{
    int index;
    struct lists *parent;
    struct lists *child;
}
list_;


struct CompareCStrings
{
    bool operator()(const char* lhs, const char* rhs) const {
        return std::strcmp(lhs, rhs) < 0;
    }
};


class Nsga2
{
  private:
	int nreal;
	int nbin;
	int nobj;
	int ncon;
	int popsize;
	double pcross_real;
	double pcross_bin;
	double pmut_real;
	double pmut_real_subs;
	double pmut_bin;
	double eta_c;
	double eta_m;
	int ngen;
	int nbinmut;
	int nrealmut;
	int nbincross;
	int nrealcross;
	int *nbits;
	double *min_realvar;
	double *max_realvar;
	double *min_binvar;
	double *max_binvar;
	int bitlength;
	int choice;
	int obj1;
	int obj2;
	int obj3;
	int angle1;
	int angle2;

	int len; // targetRNA length
	char* targetRNA;
	int* initLoopTargetRNA;
	int* endLoopTargetRNA;
	int nLoops;
	int* pointTargetRNA;
	int nPoints;
	int targetRNALoops;
	int targetRNAPoints;
	int *ssTargetRNA;    // what substructure?
	int *ssInitTargetRNA; // what position in global structure?
	int *ssEndTargetRNA; // what position in global structure?
	int nSS;
	int SimilarityGlobal;
	/* Evaluation */
	char* seq;
	//char* mfe_structure;
	char* pfold_structure;
	int* similarityVector;
	int nEvals;
	/**************/
	double* bestObj;
	bool changeBestObj;
	double minSim;
	double probabilitiesLoops[6]; /* CG GC AU UA GU UG */
	double probabilitiesUnpaired[4]; /* C G A U */
	double maxGC;

	int fast_mode;
	bool solvingFlag;

	bool printSolution;
	int* failedPos;
	double stop_crit;
	int nMinSim;
	int nEvMFE;
	int nEvPFBP;
	int nEvOracle;
	int nConstr;

	map<char*, double*, CompareCStrings > oracle;
	map<char*, double* >::iterator itr;

	population *parent_pop;
	population *child_pop;
	population *mixed_pop;
	Problem* problem;
	vector<string> results;
	int validSolutionsFound;
	Rand *rndGenerator;
	char* minSimString;

  public:
        Nsga2(int sseed, int popsize_, double stop_crit_, bool solvingFlag_, Problem* problem_);

	vector<string> run();
	bool checkStopCrit(clock_t t, int stagnation);
	int getNEvals();
	void showResults();

	/****************************************/
	/* ALLOCATE */
	/****************************************/
	void allocate_memory_pop (population *pop, int size);
	void allocate_memory_ind (individual *ind);
	void deallocate_memory_pop (population *pop, int size);
	void deallocate_memory_ind (individual *ind);

	/****************************************/
	/* AUXILIARY */
	/****************************************/
	double maximum (double a, double b);
	double minimum (double a, double b);

	/****************************************/
	/* CROSSOVER */
	/****************************************/
	void crossover (individual *parent1, individual *parent2, individual *child1, individual *child2);
	void realcross (individual *parent1, individual *parent2, individual *child1, individual *child2);
	void bincross (individual *parent1, individual *parent2, individual *child1, individual *child2);

	/****************************************/
	/* CROWDIST */
	/****************************************/
	void assign_crowding_distance_list (population *pop, list_ *lst, int front_size);
	void assign_crowding_distance_indices (population *pop, int c1, int c2);
	void assign_crowding_distance (population *pop, int *dist, int **obj_array, int front_size);

	/****************************************/
	/* SORT */
	/****************************************/
	void quicksort_front_obj(population *pop, int objcount, int obj_array[], int obj_array_size);
	void q_sort_front_obj(population *pop, int objcount, int obj_array[], int left, int right);
	void quicksort_dist(population *pop, int *dist, int front_size);
	void q_sort_dist(population *pop, int *dist, int left, int right);

	/****************************************/
	/* DECODE */
	/****************************************/
	void decode_pop (population *pop);
	void decode_ind (individual *ind);

	/****************************************/
	/* DISPLAY */
	/****************************************/
	void onthefly_display (population *pop, FILE *gp, int ii);

	/****************************************/
	/* DOMINANCE */
	/****************************************/
	int check_dominance (individual *a, individual *b);

	/****************************************/
	/* EVAL */
	/****************************************/
	void evaluate_pop (population *pop);
	void evaluate_ind (individual *ind);

	/****************************************/
	/* FILLNDS */
	/****************************************/
	void fill_nondominated_sort (population *mixed_pop, population *new_pop);
	void crowding_fill (population *mixed_pop, population *new_pop, int count, int front_size, list_ *elite);

	/****************************************/
	/* FILLNDS */
	/****************************************/
	void insert (list_ *node, int x);
	list_* del (list_ *node);


	/****************************************/
	/* INITIALIZE */
	/****************************************/
	void initialize_pop (population *pop, vector<string>* ssPop);
	void initialize_ind (individual *ind, int i);
	void initialize_ind (individual *ind);
	void initialize_ind0 (individual *ind);
	void initialize_ind1 (individual *ind);
	void initialize_ind2 (individual *ind);
	void initialize_ind3 (individual *ind);
	void read_ind(individual *ind, string seq_);

	/****************************************/
	/* MERGE */
	/****************************************/
	void merge(population *pop1, population *pop2, population *pop3);
	void copy_ind (individual *ind1, individual *ind2);

	/****************************************/
	/* MUTATION */
	/****************************************/
	void mutation_pop (population *pop, int stagnation);
	void mutation_ind (individual *ind);
	void bin_mutate_ind (individual *ind);
	void real_mutate_ind (individual *ind);

	/****************************************/
	/* PROBLEM_DEF */
	/****************************************/
	double findStringSimilarity(std::string first, std::string second);
        int getEditDistance(std::string first, std::string second);
	void test_problem (double *xreal, int *rnaMatches, char* mfe, char* sequence, double *xbin, int **gene, double *obj, double *constr);

	/****************************************/
	/* RANK */
	/****************************************/
	void assign_rank_and_crowding_distance (population *new_pop);

	/****************************************/
	/* REPORT */
	/****************************************/
	void getResults(population *pop);
	void getOracle();
	void report_pop_xreal(population* pop);
	void report_ind (individual* ind);
	void report_pop_ (population *pop, int flag);
	void report_feasible_ (population *pop);
	void report_pop2 (population *pop);
	void report_feasible (population *pop, FILE *fpt);

	/****************************************/
	/* TOURSELECT */
	/****************************************/
	void selection (population *old_pop, population *new_pop);
	individual* tournament (individual *ind1, individual *ind2);

	~Nsga2();
};

#endif /* NSGA2_H_ */
