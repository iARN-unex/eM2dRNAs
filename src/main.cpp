
#include "Decomposition.h"
#include "Rand.h"

bool found;
clock_t runtime_start;
int iPrint;
double stepTime; // Used for studying converge

int main(int argc, char *argv[])
{
/*
  for (int n=0; n < 10; n++)
  {
    int seed = ((atoi(argv[1])*31) + 111)%1235;
    Rand *aleatorio = new Rand(seed);
    cout << "------------------------------------------------" << endl;
    cout << "RND_UNI_INIT: " << aleatorio->get_rnd_uni_init() << endl;
    cout << "------------------------------------------------" << endl;

    for (int i=0; i<100; i++)
      cout << aleatorio->rndreal(0.0, 1.0) << endl;

    free(aleatorio);
    cout << endl;
  }
  exit(1);
*/
  runtime_start = clock();
  found = false;

  cout << "*****************************************************************************************************************************************************" << endl;
  cout << "*                                                                                                                                                   *" << endl;
  cout << "*                                     ENHANCED MULTIOBJECTIVE METAHEURISTIC TO DESIGN RNA SEQUENCES (eM2dRNAs)                                      *" << endl;
  cout << "*                                                           A. Rubio-Largo et al. (2022)                                                            *" << endl;
  cout << "*                                                                                                                                                   *" << endl;
  cout << "*****************************************************************************************************************************************************" << endl;

  if (argc<8)
  {
    cout << endl << "Usage " << argv[0] << " <random-seed> <target-structure> <population-size> <stopping-criteria> <model> <solving-flag>" << endl;
    cout << "\t<random-seed>: number between 1-31 to generate a random-seed used." << endl;
    cout << "\t<target-structure>: dot-bracket notation structure to solve." << endl;
    cout << "\t<population-size>: number of individuals, the algorithm will output this number of RNA sequences at most." << endl;
    cout << "\t<stopping-criteria>: positive or negative number. If positive, it indicates a number of evaluations. If negative, it indicates a number of seconds to evolve." << endl;
    cout << "\t<model>: the following two models are available: TURNER1999 and TURNER2004." << endl;
    cout << "\t<solving-flag>: If it is 1, it indicates that the algorithm stop as soon as it found a valid sequence." << endl;
    cout << "\t<decomposition-flag>: Set to 1 to decompose, Otherwise, it indicates that the algorithm will not decompose the target RNA structure (classical m2dRNAs)." << endl;
    cout << endl;
    cout << "\tExample: " << argv[0] << " 1 \"..((((....))))((((....))))((((...))))\" 50 2700 TURNER2004 0 1" << endl << endl;
    exit(1);
  }

  int randomSeed = atoi(argv[1]);
  string targetRNA(argv[2]);
  int popSize = atoi(argv[3]);
  double stopCrit = atof(argv[4]);
  bool fastMode = 0; // atoi(argv[5]);
  int model;

  for(char* c=argv[5]; *c=toupper(*c); ++c) ;
  if (strcmp(argv[5], "TURNER1999") == 0) 	model = MODEL_TURNER1999;
  else if (strcmp(argv[5], "TURNER2004") == 0) 	model = MODEL_TURNER2004;
  else
  {
    cout << endl << " **** WARNING **** Selected MODEL not recognized:" << argv[6] << ", try with: turner1999 or turner2004" << endl << endl;
    exit(1);
  }

  bool solvingFlag = false;
  if (atoi(argv[6]) > 0) solvingFlag = true;

  bool decompositionFlag = true;
  if (atoi(argv[7]) <= 0) decompositionFlag = false;

  cout << " PARAMETERS:" << endl;
  cout << "   * Random-seed:       " << randomSeed << endl;
  cout << "   * Target structure:  " << targetRNA  << endl;
  cout << "   * Population size:   " << popSize    << endl;
  if (stopCrit < 0)
  cout << "   * Stopping criteria: " << -stopCrit << " seconds " << endl;
  else
  cout << "   * Stopping criteria: " <<  stopCrit << " evaluations " << endl;
  cout << "   * Model:             " << argv[5] << endl;
  cout << "   * Solving flag:      " << solvingFlag << endl;
  cout << "   * Decomposition flag:      " << decompositionFlag << endl;
  cout << endl;
  cout << "*****************************************************************************************************************************************************" << endl;

  Decomposition d (targetRNA, randomSeed, popSize, stopCrit, fastMode, model, solvingFlag, decompositionFlag);
  d.decompose();
//  d.print();
//  d.printGraphviz(); 
//  cout << d.getProblemSize() << endl;
//  exit(1);
  d.solve();
//  d.print();

  d.printSolutions();
  clock_t t2 = clock() - runtime_start;
  printf("\n@TIME: %f secs\n", (((double) t2) / CLOCKS_PER_SEC));
  printf("\n");

  return EXIT_SUCCESS;
}
