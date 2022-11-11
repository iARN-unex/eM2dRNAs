#include "Nsga2.h"
#include "Global.h"

using namespace std;

Nsga2::Nsga2(int sseed, int popsize_, double stop_crit_, bool solvingFlag_, Problem* problem_)
{
  problem = problem_;
  rndGenerator = NULL;
//  mfe_structure = NULL;
  pfold_structure = NULL;
  targetRNA = NULL;
  initLoopTargetRNA = NULL;
  endLoopTargetRNA = NULL;
  pointTargetRNA = NULL;
  ssTargetRNA = NULL;
  ssInitTargetRNA = NULL;
  ssEndTargetRNA = NULL;
  bestObj = NULL;
  min_realvar = NULL;
  max_realvar = NULL;

  targetRNALoops = 0;
  targetRNAPoints = 0;

  nbits = NULL;
  min_binvar = NULL;
  max_binvar = NULL;
  seq = NULL;
  similarityVector = NULL;
  failedPos = NULL;

  int i;
  int seed = ((sseed*31) + 111)%1235;
  rndGenerator = new Rand(seed);

  fast_mode = 0;//fast_mode_;
  solvingFlag = solvingFlag_;
  validSolutionsFound = 0;

  nMinSim = 0;
  results.clear();
  nEvMFE = 0;
  nEvPFBP = 0;
  nEvOracle = 0;
  nConstr = 0;

  /********************************************/
  nEvals=0;
  SimilarityGlobal = 0;

  targetRNA = (char *)malloc(sizeof(char) * (problem->getTargetRNA().size() + 1));
  minSimString = (char *)malloc(sizeof(char) * (problem->getTargetRNA().size() + 1));
  strcpy(targetRNA, problem->getTargetRNA().c_str());
  strcpy(minSimString, "");

  len = strlen(targetRNA);

  probabilitiesLoops[0] = 0.17; /*0.33; /*atof(argv[5]);  /* CG */
  probabilitiesLoops[1] = 0.17; /*0.33; /*atof(argv[6]);  /* GC */
  probabilitiesLoops[2] = 0.17; /*0.14; /*atof(argv[7]);  /* AU */
  probabilitiesLoops[3] = 0.17; /*0.14; /*atof(argv[8]);  /* UA */
  probabilitiesLoops[4] = 0.16; /*0.03; /*atof(argv[9]);  /* GU */
  probabilitiesLoops[5] = 0.16; /*0.03; /*atof(argv[10]); /* UG */

  probabilitiesUnpaired[0] = 0.25; /*0.19; /* C */
  probabilitiesUnpaired[1] = 0.25; /*0.17; /* G */
  probabilitiesUnpaired[2] = 0.25; /*0.38; /* A */
  probabilitiesUnpaired[3] = 0.25; /*0.26; /* U */

  initLoopTargetRNA = (int *)malloc(len*sizeof(int));
  endLoopTargetRNA = (int *)malloc(len*sizeof(int));
  pointTargetRNA = (int *)malloc(len*sizeof(int));
  ssTargetRNA = (int*) malloc(len * sizeof(int));
  ssInitTargetRNA = (int*) malloc(len * sizeof(int));
  ssEndTargetRNA = (int*) malloc(len * sizeof(int));

  int c, j;
  nLoops = 0, nPoints=0, nSS=0;
  nreal = 0;

  int dollar=0;
  // Add substructures (if any)
  if (problem->getChildrenSize() > 0)
  {
    pair<int, Problem*> p;
    for (int i=0; i < problem->getChildrenSize(); i++)
    {
      p = problem->getChildrenPairAt(i);
      if (p.second->getSequencesSize() > 0)
      {
        ssTargetRNA[nSS] = p.second->getID();
        ssInitTargetRNA[nSS] = p.first;
        ssEndTargetRNA[nSS] = p.first + p.second->getTargetRNA().size() - 1;
        nSS++;
	//for (int i=p.first; i<(p.first + p.second->getTargetRNA().size()); i++)
	//  targetRNA[i] = '*';
      }
    }
  }
  /* Add unpairs and pairs (or loops: CG, GC, AU, UA, GU, or UG) */
  for(i=0; i<len; i++){
    if (targetRNA[i] == '.'){
      pointTargetRNA[nPoints] = i;
      nPoints++;
    }
    if (targetRNA[i] == '('){
      initLoopTargetRNA[nLoops] = i;
      c = 0;
      for(j=i; j<len; j++){
        if (targetRNA[j] == '(') c++;
	if (targetRNA[j] == ')') c--;
	if( (c <= 0) && (targetRNA[j] == ')') ){
	  endLoopTargetRNA[nLoops] = j;
	  nLoops++;
	  break;
	}
      }
    }
  }
  problem->setMask(targetRNA);
//  cout << "   MASK:       " << targetRNA << endl;
  strcpy(targetRNA, problem->getTargetRNA().c_str());
//  cout << "   TARGET RNA: " << targetRNA << endl;
  len = strlen(targetRNA);


  if (nSS == 0)
  {
    targetRNAPoints = nPoints;
    targetRNALoops = nLoops;
  }
  else
  {
    // Compute total of nLoops and nPoints (including substructures)
    for(i=0; i<len; i++){
      if (targetRNA[i] == '.'){
        targetRNAPoints++;
      }
      if (targetRNA[i] == '('){
        c = 0;
        for(j=i; j<len; j++){
          if (targetRNA[j] == '(') c++;
          if (targetRNA[j] == ')') c--;
          if( (c <= 0) && (targetRNA[j] == ')') ){
            targetRNALoops++;
            break;
          }
        }
      }
    }
  }

  /* For evaluation */
  //seq = (char*) vrna_alloc(sizeof(char) * (len + 1));
  seq = (char*) malloc(sizeof(char) * (len + 1));
//  mfe_structure = (char*) vrna_alloc(sizeof(char) * (len + 1));
//  mfe_structure = (char*) malloc(sizeof(char) * (len + 1));
//  pfold_structure = (char*) vrna_alloc(sizeof(char) * (len + 1));
  pfold_structure = (char*) malloc(sizeof(char) * (len + 1));
  strcpy(seq,"");
//  strcpy(mfe_structure,"");
  strcpy(pfold_structure,"");
  /********************************************/

  nreal = nLoops + nPoints + nSS;

  printf("   PAIRED POSITIONS:   ");
  for (i=0; i<nLoops; i++){
    if (i%5 == 0 && i>0) printf("\n                       ");
    printf("(%d, %d) ", initLoopTargetRNA[i]+1, endLoopTargetRNA[i]+1);
  }
  printf("\n");
  printf("   UNPAIRED POSITIONS: ");
  for (i=0; i<nPoints; i++){
    if (i%10 == 0 && i>0) printf("\n                       ");
    printf("(%d) ", pointTargetRNA[i]+1);
  }
  printf("\n");

  printf("   SUBSTRUCTURES:      ");
  for (i=0; i<nSS; i++){
    if (i%2 == 0 && i>0) printf("\n                       ");
    printf("(%d, %d, $%d) ", ssInitTargetRNA[i] + 1, ssEndTargetRNA[i] + 1, ssTargetRNA[i]);
  }
  printf("\n");

  popsize=popsize_; //atoi(argv[3]);
  while (popsize % 4 != 0) popsize++;
  stop_crit = stop_crit_; //atoi(argv[4]);
  nobj=3;
  if(fast_mode == 1) nobj=2;
  ncon=1;
  min_realvar = (double *)malloc(nreal*sizeof(double));
  max_realvar = (double *)malloc(nreal*sizeof(double));
  bestObj = (double *)malloc(nobj*sizeof(double));
  changeBestObj = false;
  minSim = -99999.0;
  for(i=0; i<nobj; i++)
    bestObj[i]=999999;
  for (i=0; i<nreal; i++){
    if(i < nLoops+nPoints){
      min_realvar[i] = 0;
      max_realvar[i] = 1; /* CG (0 to 0.16), GC (0.16 to 0.33), AU (0.33 to 0.49), UA (0.49 to 0.66), GU (0.66 to 0.83), or UG (0.83 to 1.00) */
    }
    else{
      min_realvar[i] = -1.0;
      max_realvar[i] = 1; /* C (0.00 to 0.25), G (0.25 to 0.50), U (0.50 to 0.75), or A (0.75 to 1.00)*/		
    }
  }

  pcross_real = 0.6;
  pmut_real = 10.0 / ((nLoops + nPoints)*1.0); if (pmut_real > 0.15) pmut_real = 0.15; //if (pmut_real < 0.01) pmut_real = 0.01;
  pmut_real_subs = 0.9; // Fixed mutation for substructures
  eta_c = 10;/*atoi(argv[3]);*/
  eta_m = 20;/*atoi(argv[4]);*/
  nbin = 0;
  choice = 0; /* no gnuplot */
/*
  printf("   Probabilities Loops={GC/CG=%.3f, AU/UA=%.3f, GU/UG=%.3f} \n", probabilitiesLoops[0], probabilitiesLoops[2], probabilitiesLoops[4], probabilitiesLoops[6]);
  printf("   ProbUnpairs={C=%.3f, G=%.3f, A=%.3f, U=%.3f} \n", probabilitiesUnpaired[0], probabilitiesUnpaired[1], probabilitiesUnpaired[2], probabilitiesUnpaired[3]);
  printf("   NSGA2: \n");
  printf("    - Population Size: %d \n", popsize);
  printf("    - NEvals (stopping criterion): %f \n", stop_crit);
  printf("    - Simulated Binary Crossover with %.2f of probabilty (eta_c=%.0f) \n", pcross_real, eta_c);
  printf("    - Polynomial Mutation with %.5f of probabilty (eta_m=%.0f) \n\n", pmut_real, eta_m);
*/
  nbinmut = 0;
  nrealmut = 0;
  nbincross = 0;
  nrealcross = 0;
  parent_pop = (population *)malloc(sizeof(population));
  child_pop = (population *)malloc(sizeof(population));
  mixed_pop = (population *)malloc(sizeof(population));
  allocate_memory_pop (parent_pop, popsize);
  allocate_memory_pop (child_pop, popsize);
  allocate_memory_pop (mixed_pop, 2*popsize);
}

vector<string> Nsga2::run()
{
  clock_t t, t2;
  t = clock();
  bool stop = false;
  int stagnation = 0;
  int i;
  int repeated = 0;

  initialize_pop (parent_pop, problem->getSequences()); //NULL);//ssPop);
  decode_pop(parent_pop);
  evaluate_pop (parent_pop);
  do
  {
    initialize_pop (child_pop, NULL);
    decode_pop(child_pop);
    evaluate_pop (child_pop);
    merge (parent_pop, child_pop, mixed_pop);
    fill_nondominated_sort (mixed_pop, parent_pop);
    stop = checkStopCrit(t, stagnation);
    repeated++;
  } while (!stop && minSim != 0 && repeated < 5);
/*
  if (!stop)
  {
    initialize_pop (parent_pop, problem->getSequences()); //NULL);//ssPop);
    decode_pop(parent_pop);
    evaluate_pop (parent_pop);
    stop = checkStopCrit(t, stagnation);
  }

  if (!stop)
  {
    initialize_pop (child_pop, NULL);
    decode_pop(child_pop);
    evaluate_pop (child_pop);
    stop = checkStopCrit(t, stagnation);
  }
*/

  int reset=0;
  int first = 1;
  i=1;

  bool report = false;

  while (!stop)
  {
    selection (parent_pop, child_pop);
    mutation_pop (child_pop, stagnation);
    decode_pop(child_pop);
    evaluate_pop(child_pop);
    merge (parent_pop, child_pop, mixed_pop);
    fill_nondominated_sort (mixed_pop, parent_pop);
    stagnation++;

    if (changeBestObj)
    {
      stagnation = 0;
      changeBestObj = false;
    }

    stop = checkStopCrit(t, stagnation);

    if (i%5 == 0 || stop) // Show best results every 5 iterations
    {
      if(minSim == 0){
        if(fast_mode == 1)
          printf("\tIT#%d (%d evals).\tBest [MFE(f1), Euclid. Dist.(f2)] = [%.2f, %.5f]  --> stagnation: %d\n ", i, nEvals, bestObj[0], bestObj[1], stagnation);
        else
          printf("\tIT#%d (%d evals).\tBest [Part. Funct.(f1), Ensemb. Div.(f2), Euclid. Dist.(f3)] = [%.2f, %.2f, %.5f]  --> stagnation: %d \n ", i, nEvals, bestObj[0], bestObj[1], bestObj[2], stagnation);
      }
      else{
        printf("\tIT#%d (%d evals).\t[No successfully RNA sequence found] Min. Sim. = %.5f  --> stagnation: %d \n ", i, nEvals, minSim, stagnation);
        printf("\t %s \n\t %s \n", targetRNA, minSimString);
      }
    }
/*
    if ( stop_crit < 0 || (problem->getID()==0 && stop) )
    {
      clock_t t2 = clock() - runtime_start;
      if ( (((double) t2) / CLOCKS_PER_SEC) >= (iPrint * stepTime) || (problem->getID()==0 && stop) )
      {
        cout << endl;
        cout << "@BEGIN i: " << iPrint << "\tEVALS: " << nEvals << "\tTIME: " << (((double) t2) / CLOCKS_PER_SEC) << endl;
        if (problem->getID() == 0)
          report_pop_(parent_pop, 1);
        cout << "@iTIME: " << (((double) t2) / CLOCKS_PER_SEC) << endl;
        cout << "@END i: " << iPrint << "\tEVALS: " << nEvals << endl;
        cout << endl;
        iPrint++;
      }
   }
*/
/*
    if (problem->getID()==0 && stop_crit < 0 || problem->getID()==0 && stop) // Convergence study 
    {
      clock_t t2 = clock() - runtime_start;
      clock_t t3 = clock() - t;
      if ( (((double) t3) / CLOCKS_PER_SEC) >= iPrint*(-stop_crit/50.0) || stop ) // Print population 50 times (convergence study)
      {
        cout << endl;
        cout << "@BEGIN i: " << iPrint << "\tEVALS: " << nEvals << "\tTIME: " << (((double) t2) / CLOCKS_PER_SEC) << endl;
        report_pop_(parent_pop,1);
        cout << "@iTIME: " << (((double) t2) / CLOCKS_PER_SEC) << endl;
        cout << "@END i: " << iPrint << "\tEVALS: " << nEvals << endl;
        cout << endl;
        iPrint++;
     }
    }
*/
    i++;
  }
/*
    if (isSubstr)// || gSS != NULL)
    {
      // If dealing with substructure or partial global
      if ( (minSim == 0 && stagnation >= 5) || (minSim != 0 && stagnation >= 10) )
        stop = true;
      }

      if (stop_crit < 0)
      {
        clock_t t2 = clock() - runtime_start;
        if ((((double) t2) / CLOCKS_PER_SEC) >= (-stop_crit))
        {
          stop = true;
          report = true;
        }
      }
      else
      {
        if (!isSubstr && nEvals >= stop_crit)
        {
          stop = true;
          report = true;
        }
      }
      i++;
    }
*/
/*
    t2 = clock() - t;
    printf(
         "\n-TIME: %f secs\n-EVALUATIONS %d \n-nEvMFE: %d \n-nEvPFBP: %d \n-nEvOracle: %d \n-nConstr: %d \n", 
         (((double) t2) / CLOCKS_PER_SEC), nEvals, nEvMFE, nEvPFBP, nEvOracle, nConstr);
    //cout << "@ORACLE: " << getSizeOracle() << " sequences found" << endl;
    printf("\n");
*/
/*
    if ( problem->getID() == 0) //!isSubstr && report && printSolution)
    {
      t2 = clock() - runtime_start;
      printf("\n@TIME: %f secs\n", (((double) t2) / CLOCKS_PER_SEC));
      printf("\n");
//      report_pop_(parent_pop, 1);
    }
*/

//    getResults(parent_pop);
    getOracle();
    return results;
}

bool Nsga2::checkStopCrit(clock_t t, int stagnation)
{
    if (stop_crit > 0) // Stopping criterion: #Evaluations
    {
      int criterion = nEvals;
      if ( (criterion >= stop_crit) \
         || (minSim >= 0 && stagnation >= 10 && problem->getID() != 0) \
         || (solvingFlag && stagnation >= 5 && validSolutionsFound >= 50) )
        return true;
    }
    else // Stopping criterion: Runtime
    {
      clock_t criterion = clock() - t;
      //cout << "Nsga stop: " << (((double) criterion * 1.0) / CLOCKS_PER_SEC) << endl;
      if ( (((double) criterion * 1.0) / CLOCKS_PER_SEC) >= -stop_crit \
         || (minSim >= 0 && stagnation >= 10 && problem->getID() != 0) \
         || (solvingFlag && stagnation >= 5 && validSolutionsFound >= 50) )
        return true;
    }
  return false;
//      || (problem->getChildrenSize() >  0 && stagnation >= 5 && nSS > 0 && criterion <= 0.8*stop_crit) )
//      || (problem->getChildrenSize() >  0 && stagnation >= 5 && nSS > 0 && (((double) criterion * 1.0) / CLOCKS_PER_SEC) <= stop_crit*-0.8) \
//      || (minSim != 0 && stagnation >= 10 && problem->getID() != 0) \
//      || (minSim != 0 && stagnation >= 10 && problem->getID() != 0) \
//      || (problem->getChildrenSize() == 0 && stagnation >= 5 && criterion <= 0.8*stop_crit)
//      || (problem->getChildrenSize() == 0 && stagnation >= 5 && (((double) criterion * 1.0) / CLOCKS_PER_SEC) <= stop_crit*-0.8 ) 
}

void Nsga2::showResults()
{
  report_pop_(parent_pop, 1);
}

int Nsga2::getNEvals()
{
  return nEvals;
}

Nsga2::~Nsga2()
{
//    cout << "Eliminando NSGA-II: " << problem->getID() << endl; 
    deallocate_memory_pop(parent_pop, popsize);
    deallocate_memory_pop(child_pop, popsize);
    deallocate_memory_pop(mixed_pop, 2*popsize);
    if (rndGenerator != NULL)
      { /*cout << "FREE rndGenerator" << endl;*/ free(rndGenerator);}
    if (targetRNA != NULL)
      { /*cout << "FREE targetRNA" << endl;*/ free(targetRNA);}
    if (minSimString != NULL)
      { /*cout << "FREE minSimString" << endl;*/ free(minSimString);}
    if (initLoopTargetRNA != NULL)
      { /*cout << "FREE initLoopTargetRNA" << endl;*/ free(initLoopTargetRNA);}
    if (endLoopTargetRNA != NULL)
      { /*cout << "FREE endLoopTargetRNA" << endl;*/ free(endLoopTargetRNA);}
    if (pointTargetRNA != NULL)
      { /*cout << "FREE pointTargetRNA" << endl;*/ free(pointTargetRNA);}
    if (ssTargetRNA != NULL)
      { /*cout << "FREE ssTargetRNA" << endl;*/ free(ssTargetRNA);}
    if (ssInitTargetRNA != NULL)
      { /*cout << "FREE ssInitTargetRNA" << endl;*/ free(ssInitTargetRNA);}
    if (ssEndTargetRNA != NULL)
      { /*cout << "FREE ssEndTargetRNA" << endl;*/ free(ssEndTargetRNA);}
    if (bestObj != NULL)
      { /*cout << "FREE bestObj" << endl;*/ free(bestObj);}
    if (min_realvar != NULL)
      { /*cout << "FREE min_realvar" << endl;*/ free(min_realvar);}
    if (max_realvar != NULL)
      { /*cout << "FREE max_realvar" << endl;*/ free(max_realvar);}

    if (nbits != NULL)
      { /*cout << "FREE nbits" << endl;*/ free(nbits);}
    if (min_binvar != NULL)
      { /*cout << "FREE min_binvar" << endl;*/ free(min_binvar);}
    if (max_binvar != NULL)
      { /*cout << "FREE similarityVector" << endl;*/ free(max_binvar);}
    if (seq != NULL)
      { /*cout << "FREE seq" << endl;*/ free(seq);}
    //if (mfe_structure != NULL)
    //  { cout << "FREE mfe_structure" << endl; free(mfe_structure);}
    if (pfold_structure != NULL)
      { /*cout << "FREE pfold_structure" << endl;*/ free(pfold_structure);}
    if (similarityVector != NULL)
      { /*cout << "FREE similarityVector" << endl;*/ free(similarityVector);}
    if (failedPos != NULL)
      { /*cout << "FREE failedPos" << endl;*/ free(failedPos);}

//    cout << oracle.size() << endl;
    for (itr = oracle.begin(); itr != oracle.end(); itr++)
    {
      char* s = itr->first;
      double* o = itr->second;
      free(s);
      free(o);
    }
    oracle.clear();
//    cout << oracle.size() << endl;
//    cout << "END DELETE" << endl;
}

/*********************************************************************************************************************************************************************************/
/* ALLOCATE */
/*********************************************************************************************************************************************************************************/
/* Function to allocate memory to a population */
void Nsga2::allocate_memory_pop (population *pop, int size)
{
    int i;
    pop->ind = (individual *)malloc(size*sizeof(individual));
    for (i=0; i<size; i++)
    {
        allocate_memory_ind (&(pop->ind[i]));
    }
    return;
}

/* Function to allocate memory to an individual */
void Nsga2::allocate_memory_ind (individual *ind)
{
    int j;
    if (nreal != 0)
    {
        ind->xreal = (double *)malloc(nreal*sizeof(double));
        ind->rnaMatches = (int *)malloc(len*sizeof(int));
        ind->mfe = (char *)malloc((len+1)*sizeof(char));
        ind->sequence = (char *)malloc((len+1)*sizeof(char));
    }
    if (nbin != 0)
    {
        ind->xbin = (double *)malloc(nbin*sizeof(double));
        ind->gene = (int **)malloc(nbin*sizeof(int *));
        for (j=0; j<nbin; j++)
        {
            ind->gene[j] = (int *)malloc(nbits[j]*sizeof(int));
        }
    }
    ind->obj = (double *)malloc(nobj*sizeof(double));
    if (ncon != 0)
    {
        ind->constr = (double *)malloc(ncon*sizeof(double));
    }
    return;
}

/* Function to deallocate memory to a population */
void Nsga2::deallocate_memory_pop (population *pop, int size)
{
    int i;
    for (i=0; i<size; i++)
    {
        deallocate_memory_ind (&(pop->ind[i]));
    }
    free (pop->ind);
    return;
}

/* Function to deallocate memory to an individual */
void Nsga2::deallocate_memory_ind (individual *ind)
{
    int j;
    if (nreal != 0)
    {
        free(ind->xreal);
//        free(ind->rnaMatches);
        free(ind->mfe);
        free(ind->sequence);
    }
    if (nbin != 0)
    {
        for (j=0; j<nbin; j++)
        {
            free(ind->gene[j]);
        }
        free(ind->xbin);
        free(ind->gene);
    }
    free(ind->obj);
    if (ncon != 0)
    {
        free(ind->constr);
    }
    return;
}

/*********************************************************************************************************************************************************************************/
/* AUXILIARY */
/*********************************************************************************************************************************************************************************/
/* Function to return the maximum of two variables */
double Nsga2::maximum (double a, double b)
{
    if (a>b)
    {
        return(a);
    }
    return (b);
}

/* Function to return the minimum of two variables */
double Nsga2::minimum (double a, double b)
{
    if (a<b)
    {
        return (a);
    }
    return (b);
}

/*********************************************************************************************************************************************************************************/
/* CROSSOVER */
/*********************************************************************************************************************************************************************************/
/* Function to cross two individuals */
void Nsga2::crossover (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
    if (nreal!=0)
    {
          realcross (parent1, parent2, child1, child2);
    }
    if (nbin!=0)
    {
        bincross (parent1, parent2, child1, child2);
    }
    return;
}

/* Routine for real variable SBX crossover */
void Nsga2::realcross (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
    int i;
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;
    if (rndGenerator->randomperc() <= pcross_real)
    {
        nrealcross++;
        for (i=0; i<nreal; i++)
        {
            if (rndGenerator->randomperc()<=0.5 )
            {
                if (fabs(parent1->xreal[i]-parent2->xreal[i]) > EPS)
                {
                    if (parent1->xreal[i] < parent2->xreal[i])
                    {
                        y1 = parent1->xreal[i];
                        y2 = parent2->xreal[i];
                    }
                    else
                    {
                        y1 = parent2->xreal[i];
                        y2 = parent1->xreal[i];
                    }
                    yl = min_realvar[i];
                    yu = max_realvar[i];
                    rand = rndGenerator->randomperc();
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(eta_c+1.0));
                    if (rand <= (1.0/alpha))
                    {
                        betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(eta_c+1.0));
                    if (rand <= (1.0/alpha))
                    {
                        betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                    if (c1<yl)
                        c1=yl;
                    if (c2<yl)
                        c2=yl;
                    if (c1>yu)
                        c1=yu;
                    if (c2>yu)
                        c2=yu;
                    if (rndGenerator->randomperc()<=0.5)
                    {
                        child1->xreal[i] = c2;
                        child2->xreal[i] = c1;
                    }
                    else
                    {
                        child1->xreal[i] = c1;
                        child2->xreal[i] = c2;
                    }
                }
                else
                {
                    child1->xreal[i] = parent1->xreal[i];
                    child2->xreal[i] = parent2->xreal[i];
                }
            }
            else
            {
                child1->xreal[i] = parent1->xreal[i];
                child2->xreal[i] = parent2->xreal[i];
            }
        }
    }
    else
    {
        for (i=0; i<nreal; i++)
        {
            child1->xreal[i] = parent1->xreal[i];
            child2->xreal[i] = parent2->xreal[i];
        }
    }
    return;
}


/* Routine for two point binary crossover */
void Nsga2::bincross (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
    int i, j;
    double rand;
    int temp, site1, site2;
    for (i=0; i<nbin; i++)
    {
        rand = rndGenerator->randomperc();
        if (rand <= pcross_bin)
        {
            nbincross++;
            site1 = rndGenerator->rnd(0,nbits[i]-1);
            site2 = rndGenerator->rnd(0,nbits[i]-1);
            if (site1 > site2)
            {
                temp = site1;
                site1 = site2;
                site2 = temp;
            }
            for (j=0; j<site1; j++)
            {
                child1->gene[i][j] = parent1->gene[i][j];
                child2->gene[i][j] = parent2->gene[i][j];
            }
            for (j=site1; j<site2; j++)
            {
                child1->gene[i][j] = parent2->gene[i][j];
                child2->gene[i][j] = parent1->gene[i][j];
            }
            for (j=site2; j<nbits[i]; j++)
            {
                child1->gene[i][j] = parent1->gene[i][j];
                child2->gene[i][j] = parent2->gene[i][j];
            }
        }
        else
        {
            for (j=0; j<nbits[i]; j++)
            {
                child1->gene[i][j] = parent1->gene[i][j];
                child2->gene[i][j] = parent2->gene[i][j];
            }
        }
    }
    return;
}

/*********************************************************************************************************************************************************************************/
/* CROWDIST */
/*********************************************************************************************************************************************************************************/
/* Routine to compute crowding distance based on ojbective function values when the population in in the form of a list */
void Nsga2::assign_crowding_distance_list (population *pop, list_ *lst, int front_size)
{
    int **obj_array;
    int *dist;
    int i, j;
    list_ *temp;
    temp = lst;
    if (front_size==1)
    {
        pop->ind[lst->index].crowd_dist = INF;
        return;
    }
    if (front_size==2)
    {
        pop->ind[lst->index].crowd_dist = INF;
        pop->ind[lst->child->index].crowd_dist = INF;
        return;
    }
    obj_array = (int **)malloc(nobj*sizeof(int*));
    dist = (int *)malloc(front_size*sizeof(int));
    for (i=0; i<nobj; i++)
    {
        obj_array[i] = (int *)malloc(front_size*sizeof(int));
    }
    for (j=0; j<front_size; j++)
    {
        dist[j] = temp->index;
        temp = temp->child;
    }
    assign_crowding_distance (pop, dist, obj_array, front_size);
    free (dist);
    for (i=0; i<nobj; i++)
    {
        free (obj_array[i]);
    }
    free (obj_array);
    return;
}

/* Routine to compute crowding distance based on objective function values when the population in in the form of an array */
void Nsga2::assign_crowding_distance_indices (population *pop, int c1, int c2)
{
    int **obj_array;
    int *dist;
    int i, j;
    int front_size;
    front_size = c2-c1+1;
    if (front_size==1)
    {
        pop->ind[c1].crowd_dist = INF;
        return;
    }
    if (front_size==2)
    {
        pop->ind[c1].crowd_dist = INF;
        pop->ind[c2].crowd_dist = INF;
        return;
    }
    obj_array = (int **)malloc(nobj*sizeof(int*));
    dist = (int *)malloc(front_size*sizeof(int));
    for (i=0; i<nobj; i++)
    {
        obj_array[i] = (int *)malloc(front_size*sizeof(int));
    }
    for (j=0; j<front_size; j++)
    {
        dist[j] = c1++;
    }
    assign_crowding_distance (pop, dist, obj_array, front_size);
    free (dist);
    for (i=0; i<nobj; i++)
    {
        free (obj_array[i]);
    }
    free (obj_array);
    return;
}

/* Routine to compute crowding distances */
void Nsga2::assign_crowding_distance (population *pop, int *dist, int **obj_array, int front_size)
{
    int i, j;
    for (i=0; i<nobj; i++)
    {
        for (j=0; j<front_size; j++)
        {
            obj_array[i][j] = dist[j];
        }
        quicksort_front_obj (pop, i, obj_array[i], front_size);
    }
    for (j=0; j<front_size; j++)
    {
        pop->ind[dist[j]].crowd_dist = 0.0;
    }
    for (i=0; i<nobj; i++)
    {
        pop->ind[obj_array[i][0]].crowd_dist = INF;
    }
    for (i=0; i<nobj; i++)
    {
        for (j=1; j<front_size-1; j++)
        {
            if (pop->ind[obj_array[i][j]].crowd_dist != INF)
            {
                if (pop->ind[obj_array[i][front_size-1]].obj[i] == pop->ind[obj_array[i][0]].obj[i])
                {
                    pop->ind[obj_array[i][j]].crowd_dist += 0.0;
                }
                else
                {
                    pop->ind[obj_array[i][j]].crowd_dist += (pop->ind[obj_array[i][j+1]].obj[i] - pop->ind[obj_array[i][j-1]].obj[i])/(pop->ind[obj_array[i][front_size-1]].obj[i] - pop->ind[obj_array[i][0]].obj[i]);
                }
            }
        }
    }
    for (j=0; j<front_size; j++)
    {
        if (pop->ind[dist[j]].crowd_dist != INF)
        {
            pop->ind[dist[j]].crowd_dist = (pop->ind[dist[j]].crowd_dist)/nobj;
        }
    }
    return;
}

/*********************************************************************************************************************************************************************************/
/* SORT */
/*********************************************************************************************************************************************************************************/
/* Randomized quick sort routine to sort a population based on a particular objective chosen */
void Nsga2::quicksort_front_obj(population *pop, int objcount, int obj_array[], int obj_array_size)
{
    q_sort_front_obj (pop, objcount, obj_array, 0, obj_array_size-1);
    return;
}

/* Actual implementation of the randomized quick sort used to sort a population based on a particular objective chosen */
void Nsga2::q_sort_front_obj(population *pop, int objcount, int obj_array[], int left, int right)
{
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rndGenerator->rnd (left, right);
        temp = obj_array[right];
        obj_array[right] = obj_array[index];
        obj_array[index] = temp;
        pivot = pop->ind[obj_array[right]].obj[objcount];
        i = left-1;
        for (j=left; j<right; j++)
        {
            if (pop->ind[obj_array[j]].obj[objcount] <= pivot)
            {
                i+=1;
                temp = obj_array[j];
                obj_array[j] = obj_array[i];
                obj_array[i] = temp;
            }
        }
        index=i+1;
        temp = obj_array[index];
        obj_array[index] = obj_array[right];
        obj_array[right] = temp;
        q_sort_front_obj (pop, objcount, obj_array, left, index-1);
        q_sort_front_obj (pop, objcount, obj_array, index+1, right);
    }
    return;
}

/* Randomized quick sort routine to sort a population based on crowding distance */
void Nsga2::quicksort_dist(population *pop, int *dist, int front_size)
{
    q_sort_dist (pop, dist, 0, front_size-1);
    return;
}

/* Actual implementation of the randomized quick sort used to sort a population based on crowding distance */
void Nsga2::q_sort_dist(population *pop, int *dist, int left, int right)
{
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rndGenerator->rnd (left, right);
        temp = dist[right];
        dist[right] = dist[index];
        dist[index] = temp;
        pivot = pop->ind[dist[right]].crowd_dist;
        i = left-1;
        for (j=left; j<right; j++)
        {
            if (pop->ind[dist[j]].crowd_dist <= pivot)
            {
                i+=1;
                temp = dist[j];
                dist[j] = dist[i];
                dist[i] = temp;
            }
        }
        index=i+1;
        temp = dist[index];
        dist[index] = dist[right];
        dist[right] = temp;
        q_sort_dist (pop, dist, left, index-1);
        q_sort_dist (pop, dist, index+1, right);
    }
    return;
}

/*********************************************************************************************************************************************************************************/
/* DECODE */
/*********************************************************************************************************************************************************************************/
/* Function to decode a population to find out the binary variable values based on its bit pattern */
void Nsga2::decode_pop (population *pop)
{
    int i;
    if (nbin!=0)
    {
        for (i=0; i<popsize; i++)
        {
            decode_ind (&(pop->ind[i]));
        }
    }
    return;
}

/* Function to decode an individual to find out the binary variable values based on its bit pattern */
void Nsga2::decode_ind (individual *ind)
{
    int j, k;
    int sum;
    if (nbin!=0)
    {
        for (j=0; j<nbin; j++)
        {
            sum=0;
            for (k=0; k<nbits[j]; k++)
            {
                if (ind->gene[j][k]==1)
                {
                    sum += pow(2,nbits[j]-1-k);
                }
            }
            ind->xbin[j] = min_binvar[j] + (double)sum*(max_binvar[j] - min_binvar[j])/(double)(pow(2,nbits[j])-1);
        }
    }
    return;
}

/*********************************************************************************************************************************************************************************/
/* DISPLAY */
/*********************************************************************************************************************************************************************************/
/* Function to display the current population for the subsequent generation */
void Nsga2::onthefly_display (population *pop, FILE *gp, int ii)
{
    int i;
    int flag;
    FILE *fpt;
    fpt = fopen("plot.out","w");
    flag = 0;
    for (i=0; i<popsize; i++)
    {
        if (pop->ind[i].constr_violation==0)
        {
            if (choice!=3)
                fprintf(fpt,"%e\t%e\n",pop->ind[i].obj[obj1-1],pop->ind[i].obj[obj2-1]);
            else
                fprintf(fpt,"%e\t%e\t%e\n",pop->ind[i].obj[obj1-1],pop->ind[i].obj[obj2-1],pop->ind[i].obj[obj3-1]);
            fflush(fpt);
            flag = 1;
        }
    }
    if (flag==0)
    {
        printf("\n No feasible soln in this pop, hence no display");
    }
    else
    {
        if (choice!=3)
            fprintf(gp,"set title 'Generation #%d'\n unset key\n plot 'plot.out' w points pointtype 6 pointsize 1\n",ii);
        else
            fprintf(gp,"set title 'Generation #%d'\n set view %d,%d\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n",ii,angle1,angle2);
        fflush(gp);
    }
    fclose(fpt);
    return;
}

/*********************************************************************************************************************************************************************************/
/* DOMINANCE */
/*********************************************************************************************************************************************************************************/
/* Routine for usual non-domination checking
   It will return the following values
   1 if a dominates b
   -1 if b dominates a
   0 if both a and b are non-dominated */

int Nsga2::check_dominance (individual *a, individual *b)
{
    int i;
    int flag1;
    int flag2;
    flag1 = 0;
    flag2 = 0;
	
	/* Aniadido para poder guardar todas las RNAinverse folding sequences con similaridad igual a 1 */
	 /*if ( a->obj[1] == b->obj[1] ) return 0; */
	
    if (a->constr_violation<0 && b->constr_violation<0)
    {
        if (a->constr_violation > b->constr_violation)
        {
            return (1);
        }
        else
        {
            if (a->constr_violation < b->constr_violation)
            {
                return (-1);
            }
            else
            {
                for (i=0; i<nobj; i++)
                {
                    if (a->obj[i] < b->obj[i])
                    {
                        flag1 = 1;

                    }
                    else
                    {
                        if (a->obj[i] > b->obj[i])
                        {
                            flag2 = 1;
                        }
                    }
                }
                if (flag1==1 && flag2==0)
                {
                    return (1);
                }
                else
                {
                    if (flag1==0 && flag2==1)
                    {
                        return (-1);
                    }
                    else
                    {
                        return (0);
                    }
                }
                //return (0);
            }
        }
    }
    else
    {
        if (a->constr_violation < 0 && b->constr_violation == 0)
        {
            return (-1);
        }
        else
        {
            if (a->constr_violation == 0 && b->constr_violation <0)
            {
                return (1);
            }
            else
            {
                for (i=0; i<nobj; i++)
                {
                    if (a->obj[i] < b->obj[i])
                    {
                        flag1 = 1;

                    }
                    else
                    {
                        if (a->obj[i] > b->obj[i])
                        {
                            flag2 = 1;
                        }
                    }
                }
                if (flag1==1 && flag2==0)
                {
                    return (1);
                }
                else
                {
                    if (flag1==0 && flag2==1)
                    {
                        return (-1);
                    }
                    else
                    {
                        return (0);
                    }
                }
            }
        }
    }
}

/*********************************************************************************************************************************************************************************/
/* EVAL */
/*********************************************************************************************************************************************************************************/
/* Routine to evaluate objective function values and constraints for a population */
void Nsga2::evaluate_pop (population *pop)
{
    int i;
    for (i=0; i<popsize; i++)
    {
        evaluate_ind (&(pop->ind[i]));
	if (minSim == 0 && problem->getID() == 0 && !found)
	{
	  if (solvingFlag) printf("#");
	  report_ind(&pop->ind[i]);
          printf(
              "\n@FOUND: %f secs\t", (((double) (clock() - runtime_start)) / CLOCKS_PER_SEC) );
              //cout << "@ORACLE: " << getSizeOracle() << " sequences found" << endl;
          printf("\n");
	  cout << solvingFlag << endl;
	  found = true;
	  if (solvingFlag) exit(1);
	}
    }
    return;
}

/* Routine to evaluate objective function values and constraints for an individual */
void Nsga2::evaluate_ind (individual *ind)
{
    int j;
    test_problem (ind->xreal, ind->rnaMatches, ind->mfe, ind->sequence, ind->xbin, ind->gene, ind->obj, ind->constr);
//    nEvals++;
    if (ncon==0)
    {
        ind->constr_violation = 0.0;
    }
    else
    {
        ind->constr_violation = 0.0;
        for (j=0; j<ncon; j++)
        {
            if (ind->constr[j]<0.0)
            {
                ind->constr_violation += ind->constr[j];
            }
        }
    }
    return;
}

/*********************************************************************************************************************************************************************************/
/* FILLNDS */
/*********************************************************************************************************************************************************************************/
/* Routine to perform non-dominated sorting */
void Nsga2::fill_nondominated_sort (population *mixed_pop, population *new_pop)
{
    int flag;
    int i, j;
    int end;
    int front_size;
    int archieve_size;
    int rank=1;
    list_ *pool;
    list_ *elite;
    list_ *temp1, *temp2;
    pool = (list_ *)malloc(sizeof(list_));
    elite = (list_ *)malloc(sizeof(list_));
    front_size = 0;
    archieve_size=0;
    pool->index = -1;
    pool->parent = NULL;
    pool->child = NULL;
    elite->index = -1;
    elite->parent = NULL;
    elite->child = NULL;
    temp1 = pool;
    for (i=0; i<2*popsize; i++)
    {
        insert (temp1,i);
        temp1 = temp1->child;
    }
    i=0;
    do
    {
        temp1 = pool->child;
        insert (elite, temp1->index);
        front_size = 1;
        temp2 = elite->child;
        temp1 = del (temp1);
        temp1 = temp1->child;
        do
        {
            temp2 = elite->child;
            if (temp1==NULL)
            {
                break;
            }
            do
            {
                end = 0;
                flag = check_dominance (&(mixed_pop->ind[temp1->index]), &(mixed_pop->ind[temp2->index]));
                if (flag == 1)
                {
                    insert (pool, temp2->index);
                    temp2 = del (temp2);
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0)
                {
                    temp2 = temp2->child;
                }
                if (flag == -1)
                {
                    end = 1;
                }
            }
            while (end!=1 && temp2!=NULL);
            if (flag == 0 || flag == 1)
            {
                insert (elite, temp1->index);
                front_size++;
                temp1 = del (temp1);
            }
            temp1 = temp1->child;
        }
        while (temp1 != NULL);
        temp2 = elite->child;
        j=i;
        if ( (archieve_size+front_size) <= popsize)
        {
            do
            {
                copy_ind (&mixed_pop->ind[temp2->index], &new_pop->ind[i]);
                new_pop->ind[i].rank = rank;
                archieve_size+=1;
                temp2 = temp2->child;
                i+=1;
            }
            while (temp2 != NULL);
            assign_crowding_distance_indices (new_pop, j, i-1);
            rank+=1;
        }
        else
        {
            crowding_fill (mixed_pop, new_pop, i, front_size, elite);
            archieve_size = popsize;
            for (j=i; j<popsize; j++)
            {
                new_pop->ind[j].rank = rank;
            }
        }
        temp2 = elite->child;
        do
        {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
        while (elite->child !=NULL);
    }
    while (archieve_size < popsize);
    while (pool!=NULL)
    {
        temp1 = pool;
        pool = pool->child;
        free (temp1);
    }
    while (elite!=NULL)
    {
        temp1 = elite;
        elite = elite->child;
        free (temp1);
    }
    return;
}

/* Routine to fill a population with individuals in the decreasing order of crowding distance */
void Nsga2::crowding_fill (population *mixed_pop, population *new_pop, int count, int front_size, list_ *elite)
{
    int *dist;
    list_ *temp;
    int i, j;
    int equal=0;
    assign_crowding_distance_list (mixed_pop, elite->child, front_size);
    dist = (int *)malloc(front_size*sizeof(int));
    temp = elite->child;
    for (j=0; j<front_size; j++)
    {
        dist[j] = temp->index;
        temp = temp->child;
    }
    quicksort_dist (mixed_pop, dist, front_size);
    for (i=count, j=front_size-1; i<popsize; i++, j--)
    {
        copy_ind(&mixed_pop->ind[dist[j]], &new_pop->ind[i]);
    }
    free (dist);
    return;
}

/*********************************************************************************************************************************************************************************/
/* LIST */
/*********************************************************************************************************************************************************************************/
/* A custom doubly linked list implemenation */

/* Insert an element X into the list at location specified by NODE */
void Nsga2::insert (list_ *node, int x)
{
    list_ *temp;
    if (node==NULL)
    {
        printf("\n Error!! asked to enter after a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = (list_ *)malloc(sizeof(list_));
    temp->index = x;
    temp->child = node->child;
    temp->parent = node;
    if (node->child != NULL)
    {
        node->child->parent = temp;
    }
    node->child = temp;
    return;
}

/* Delete the node NODE from the list */
list_* Nsga2::del (list_ *node)
{
    list_ *temp;
    if (node==NULL)
    {
        printf("\n Error!! asked to delete a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = node->parent;
    temp->child = node->child;
    if (temp->child!=NULL)
    {
        temp->child->parent = temp;
    }
    free (node);
    return (temp);
}

/*********************************************************************************************************************************************************************************/
/* INITIALIZE */
/*********************************************************************************************************************************************************************************/
/* Function to initialize a population randomly */
void Nsga2::initialize_pop (population *pop, vector<string>* ssPop)
{
    int i;
    if (ssPop != NULL)
    {
      for (i=0; i < ssPop->size(); i++)
      {
	if (i >= popsize) break;
        cout << "**** " << ssPop->at(i) << " ****" << endl;
        read_ind(&(pop->ind[i]), ssPop->at(i));
      }
      for (i=ssPop->size(); i<popsize; i++)
	initialize_ind(&pop->ind[i], i+1);
    }
    else
    {
      for (i=0; i<popsize; i++)
	initialize_ind(&pop->ind[i], i+1);
    }
//    report_pop_xreal(pop);
    return;
}

void Nsga2::initialize_ind (individual *ind, int i)
{
        if(i%10 == 0)
            initialize_ind0 (ind);
        else
          if(i%5==0)
            initialize_ind3 (ind);
          else
            if(i%2==0)
              initialize_ind1 (ind);
            else
              initialize_ind2 (ind);
}

void Nsga2::initialize_ind (individual *ind)
{
      double lb, ub, v, randomN;
      randomN = rndGenerator->rndreal(0.0, 1.0);

      if (randomN < 0.2) {
          lb = 0.66;
          ub = 1.00;
          v = 0.9;
      }
      else {
          if (randomN < 0.6) {
              lb = 0.00;
              ub = 0.32;
              v = 0.6;
           } else {
              lb = 0.33;
              ub = 0.65;
              v = 0.2;
           }
      }
      for (int n = 0; n < nLoops; n++)
          ind->xreal[n] = rndGenerator->rndreal (lb, ub);

      for (int n = 0; n < nPoints; n++)
          ind->xreal[nLoops + n] = v;

      for (int n = 0; n < nSS; n++)
          ind->xreal[nLoops + nPoints + n] = rndGenerator->rndreal (-0.25, 1.0); // If negative, dismiss subsequence
}

void Nsga2::initialize_ind0 (individual *ind)
{
        int i;
        for (i=0; i<nLoops; i++){
          ind->xreal[i] = rndGenerator->rndreal (0, 1.0);
        }
        for(i=0; i<nPoints; i++){
	  ind->xreal[nLoops+i] = rndGenerator->rndreal(0.0, 1.0);
        }
        for (i = 0; i < nSS; i++)
          ind->xreal[nLoops + nPoints + i] = rndGenerator->rndreal (-0.25, 1.0); // If negative, dismiss subsequence
}

void Nsga2::initialize_ind1 (individual *ind)
{
        int i;
        for (i=0; i<nLoops; i++){ // GC o CG
          ind->xreal[i] = rndGenerator->rndreal (0, probabilitiesLoops[0]+probabilitiesLoops[1]); //0.32);
        }
        for(i=0; i<nPoints; i++){ // A
	  ind->xreal[nLoops+i] = rndGenerator->rndreal(probabilitiesUnpaired[0]+probabilitiesUnpaired[1], probabilitiesUnpaired[0] + probabilitiesUnpaired[1] + probabilitiesUnpaired[2]); //0.6;
        }
        for (i = 0; i < nSS; i++)
          ind->xreal[nLoops + nPoints + i] = rndGenerator->rndreal (-0.25, 1.0); // If negative, dismiss subsequence
}

void Nsga2::initialize_ind2 (individual *ind)
{
        int i;
        for (i=0; i<nLoops; i++){ // UA o AU
          ind->xreal[i] = rndGenerator->rndreal (probabilitiesLoops[0]+probabilitiesLoops[1], probabilitiesLoops[0]+probabilitiesLoops[1]+probabilitiesLoops[2]+probabilitiesLoops[3]);//0.33, 0.65);
        }
        for(i=0; i<nPoints; i++){ // C
          ind->xreal[nLoops+i] = rndGenerator->rndreal(0, probabilitiesUnpaired[0]); //0.2;
        }
        for (i = 0; i < nSS; i++)
          ind->xreal[nLoops + nPoints + i] = rndGenerator->rndreal (-0.25, 1.0); // If negative, dismiss subsequence
}

void Nsga2::initialize_ind3 (individual *ind)
{
        int i;
        for (i=0; i<nLoops; i++){ // GU o UG
          ind->xreal[i] = rndGenerator->rndreal (probabilitiesLoops[0]+probabilitiesLoops[1]+probabilitiesLoops[2]+probabilitiesLoops[3], 1.00); //0.66, 1.00);
        }
        for(i=0; i<nPoints; i++){ // U
          ind->xreal[nLoops+i] = rndGenerator->rndreal(probabilitiesUnpaired[0]+probabilitiesUnpaired[1]+probabilitiesUnpaired[2], 1.0);//0.9;
        }
        for (i = 0; i < nSS; i++)
          ind->xreal[nLoops + nPoints + i] = rndGenerator->rndreal (-0.25, 1.0); // If negative, dismiss subsequence
}

/**
 * Inicialización NO ALEATORIA:
 * inicializa el cromosoma (x_var) del individuo de acuerdo al string de entrada
 */
void Nsga2::read_ind(individual *ind, string seq_)
{

  double pCG = probabilitiesLoops[0] / 2;
  double pGC = probabilitiesLoops[0] + probabilitiesLoops[1] / 2;
  double pAU = probabilitiesLoops[0] + probabilitiesLoops[1] + probabilitiesLoops[2] / 2;
  double pUA = probabilitiesLoops[0] + probabilitiesLoops[1] + probabilitiesLoops[2] + probabilitiesLoops[3] / 2;
  double pGU = probabilitiesLoops[0] + probabilitiesLoops[1] + probabilitiesLoops[2] + probabilitiesLoops[3] + probabilitiesLoops[4] / 2;
  double pUG = probabilitiesLoops[0] + probabilitiesLoops[1] + probabilitiesLoops[2] + probabilitiesLoops[3] + probabilitiesLoops[4] + probabilitiesLoops[5] / 2;

  double pC = probabilitiesUnpaired[0] / 2;
  double pG = probabilitiesUnpaired[0] + probabilitiesUnpaired[1] / 2;
  double pA = probabilitiesUnpaired[0] + probabilitiesUnpaired[1] + probabilitiesUnpaired[2] / 2;
  double pU = probabilitiesUnpaired[0] + probabilitiesUnpaired[1] + probabilitiesUnpaired[2] + probabilitiesUnpaired[3] / 2;

  for (int n = 0; n < nLoops; n++)
  {
    if      ( seq_[initLoopTargetRNA[n]] == 'C' && seq_[endLoopTargetRNA[n]] == 'G' ) ind->xreal[n] = pCG;
    else if ( seq_[initLoopTargetRNA[n]] == 'G' && seq_[endLoopTargetRNA[n]] == 'C' ) ind->xreal[n] = pGC;
    else if ( seq_[initLoopTargetRNA[n]] == 'A' && seq_[endLoopTargetRNA[n]] == 'U' ) ind->xreal[n] = pAU;
    else if ( seq_[initLoopTargetRNA[n]] == 'U' && seq_[endLoopTargetRNA[n]] == 'A' ) ind->xreal[n] = pUA;
    else if ( seq_[initLoopTargetRNA[n]] == 'G' && seq_[endLoopTargetRNA[n]] == 'U' ) ind->xreal[n] = pGU;
    else if ( seq_[initLoopTargetRNA[n]] == 'U' && seq_[endLoopTargetRNA[n]] == 'G' ) ind->xreal[n] = pUG;
    else cout << "ERROR 1" << endl;
  }

  for (int n = 0; n < nPoints; n++)
  {
    if      ( seq_[pointTargetRNA[n]] == 'A' ) ind->xreal[nLoops + n] = pA;
    else if ( seq_[pointTargetRNA[n]] == 'C' ) ind->xreal[nLoops + n] = pC;
    else if ( seq_[pointTargetRNA[n]] == 'G' ) ind->xreal[nLoops + n] = pG;
    else if ( seq_[pointTargetRNA[n]] == 'U' ) ind->xreal[nLoops + n] = pU;
    else cout << "ERROR 2" << endl;
  }

  for (int i = 0; i < nSS; i++)
  {
    cout << ssTargetRNA[i] << endl;
    //cout << ssInitTargetRNA[i] << endl;
    Problem* p = problem->getChild(ssTargetRNA[i]);
    string s = seq_.substr(ssInitTargetRNA[i], p->getTargetRNA().size());
    //cout << s << endl << endl;
    int k;
    for (k = 0; k < p->getSequences()->size(); k++)
    {
      //cout << p->getSequenceAt(k) << endl;
      //cout << s << endl << endl;
      if (p->getSequenceAt(k).compare(s) == 0) break;
    }
    cout << k << " / " << p->getSequences()->size() << endl << endl;
    if (k == p->getSequences()->size())
      ind->xreal[nLoops + nPoints + i] = rndGenerator->rndreal (0.0, 1.0);
    else
      ind->xreal[nLoops + nPoints + i] = (k*1.0) / p->getSequences()->size();
  }
/*
 cout << "----" << endl;
  for (int i=0; i<nreal; i++)
  {
   if (i==nLoops) cout << endl;
   if (i==nLoops+nPoints) cout << endl;
   cout << ind->xreal[i] << endl;
  }
 cout << "----" << endl;
*/
}

/*********************************************************************************************************************************************************************************/
/* MERGE */
/*********************************************************************************************************************************************************************************/
/* Routine to merge two populations into one */
void Nsga2::merge(population *pop1, population *pop2, population *pop3)
{
    int i, k;
    for (i=0; i<popsize; i++)
    {
        copy_ind (&(pop1->ind[i]), &(pop3->ind[i]));
    }
    for (i=0, k=popsize; i<popsize; i++, k++)
    {
        copy_ind (&(pop2->ind[i]), &(pop3->ind[k]));
    }
    return;
}

/* Routine to copy an individual 'ind1' into another individual 'ind2' */
void Nsga2::copy_ind (individual *ind1, individual *ind2)
{
    int i, j;
    ind2->rank = ind1->rank;
    ind2->constr_violation = ind1->constr_violation;
    ind2->crowd_dist = ind1->crowd_dist;
    strcpy(ind2->mfe,ind1->mfe);
    strcpy(ind2->sequence,ind1->sequence);

    if (nreal!=0)
    {
        for (i=0; i<nreal; i++)
        {
            ind2->xreal[i] = ind1->xreal[i];
	    ind2->rnaMatches[i] = ind1->rnaMatches[i];
	}
/*			ind2->rnaMatches[i] = ind1->rnaMatches[i];

        }
	for (i=nreal; i<strlen(targetRNA); i++){
		ind2->rnaMatches[i] = ind1->rnaMatches[i];
	}		*/
    }
    if (nbin!=0)
    {
        for (i=0; i<nbin; i++)
        {
            ind2->xbin[i] = ind1->xbin[i];
            for (j=0; j<nbits[i]; j++)
            {
                ind2->gene[i][j] = ind1->gene[i][j];
            }
        }
    }
    for (i=0; i<nobj; i++)
    {
        ind2->obj[i] = ind1->obj[i];
    }
    if (ncon!=0)
    {
        for (i=0; i<ncon; i++)
        {
            ind2->constr[i] = ind1->constr[i];
        }
    }
    return;
}

/*********************************************************************************************************************************************************************************/
/* MUTATION */
/*********************************************************************************************************************************************************************************/
/* Function to perform mutation in a population */
void Nsga2::mutation_pop (population *pop, int stagnation)
{
    int i;
    for (i=0; i<popsize; i++)
    {
        //double rndN = rndGenerator->rndreal(0.0, 1.0);
	//if (rndN < (0.3 + stagnation*0.1) && minSim != 0)// || (rndN < (0.3 + stagnation*0.01) && minSim == 0) )
        //  initialize_ind(&pop->ind[i], i+1);
        //else
          mutation_ind(&(pop->ind[i]));
    }
    return;
}

/* Function to perform mutation of an individual */
void Nsga2::mutation_ind (individual *ind)
{
    if (nreal!=0)
    {
          real_mutate_ind(ind);
    }
    if (nbin!=0)
    {
        bin_mutate_ind(ind);
    }
    return;
}

/* Routine for binary mutation of an individual */
void Nsga2::bin_mutate_ind (individual *ind)
{
    int j, k;
    double prob;
    for (j=0; j<nbin; j++)
    {
        for (k=0; k<nbits[j]; k++)
        {
            prob = rndGenerator->randomperc();
            if (prob <=pmut_bin)
            {
                if (ind->gene[j][k] == 0)
                {
                    ind->gene[j][k] = 1;
                }
                else
                {
                    ind->gene[j][k] = 0;
                }
                nbinmut+=1;
            }
        }
    }
    return;
}

/* Routine for real polynomial mutation of an individual */
void Nsga2::real_mutate_ind (individual *ind)
{
    int j;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;

//    cout << "Failed: ";
//    for (int i=0; i<len; i++) cout << "(" << i << ", " << failedPos[i] << ") ";
//    cout << endl << endl;
/*
    std::priority_queue<std::pair<int, int>> q;
    for (int i=0; i<len; i++)
      q.push(std::pair<int,int>(failedPos[i],i));

    // MUTATE THE FIVE MOST FAILED POSITIONS
    int k = 5; // number of indices we need
    for (int i = 0; i < k; ++i) {
      int ki = q.top().second;
//      std::cout << "index[" << i << "] = " << ki << std::endl;
//          cout << "ANTES: " << ind->xreal[j] << endl;
      int j = -1;
      for (int k=0; k<nLoops && j==-1; k++)
        if (initLoopTargetRNA[k] == ki || endLoopTargetRNA[k] == ki)
	  j = k;
      for (int k=0; k<nPoints && j == -1; k++)
        if (pointTargetRNA[k] == ki)
          j= nLoops+k;
      for (int k=0; k<nSS && j==-1; k++)
        if (ki >= ssInitTargetRNA[k] && ki <= ssEndTargetRNA[k])
          j = nLoops+nPoints+k;
//      cout << j << endl;
      if (j != -1)
      {
        y = ind->xreal[j];
        yl = min_realvar[j];
        yu = max_realvar[j];
        delta1 = (y-yl)/(yu-yl);
        delta2 = (yu-y)/(yu-yl);
        rnd = randomperc();
        mut_pow = 1.0/(eta_m+1.0);
        if (rnd <= 0.5)
        {
           xy = 1.0-delta1;
           val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
           deltaq =  pow(val,mut_pow) - 1.0;
        }
        else
        {
          xy = 1.0-delta2;
          val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
          deltaq = 1.0 - (pow(val,mut_pow));
        }
        y = y + deltaq*(yu-yl);
        if (y<yl)
           y = yl;
        if (y>yu)
        y = yu;
        ind->xreal[j] = y;
//          cout << "DESPUÉS: " << ind->xreal[j] << endl;
        nrealmut+=1;
      }
      q.pop();
    }
//    cout << endl << endl;
*/
//    cout << "MUTACIÓN" << endl;
    for (j=0; j<nreal; j++)
    {
        if ( (j < (nLoops + nPoints) && rndGenerator->randomperc() <= pmut_real) || (j >= (nLoops + nPoints) && rndGenerator->randomperc() <= pmut_real_subs) )
        {
//	    cout << "ANTES: " << ind->xreal[j] << endl;
            y = ind->xreal[j];
            yl = min_realvar[j];
            yu = max_realvar[j];
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rnd = rndGenerator->randomperc();
            mut_pow = 1.0/(eta_m+1.0);
            if (rnd <= 0.5)
            {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            }
            else
            {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl)
                y = yl;
            if (y>yu)
                y = yu;
            ind->xreal[j] = y;
//	    cout << "DESPUÉS: " << ind->xreal[j] << endl;
            nrealmut+=1;
        }
    }
//    cout << endl;
    return;
}

int Nsga2::getEditDistance(std::string first, std::string second)
{
    int m = first.length();
    int n = second.length();
 
    int T[m + 1][n + 1];
    T[0][0] = 0;
    for (int i = 1; i <= m; i++) {
        T[i][0] = i;
    }
 
    for (int j = 1; j <= n; j++) {
        T[0][j] = j;
    }
 
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int weight = first[i - 1] == second[j - 1] ? 0: 1;
            T[i][j] = std::min(std::min(T[i-1][j] + 1, T[i][j-1] + 1), T[i-1][j-1] + weight);
        }
    }
 
/*    for (int i=0; i<=m; i++)
    {
      for (int j=0; j<=n; j++)
        std::cout << T[i][j] << " ";
      std::cout << std::endl;
    }*/

    return T[m][n];
}

double Nsga2::findStringSimilarity(std::string first, std::string second) {
    double max_length = std::max(first.length(), second.length());
    if (max_length > 0) {
        return (max_length - getEditDistance(first, second)) / max_length;
    }
    return 1.0;
}

/*********************************************************************************************************************************************************************************/
/* PROBLEMDEF */
/*********************************************************************************************************************************************************************************/
void Nsga2::test_problem (double *xreal, int *rnaMatches, char* mfe, char* sequence, double *xbin, int **gene, double *obj, double *constr)
{
  strcpy(seq,"");
  strcpy(sequence,"");
  strcpy(mfe,"");
  seq[len]='\0';
  sequence[len]='\0';
  mfe[len]='\0';
  int i,j, c;
  double value;
  int len = strlen(targetRNA);
  for (i=0; i<nLoops; i++){
    value = xreal[i];
//cout << "value: " << value << endl;
    rnaMatches[i]=0;
	/* CG (0 to 0.16), GC (0.16 to 0.33), AU (0.33 to 0.49), UA (0.49 to 0.66), GU (0.66 to 0.83), or UG (0.83 to 1.00) */
	/*if ( (value >= 0) && (value < 0.16) )*/
	double p = 0;
	double p1 = p + probabilitiesLoops[0];
//cout << "p: " << p << " p1: " << p1 << endl;
	if ( (value >= p) && (value < p1) )
	{
		seq[initLoopTargetRNA[i]] = 'C';
		seq[endLoopTargetRNA[i]]  = 'G';
	}
	p = p + probabilitiesLoops[0];
	p1 = p + probabilitiesLoops[1];
//cout << "p: " << p << " p1: " << p1 << endl;
	/*if ( (value >= 0.16) && (value < 0.33) )*/
	if ( (value >= p) && (value < p1) )
	{
		seq[initLoopTargetRNA[i]] = 'G';
		seq[endLoopTargetRNA[i]]  = 'C';	
	}	
        p = p + probabilitiesLoops[1];
        p1 = p + probabilitiesLoops[2];
//cout << "p: " << p << " p1: " << p1 << endl;
	/*if ( (value >= 0.33) && (value < 0.49) )*/
	if ( (value >= p) && (value < p1) )
	{
		seq[initLoopTargetRNA[i]] = 'A';
		seq[endLoopTargetRNA[i]]  = 'U';	
	}	
        p = p + probabilitiesLoops[2];
        p1 = p + probabilitiesLoops[3];
//cout << "p: " << p << " p1: " << p1 << endl;
	/*if ( (value >= 0.49) && (value < 0.66) )*/
	if ( (value >= p) && (value < p1) )
	{
		seq[initLoopTargetRNA[i]] = 'U';
		seq[endLoopTargetRNA[i]]  = 'A';	
	}		
        p = p + probabilitiesLoops[3];
        p1 = p + probabilitiesLoops[4];
//cout << "p: " << p << " p1: " << p1 << endl;
	/*if ( (value >= 0.7) && (value < 0.83) )*/
	if ( (value >= p) && (value < p1) )
	{
		seq[initLoopTargetRNA[i]] = 'G';
		seq[endLoopTargetRNA[i]]  = 'U';	
	}
        p = p + probabilitiesLoops[4];
        p1 = p + probabilitiesLoops[5];
//cout << "p: " << p << " p1: " << p1 << endl;
	/*if ( (value >= 0.83) && (value <= 1.00) )*/
	if ( (value >= p) && (value <= p1) )
	{
		seq[initLoopTargetRNA[i]] = 'U';
		seq[endLoopTargetRNA[i]]  = 'G';	
	}
//cout << "[" << seq[initLoopTargetRNA[i]] << "]" << "[" << seq[endLoopTargetRNA[i]] << "]" << endl;
  }
  for(i=0; i<nPoints; i++){
    value = xreal[nLoops+i];
    rnaMatches[nLoops+i]=0;
//cout << "value: " << value << endl;
        double p = 0; 
        double p1 = p + probabilitiesUnpaired[0];
//cout << "p: " << p << " p1: " << p1 << endl;
	if ( (value >= p) && (value < p1) )
	{
		seq[pointTargetRNA[i]] = 'C';
	}
        p = p + probabilitiesUnpaired[0];
        p1 = p + probabilitiesUnpaired[1];
//cout << "p: " << p << " p1: " << p1 << endl;
	if ( (value >= p) && (value < p1) )
	{
		seq[pointTargetRNA[i]] = 'G';
	}
        p = p + probabilitiesUnpaired[1];
        p1 = p + probabilitiesUnpaired[2];
//cout << "p: " << p << " p1: " << p1 << endl;
	if ( (value >= p) && (value < p1) )
	{
		seq[pointTargetRNA[i]] = 'A';
	}
        p = p + probabilitiesUnpaired[2];
        p1 = p + probabilitiesUnpaired[3];
//cout << "p: " << p << " p1: " << p1 << endl;
	if ( (value >= p) && (value <= p1) )
	{
		seq[pointTargetRNA[i]] = 'U';
	}
//cout << "[" <<seq[pointTargetRNA[i]] << "]" <<endl;
  }
/*if(gSS != NULL)
{
            for (int i=0; i < gSS->n; i++)
            {
                cout << gSS->substr[i].structure << endl;
                for (int j=0; j < gSS->substr[i].sequences.size(); j++)
                {
                  cout << "\t" << gSS->substr[i].sequences.at(j) << endl;
                }
                cout << endl;
            }
}*/
  int k;
  for (int c1=0; c1 < nSS; c1++)
  {
    if (xreal[nLoops + nPoints + c1] >= 0) // If negative value, dismiss subsequence found
    {
      k = (int) (xreal[nLoops + nPoints + c1] * (problem->getChild(ssTargetRNA[c1])->getSequencesSize() - 1));// (gSS->substr[ssTargetRNA[c1]].sequences.size() - 1));
      //cout << k << "\t" << ssTargetRNA[c1] << "\t" << nSS << endl;
      //cout << "@@" << problem->getChildrenAt(ssTargetRNA[c1])->getSequenceAt(k) << "@@" << endl;
//    cout << "@@" << gSS->substr[ssTargetRNA[c1]].sequences[k] << "@@" << endl;
      for (int s1=0; s1 < problem->getChild(ssTargetRNA[c1])->getSequenceAt(k).size(); s1++)
        seq[s1+ssInitTargetRNA[c1]] = problem->getChild(ssTargetRNA[c1])->getSequenceAt(k).at(s1);
    }
  }
//  if (nSS > 0) exit(1);
  seq[len]='\0';
  strcpy(sequence, seq);
  //cout << "****  " << sequence << "  ****" << endl;
  //cout << "****  " << seq << "  ****" <<  endl << endl;
  itr = oracle.find(sequence);
  if ( itr == oracle.end() )
  {
	double sim1, sim2, sim3;
	nEvals++;
	/*-------------------------------------------------------------------------------------------------------------------------*/
	/* get a vrna_fold_compound with default settings */
	vrna_fold_compound_t *vc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_MFE);
	/* call MFE function */
	double _mfe = (double) vrna_mfe(vc, mfe);
	//cout << seq << endl;
	//cout << mfe << endl << endl;
        //strcpy(sequence, seq);
	c=0;
        for(i=0; i<len; i++){
          if(mfe[i] != targetRNA[i])
	  {
//	    failedPos[i]++;
	    c++;
	  }
        }
	int c1 = c;
        double distance =  1.0 - ( (c*1.0)/(len) ) ;
        constr[0] = distance - 1.0;
	sim1 = distance - 1.0;

//	cout << targetRNA << endl;
//	cout << mfe << endl;
        sim2 = - (1.0 - findStringSimilarity(targetRNA, mfe) );
//	cout << " --> " << constr[0] << endl;
//        cout << " --> " << distance - 1.0 << endl << endl;
	/***********************************************************************************************************************************************************************************************************************/
//	cout << "MFE: \t\t" << mfe << endl;
//	cout << "TARGET: \t" << targetRNA << endl;
	int d = 0;
	for (i = 0; i < len; i++)
	{
	  if      (mfe[i] == '.' && targetRNA[i] != '.'){ d++; }
	  else if (mfe[i] == '(' && targetRNA[i] != '('){ d++; }
	  else if (mfe[i] == '(' && targetRNA[i] == '(')
	  {
	    c = 0;
	    int endMFE;
	    for (int j = i; j < len; j++)
	    {
	      if      (mfe[j] == '(') c++;
	      else if (mfe[j] == ')') c--;
	      else if ( c <= 0 && mfe[j] == ')' )
	      {
 		endMFE = j;
		break;
	      }
	    }
	    c=0;
	    int endTarget;
            for (int j = i; j < len; j++)
            {
              if      (targetRNA[j] == '(') c++;
              else if (targetRNA[j] == ')') c--;
              else if ( c <= 0 && targetRNA[j] == ')' ) 
              {
                endTarget = j;
                break;
              }
            }
	    if (endMFE != endTarget){ d++; }
   	  }
	}
	int c2 = d;
//	cout << d << endl;
//        cout << constr[0] << " -- " << d/(-1.0*(nLoops+nPoints)) << endl << endl;

        sim3 = d / (-1.0*(targetRNALoops+targetRNAPoints));

	constr[0] = sim3;
/*	if (constr[0] < sim2)
          constr[0] = sim2;
	if (constr[0] < sim3)
          constr[0] = sim3;*/
//	constr[0] = sim1 + sim2 + sim3;

/*	double constr2 = d/(-1.0*(nLoops+nPoints));
	if (constr[0] != 0)
	  constr[0] = constr2;*/
	/***********************************************************************************************************************************************************************************************************************/
	   double A, C, G, U;
	   double uA, uC, uG, uU;
	   double bGC, bAU, bGU;
	   double unpaired;
	   double paired;
	   A = C = G = U = 0;
	   unpaired = uA = uC = uG = uU = 0;
	   paired = bGC = bAU = bGU = 0;
	   
           char anterior='#';
           int repeated = 0;
           int maxRepeated = -1;
	   for (i=0; i<len; i++){
             if( anterior == seq[i]){
	        repeated++;
	     }
	     else{
                if(repeated > maxRepeated) maxRepeated = repeated; 
	        anterior = seq[i];
                repeated = 0;
	     }
	     if(seq[i] == 'A'){
	       A++;
	       if(mfe[i] == '.'){
	         uA++;
	       }
	     }
	     if(seq[i] == 'C'){
	       C++;
	       if(mfe[i] == '.'){
	         uC++;
	       }
	     }
	     if(seq[i] == 'G'){
	       G++;
	       if(mfe[i] == '.'){
	         uG++;
	       }
	     }
	     if(seq[i] == 'U'){
	       U++;
	       if(mfe[i] == '.'){
	         uU++;
	       }
	     }
	     if(mfe[i] == '.') unpaired++;
	     if(mfe[i] == '('){
	       c=1;
	       for(j=i+1; j<len; j++){
	         if(mfe[j] == '(') c++;
	         if(mfe[j] == ')') c--;
	         if(c == 0){
	           if( ((seq[i] == 'G') && (seq[j] == 'C')) || ((seq[i] == 'C') && (seq[j] == 'G')) ) bGC++;
	           if( ((seq[i] == 'A') && (seq[j] == 'U')) || ((seq[i] == 'U') && (seq[j] == 'A')) ) bAU++;
	           if( ((seq[i] == 'G') && (seq[j] == 'U')) || ((seq[i] == 'U') && (seq[j] == 'G')) ) bGU++;
	           paired++;
	           break;
	         }
	       }
	     }
	   }
	   if(paired > 0){
	     bGC = bGC / (1.0 * paired);
  	     bAU = bAU / (1.0 * paired);
	     bGU = bGU / (1.0 * paired);
	   }
	   if(unpaired > 0){
	     uA = uA / (1.0 * unpaired );
	     uC = uC / (1.0 * unpaired );
	     uG = uG / (1.0 * unpaired );
	     uU = uU / (1.0 * unpaired );
	   }
           A = A / (1.0 * len);
	   C = C / (1.0 * len);
	   G = G / (1.0 * len);
	   U = U / (1.0 * len);

	   double m1 = bGC; if(bAU > m1) m1=bAU; if(bGU > m1) m1=bGU;
	   double m2 = uA; if(uC > m2) m2=uC; if(uG > m2) m2=uG; if(uU > m2) m2=uU;
	   double m3 = A; if(C > m3) m3=C; if(G > m3) m3=G; if(U > m3) m3=U;


           /*double EuclideanDistance = sqrt(  ((0.57 - bGC) * (0.57 - bGC)) +
                        ((0.30 - bAU) * (0.30 - bAU)) +
                        ((0.13 - bGU) * (0.13 - bGU)) +
                        ((0.30 - uA) * (0.30 - uA)) +
                        ((0.20 - uC) * (0.20 - uC)) +
                        ((0.23 - uG) * (0.23 - uG)) +
                        ((0.27 - uU) * (0.27 - uU)) +
                        ((0.23 - A) * (0.23 - A)) +
                        ((0.24 - C) * (0.24 - C)) +
                        ((0.28 - G) * (0.28 - G)) +
                        ((0.24 - U) * (0.24 - U))
                  );
/*
	    m1 = sqrt(  ((0.57 - bGC) * (0.57 - bGC)) +
                        ((0.30 - bAU) * (0.30 - bAU)) +
                        ((0.13 - bGU) * (0.13 - bGU))
		);
            m2 = sqrt(  ((0.30 - uA) * (0.30 - uA)) +
                        ((0.20 - uC) * (0.20 - uC)) +
                        ((0.23 - uG) * (0.23 - uG)) +
                        ((0.27 - uU) * (0.27 - uU))
		);
	    m3 = sqrt(  ((0.23 - A) * (0.23 - A)) +
                        ((0.24 - C) * (0.24 - C)) +
                        ((0.28 - G) * (0.28 - G)) +
                        ((0.24 - U) * (0.24 - U))
		);
	    */
	    double EuclideanDistance = (m1 + m2 + m3)*100.0;
	    /* CONSTRAINT IS NOT SATISFIED */
            if(constr[0] != 0){
              obj[0] = -sim1; //c1; //constr[0]; // _mfe;/*1.0 - distance;*/
	      obj[1] = -sim2; //1; //c2; //-constr2; // EuclideanDistance; /*1.0 - distance;*/
              if(fast_mode == 0) obj[2] = -sim3; //1; //EuclideanDistance; //-constr2;//1.0 - distance;
              if(constr[0] > minSim){
                 minSim = constr[0];
		 changeBestObj = true;
        	 for(i=0; i<len; i++){
		   minSimString[i] = '_';
        	   if(mfe[i] != targetRNA[i])
		   {
		     minSimString[i] = '*';
//		     failedPos[i]++;
		   }
	         }
                 minSimString[len] = '\0';
              }
	      constr[0] = -1.0;
              double* o = new double[nobj+ncon];
              char* s = new char[strlen(sequence)+1];
              strcpy(s, sequence);
              o[0] = obj[0];
              o[1] = obj[1];
	      if (fast_mode == 1)
	      {
		o[2] = constr[0];
	      }
	      else
	      {
	        o[2] = obj[2];
                o[3] = constr[0];
	      }
	      if (oracle.size() > 100000)// prevent RAM in experiments eterna COMMENT
              {
     	        for (itr = oracle.begin(); itr != oracle.end(); itr++)
    		{
      		  char* s = itr->first;
      		  double* o = itr->second;
      		  free(s);
      		  free(o);
    		}
		oracle.clear();
	      }
              oracle.insert( pair<char*, double* > (s, o) );
	      vrna_fold_compound_free(vc);
	      nEvMFE++;
              return;
            }

	    nConstr++;
	    validSolutionsFound++;

	    /* OTHERWISE */
            if(fast_mode == 1){
	      nEvMFE++;
              obj[0] = _mfe;
              obj[1] = EuclideanDistance;
              if(obj[0] < bestObj[0]){ changeBestObj = true; bestObj[0] = obj[0]; }
              if(obj[1] < bestObj[1]){ changeBestObj = true; bestObj[1] = obj[1]; }
              double* o = new double[nobj+ncon];
              char* s = new char[strlen(sequence)+1];
              strcpy(s, sequence);
              o[0] = obj[0];
              o[1] = obj[1];
              o[2] = constr[0];
	      if (oracle.size() > 100000)// prevent RAM in experiments eterna COMMENT
              {
     	        for (itr = oracle.begin(); itr != oracle.end(); itr++)
    		{
      		  char* s = itr->first;
      		  double* o = itr->second;
      		  free(s);
      		  free(o);
    		}
		oracle.clear();
	      }
              oracle.insert( pair<char*, double* >(s, o) );
	    }
	    else
            {
	      nEvPFBP++;
              double kT = (temperature+273.15)*1.98717/1000.;
              pf_scale = exp(-_mfe/kT/strlen(sequence));
              double _pfe = (double) pf_fold (sequence, pfold_structure );
              FLT_OR_DBL *pr = export_bppm();
              double mbp = vrna_mean_bp_distance_pr (strlen(sequence), pr );
              obj[0] = _pfe;
              obj[1] = mbp; /* because obj0 y obj1 are for minimization */
              obj[2] = EuclideanDistance;
              if(obj[0] < bestObj[0]){ changeBestObj = true; bestObj[0] = obj[0]; }
              if(obj[1] < bestObj[1]){ changeBestObj = true; bestObj[1] = obj[1]; }
              if(obj[2] < bestObj[2]){ changeBestObj = true; bestObj[2] = obj[2]; }
              double* o = new double[nobj+1];
	      //string s(seq);
              char* s = new char[len+1];
              strcpy(s, sequence);
              o[0] = obj[0];
              o[1] = obj[1];
              o[2] = obj[2];
              o[3] = constr[0];
	      if (oracle.size() > 100000)// prevent RAM in experiments eterna COMMENT
              {
     	        for (itr = oracle.begin(); itr != oracle.end(); itr++)
    		{
      		  char* s = itr->first;
      		  double* o = itr->second;
      		  free(s);
      		  free(o);
    		}
		oracle.clear();
	      }
              oracle.insert( pair<char*, double* >(s, o) );
	   }
           if(constr[0] > minSim){
              minSim = constr[0];
           }
	   nMinSim++;
           vrna_fold_compound_free(vc);
   }
   else
   {
     nEvOracle++;
//     cout << "No evaluada" << endl;
     if (fast_mode == 1)
     {
       obj[0] = itr->second[0];
       obj[1] = itr->second[1];
       constr[0] = itr->second[2];
       if (constr[0] == 0)
       {
	  strcpy(mfe, targetRNA);
          if(obj[0] < bestObj[0]){ changeBestObj = true; bestObj[0] = obj[0]; }
          if(obj[1] < bestObj[1]){ changeBestObj = true; bestObj[1] = obj[1]; }
       }
       else
          strcpy(mfe, "N/A");
     }
     else
     {
       obj[0] = itr->second[0];
       obj[1] = itr->second[1];
       obj[2] = itr->second[2];
       constr[0] = itr->second[3];
       if (constr[0] == 0)
       {
	  strcpy(mfe, targetRNA);
          if(obj[0] < bestObj[0]){ changeBestObj = true; bestObj[0] = obj[0]; }
          if(obj[1] < bestObj[1]){ changeBestObj = true; bestObj[1] = obj[1]; }
          if(obj[2] < bestObj[2]){ changeBestObj = true; bestObj[2] = obj[2]; }
       }
       else
          strcpy(mfe, "N/A");
     }
     if(constr[0] > minSim)
       minSim = constr[0];
   }
   return;
/*-------------------------------------------------------------------------------------------------------------------------*/
}

/*********************************************************************************************************************************************************************************/
/* RANK */
/*********************************************************************************************************************************************************************************/
/* Function to assign rank and crowding distance to a population of size pop_size*/
void Nsga2::assign_rank_and_crowding_distance (population *new_pop)
{
    int flag;
    int i;
    int end;
    int front_size;
    int rank=1;
    list_ *orig;
    list_ *cur;
    list_ *temp1, *temp2;
    orig = (list_ *)malloc(sizeof(list_));
    cur = (list_ *)malloc(sizeof(list_));
    front_size = 0;
    orig->index = -1;
    orig->parent = NULL;
    orig->child = NULL;
    cur->index = -1;
    cur->parent = NULL;
    cur->child = NULL;
    temp1 = orig;
    for (i=0; i<popsize; i++)
    {
        insert (temp1,i);
        temp1 = temp1->child;
    }
    do
    {
        if (orig->child->child == NULL)
        {
            new_pop->ind[orig->child->index].rank = rank;
            new_pop->ind[orig->child->index].crowd_dist = INF;
            break;
        }
        temp1 = orig->child;
        insert (cur, temp1->index);
        front_size = 1;
        temp2 = cur->child;
        temp1 = del (temp1);
        temp1 = temp1->child;
        do
        {
            temp2 = cur->child;
            do
            {
                end = 0;
                flag = check_dominance (&(new_pop->ind[temp1->index]), &(new_pop->ind[temp2->index]));
                if (flag == 1)
                {
                    insert (orig, temp2->index);
                    temp2 = del (temp2);
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0)
                {
                    temp2 = temp2->child;
                }
                if (flag == -1)
                {
                    end = 1;
                }
            }
            while (end!=1 && temp2!=NULL);
            if (flag == 0 || flag == 1)
            {
                insert (cur, temp1->index);
                front_size++;
                temp1 = del (temp1);
            }
            temp1 = temp1->child;
        }
        while (temp1 != NULL);
        temp2 = cur->child;
        do
        {
            new_pop->ind[temp2->index].rank = rank;
            temp2 = temp2->child;
        }
        while (temp2 != NULL);
        assign_crowding_distance_list (new_pop, cur->child, front_size);
        temp2 = cur->child;
        do
        {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
        while (cur->child !=NULL);
        rank+=1;
    }
    while (orig->child!=NULL);
    free (orig);
    free (cur);
    return;
}

/*********************************************************************************************************************************************************************************/
/* REPORT */
/*********************************************************************************************************************************************************************************/
void Nsga2::getResults(population *pop)
{
  int i,j;
  for (i=0; i<popsize; i++)
  {
    if (pop->ind[i].constr[0] == 0)// && pop->ind[i].rank == 1)
    {
      string s(pop->ind[i].sequence);
      cout << " ** " << s << "**" << endl;
      results.push_back(s);
    }
  }
}

void Nsga2::getOracle()
{
  for (itr = oracle.begin(); itr != oracle.end(); itr++)
  {
    if (itr->second[3] == 0) // constr[0] == 0
      results.push_back(itr->first);
  }
}

void Nsga2::report_pop_xreal(population* pop)
{
  for (int i=0; i<popsize; i++)
  {
    cout << "IND" << i  << ":\t";
    for (int j=0; j<nreal; j++)
      cout << pop->ind[i].xreal[j] << "\t";
    cout << endl;
  }
}

void Nsga2::report_ind (individual* ind)
{
    cout << endl;
        if(ind->constr[0] == 0){
        if(stop_crit < 0)
	  printf("#");
        printf("%f\t",ind->obj[0]);
        printf("%f\t",ind->obj[1]);
        if(fast_mode == 0) 
          printf("%f\t",ind->obj[2]);
        }
        else{
          printf("#[No Successful RNA Sequence]\t");
        }
        printf("%f\t",ind->constr[0]);
        printf("%s\t%s\n", ind->sequence, ind->mfe);
}

/* Function to print the information of a population in a file */
void Nsga2::report_pop_ (population *pop, int flag)
{
    int i,j;
    double cons = 0;
    for (i=0; i<popsize; i++)
    {
      if(flag == 1){
        if(pop->ind[i].constr[0] == 0){
        printf("#%f\t",pop->ind[i].obj[0]);
        printf("%f\t",pop->ind[i].obj[1]);
        if(fast_mode == 0) 
          printf("%f\t",pop->ind[i].obj[2]);
	}
	else{
	  printf("#[No Successful RNA Sequence]\t");
	}
	printf("%f\t",pop->ind[i].constr[0]);
	printf("%s\t%s\n", pop->ind[i].sequence, pop->ind[i].mfe);
     }
     else{
        if(pop->ind[i].rank==1){
	        printf("#%f\t",pop->ind[i].obj[0]);
	        printf("%f\t",pop->ind[i].obj[1]);
	        if(fast_mode == 0)
        	  printf("%f\t",pop->ind[i].obj[2]);
	        printf("%f\t",pop->ind[i].constr[0]);
	        printf("%s\t%s\n", pop->ind[i].sequence, pop->ind[i].mfe);
        }
     }
    }
}

void Nsga2::report_feasible_ (population *pop)
{
	int i, j, k;
	printf("\n \t");
    for (i=0; i<popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
        {
            /*for (j=0; j<nobj; j++)
            {*/
                /*printf("%.3f\t",pop->ind[i].obj[0]/(-pop->ind[i].obj[1]));*/
                printf("%.3f\t",pop->ind[i].obj[0]);
                printf("%.3f\t",pop->ind[i].obj[1]);
                /*printf("%.3f\t",pop->ind[i].obj[2]);*/
				/*printf("%.3f\t",pop->ind[i].constr[0]);*/
            /*}*/
			printf("\n \t");
		}
	}
}

void Nsga2::report_pop2 (population *pop)
{
	int i, j, k;
	printf("\n \t");
    for (i=0; i<popsize; i++)
    {
            printf("%.3f\t",pop->ind[i].obj[0]);
            printf("%.3f\t",pop->ind[i].obj[1]);
            /*printf("%.3f\t",pop->ind[i].obj[2]);*/
			printf("%.3f\t",pop->ind[i].constr[0]);
			printf("\n \t");
	}
}

/* Function to print the information of feasible and non-dominated population in a file */
void Nsga2::report_feasible (population *pop, FILE *fpt)
{
    int i, j, k;
	double value;
    for (i=0; i<popsize; i++)
    {
        /*if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
        {*/
            for (j=0; j<nobj; j++)
            {
                fprintf(fpt,"%.3f\t",pop->ind[i].obj[j]);
            }
            if (ncon!=0)
            {
                for (j=0; j<ncon; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
                }
            }
            if (nreal!=0)
            {
		    /*char seq[nreal+1];
			strcpy(seq,"");
		    for(j=0; j<nreal; j++){
			 if( (pop->ind[i].xreal[j] >= 0.00) && (pop->ind[i].xreal[j] <  0.25) ) seq[j]='G';
			 if( (pop->ind[i].xreal[j] >= 0.25) && (pop->ind[i].xreal[j] <  0.50) ) seq[j]='C';
			 if( (pop->ind[i].xreal[j] >= 0.50) && (pop->ind[i].xreal[j] <  0.75) ) seq[j]='U';
			 if( (pop->ind[i].xreal[j] >= 0.75) && (pop->ind[i].xreal[j] <= 1.00) ) seq[j]='A';	 
		    }
			seq[nreal]='\0';*/
			  int len = strlen(targetRNA);
			  /*char *seq = (char*) malloc (sizeof(char) * (len));*/
			  for (j=0; j<nLoops; j++){
				value = pop->ind[i].xreal[j];
				/* CG (0 to 0.16), GC (0.16 to 0.33), AU (0.33 to 0.49), UA (0.49 to 0.66), GU (0.66 to 0.83), or UG (0.83 to 1.00) */
				/*if ( (value >= 0) && (value < 0.16) )*/
				if ( (value >= 0) && (value < 0.16) )
				{
					seq[initLoopTargetRNA[j]] = 'C';
					seq[endLoopTargetRNA[j]]  = 'G';	
				}
				/*if ( (value >= 0.16) && (value < 0.33) )*/
				if ( (value >= 0.16) && (value < 0.33) )
				{
					seq[initLoopTargetRNA[j]] = 'G';
					seq[endLoopTargetRNA[j]]  = 'C';	
				}	
				/*if ( (value >= 0.33) && (value < 0.49) )*/
				if ( (value >= 0.33) && (value < 0.49) )
				{
					seq[initLoopTargetRNA[j]] = 'A';
					seq[endLoopTargetRNA[j]]  = 'U';	
				}	
				/*if ( (value >= 0.49) && (value < 0.66) )*/
				if ( (value >= 0.49) && (value < 0.66) )
				{
					seq[initLoopTargetRNA[j]] = 'U';
					seq[endLoopTargetRNA[j]]  = 'A';	
				}		
				/*if ( (value >= 0.7) && (value < 0.83) )*/
				if ( (value >= 0.66) && (value < 0.83) )
				{
					seq[initLoopTargetRNA[j]] = 'G';
					seq[endLoopTargetRNA[j]]  = 'U';	
				}		
				/*if ( (value >= 0.83) && (value <= 1.00) )*/
				if ( (value >= 0.83) && (value <= 1.00) )
				{
					seq[initLoopTargetRNA[j]] = 'U';
					seq[endLoopTargetRNA[j]]  = 'G';	
				}		
			  }
			  for(j=0; j<nPoints; j++){
				value = pop->ind[i].xreal[nLoops+j];
				if ( (value >= 0) && (value < 0.25) )
				{
					seq[pointTargetRNA[j]] = 'C';
				}	
				if ( (value >= 0.25) && (value < 0.5) )
				{
					seq[pointTargetRNA[j]] = 'G';
				}		
				if ( (value >= 0.5) && (value < 0.75) )
				{
					seq[pointTargetRNA[j]] = 'A';
				}		
				if ( (value >= 0.75) && (value <= 1.00) )
				{
					seq[pointTargetRNA[j]] = 'U';
				}
			  }
			  /*seq[len]='\0';*/
			
			fprintf(fpt, "%s\t%s", seq);
            }
            if (nbin!=0)
            {
                for (j=0; j<nbin; j++)
                {
                    for (k=0; k<nbits[j]; k++)
                    {
                        fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                    }
                }
            }
            fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
            fprintf(fpt,"%d\t",pop->ind[i].rank);
            fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
        /*}*/
    }
    return;
}

/*********************************************************************************************************************************************************************************/
/* TOURSELECT */
/*********************************************************************************************************************************************************************************/
/* Routine for tournament selection, it creates a new_pop from old_pop by performing tournament selection and the crossover */
void Nsga2::selection (population *old_pop, population *new_pop)
{
    int *a1, *a2;
    int temp;
    int i;
    int rand;
    individual *parent1, *parent2;
    a1 = (int *)malloc(popsize*sizeof(int));
    a2 = (int *)malloc(popsize*sizeof(int));
    for (i=0; i<popsize; i++)
    {
        a1[i] = a2[i] = i;
    }
    for (i=0; i<popsize; i++)
    {
        rand = rndGenerator->rnd (i, popsize-1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = rndGenerator->rnd (i, popsize-1);
        temp = a2[rand];
        a2[rand] = a2[i];
        a2[i] = temp;
    }

    for (i=0; i<popsize; i+=4)
    {
	if (rndGenerator->rndreal(0.0,1.0) <= 0.7)
	{
            parent1 = tournament (&old_pop->ind[a1[i]], &old_pop->ind[a1[i+1]]);
            parent2 = tournament (&old_pop->ind[a1[i+2]], &old_pop->ind[a1[i+3]]);
            crossover (parent1, parent2, &new_pop->ind[i], &new_pop->ind[i+1]);
            parent1 = tournament (&old_pop->ind[a2[i]], &old_pop->ind[a2[i+1]]);
            parent2 = tournament (&old_pop->ind[a2[i+2]], &old_pop->ind[a2[i+3]]);
            crossover (parent1, parent2, &new_pop->ind[i+2], &new_pop->ind[i+3]);
	}
	else
	{
          initialize_ind0(&new_pop->ind[i+0]);
          initialize_ind1(&new_pop->ind[i+1]);
          initialize_ind2(&new_pop->ind[i+2]);
          initialize_ind3(&new_pop->ind[i+3]);
	}
    }
    free (a1);
    free (a2);
    return;
}

/* Routine for binary tournament */
individual* Nsga2::tournament (individual *ind1, individual *ind2)
{
    int flag;
    flag = check_dominance (ind1, ind2);
    if (flag==1)
    {
        return (ind1);
    }
    if (flag==-1)
    {
        return (ind2);
    }
    if (ind1->crowd_dist > ind2->crowd_dist)
    {
        return(ind1);
    }
    if (ind2->crowd_dist > ind1->crowd_dist)
    {
        return(ind2);
    }
    if ((rndGenerator->randomperc()) <= 0.5)
    {
        return(ind1);
    }
    else
    {
        return(ind2);
    }
}

/*********************************************************************************************************************************************************************************/
/* OLD DECOMPOSITION (TO DELETE) */
/*********************************************************************************************************************************************************************************/
/*
void Nsga2::Pairs_(char* input, Substructure* ss)
{
  int c;
  int l = strlen(input);
  Pair_T* pair = new Pair_T[l+1];
  int n = 0;
  for (int i = 0; i < l; i++) {
    if (input[i] == '(') {
      c = 0;
      for (int j = i; j < l; j++) {
        if (input[j] == '(')
          c++;
        if (input[j] == ')')
          c--;
        if ((c <= 0) && (input[j] == ')')) {
          pair[n].i = i;
          pair[n].j = j;
          n++;
          break;
        }
      }
    }
  }

  int size = 0;
  int* includedPair = new int[n];
  ss->global = true;
  for (int p=0; p<n; p++)
  {
    includedPair[p] = -1;
    int i = pair[p].i;
    int j = pair[p].j;
    bool condition1, condition2;
    bool cBracket = false;
    condition1 = true; // CONDITION 1: NOT contains substructures
    for (int q=i+1; (q < j && condition1); q++)
    {
      if (input[q] == ')') cBracket=true;
      if (input[q] == '(' && cBracket) condition1=false;
    }
    condition2 = true; // CONDITION 2: is not OVERLAPPED
    for (int q=0; (q < ss->n && condition2); q++)
    {
      for (int r=0; (r < ss->substr[q].n && condition2); r++)
      {
        if (i > ss->substr[q].pair[r].i && j < ss->substr[q].pair[r].j) // isOverlapped
          if (ss->substr[q].containsSubstr && condition1)
            condition2 = true;
          else
            condition2 = false;
      }
    }
    if (condition2)
    {
        if (j-i+1 == l) ss->global=false;
//      {
        char* str = new char[j-i+2];
        for (int s=i; s<=j; s++) str[s-i] = input[s];
        str[j-i+1] = '\0';
        bool found = false;
        for (int s=0; (s < ss->n && !found); s++)
        {
          if (strcmp(ss->substr[s].structure, str) == 0)
	  {
	    ss->substr[s].pair[ss->substr[s].n].i = i;
	    ss->substr[s].pair[ss->substr[s].n].j = j;
	    ss->substr[s].n++;
	    found = true;
	    includedPair[p] = s;
	  }
        }
        if (!found)
        {
          includedPair[p] = ss->n;
	  ss->substr[ss->n].structure = new char[j-i+2];
	  ss->substr[ss->n].structure_uncompressed = new char[j-i+2];
	  strcpy(ss->substr[ss->n].structure, str);
	  strcpy(ss->substr[ss->n].structure_uncompressed, str);
	  ss->substr[ss->n].pair = new Pair_T[n];
	  ss->substr[ss->n].pair[0].i = i;
	  ss->substr[ss->n].pair[0].j = j;
	  ss->substr[ss->n].n = 1;
          if (condition1 == true)
          {
            ss->substr[ss->n].containsSubstr = false;
          }
          else
          {
            ss->substr[ss->n].containsSubstr = true;
            ss->n_containsSubstr++;
          }
          ss->n++;
        }
        free(str);
//      }
    }
  }
  if (ss->global)
  {
    ss->substr[ss->n].structure = new char[l+1];
    ss->substr[ss->n].structure_uncompressed = new char[l+1];
    strcpy(ss->substr[ss->n].structure, input);
    strcpy(ss->substr[ss->n].structure_uncompressed, input);
    ss->substr[ss->n].pair = new Pair_T[n];
    ss->substr[ss->n].pair[0].i = -1;
    ss->substr[ss->n].pair[0].j = -1;
    ss->substr[ss->n].n = 1;
    ss->substr[ss->n].containsSubstr = true;
    ss->n_containsSubstr++;
    ss->n++;
  }

  for (int i=0; i < ss->n; i++)
  {
    if (ss->substr[i].containsSubstr)
    {
      int lSubstr = strlen(ss->substr[i].structure_uncompressed);
      strcpy(ss->substr[i].structure, "");
      int c1=0;
      for (int j=0; j < lSubstr; j++)
      {
        if (ss->substr[i].structure_uncompressed[j] == '.')
	{
          ss->substr[i].structure[c1] = '.';
	  c1++;
	}
	else if (ss->substr[i].structure_uncompressed[j] == ')')
	{
          ss->substr[i].structure[c1] = ')';
          c1++;
	}
	else if (ss->substr[i].structure_uncompressed[j] == '(')
        {
          int c = 0;
	  bool found = false;
	  int contSub;
          for (int k=j; k < lSubstr; k++)
	  {
	    if (ss->substr[i].structure_uncompressed[k] == '(')        c++;
            else if (ss->substr[i].structure_uncompressed[k] == ')')   c--;
            if ((c <= 0) && (ss->substr[i].structure_uncompressed[k] == ')'))
	    {
	      int init = j;
              int end = k;
	      char* str = new char[k-j+2];
              for (int s=j; s<=k; s++) str[s-j] = ss->substr[i].structure_uncompressed[s];
              str[k-j+1] = '\0';
	      int s = -1;
              for (int s1 = 0; (s1 < ss->n && s == -1); s1++)
	      {
	        if (strcmp(str,ss->substr[s1].structure_uncompressed) == 0 && s1 != i) s = s1;
	      }
              if (s != -1)
              {
		ss->substr[i].substructureID.push_back(s);
		ss->substr[i].structure[c1] = '$';
		c1++;
		found = true;
		contSub = k;
	      }
	      free(str);
	      break;
	    }
	  }
          if (!found)
          {
            ss->substr[i].structure[c1] = '(';
            c1++;
          }
          else
          {
            j=contSub;
          }
        }
      }
      ss->substr[i].structure[c1] = '\0';
    }
  }
  free(pair);
  free(includedPair);
}

Substructure* Nsga2::Pairs(char* input)
{
   // Wrap method for recursion
   Substructure* ss = new Substructure;
   ss->substr = new SubStructure_T[strlen(input)+1];
   ss->n = 0;
   ss->n_containsSubstr = 0;
   Pairs_(input, ss);
   return ss;
}

char* Nsga2::getGlobalSequence(Substructure* ss, int lengthT, int str)
{
        int c2 = strlen(ss->substr[str].structure);
        mask = new char[lengthT+2]; // GLOBAL VARIABLE MASK
	int rm = 0;
	int dollar=0;
        char* resultant = new char[lengthT+2];
        int rc=0;
        strcpy(resultant,"");
        strcpy(mask, "");
        for (int i=0; i < c2; i++)
        {
          if (ss->substr[str].structure[i] == '$')
          {
	    int subseq = ss->substr[str].substructureID.at(dollar);
            if (ss->substr[subseq].sequences.size() == 0)
            {
              int l = strlen(ss->substr[subseq].structure);
	      int dollar_aux = 0;
              for (int k = 0; k < l; k++)
              {
	        if (ss->substr[subseq].structure[k] == '$')
	        {
	          int subseq_aux = ss->substr[subseq].substructureID.at(dollar_aux);
                  int l_aux = strlen(ss->substr[subseq_aux].structure);
                  for (int k = 0; k < l_aux; k++)
                  {
                    mask[rm] = '*';;
                    rm++;
                  }
		  dollar_aux++;
		}
		else
		{
		  mask[rm] = ss->substr[subseq].structure[k];
		  rm++;
		}
                resultant[rc] = ss->substr[subseq].structure[k];
                rc++;
              }
	      ss->substr[str].substructureID.erase(ss->substr[str].substructureID.begin() + dollar);
	      for (int k=0; k < ss->substr[subseq].substructureID.size(); k++)
		ss->substr[str].substructureID.insert(ss->substr[str].substructureID.begin() + dollar + k, ss->substr[subseq].substructureID.at(k));
            }
            else
            {
	        dollar++;
		resultant[rc] = '$';
		rc++;
                int l = strlen(ss->substr[subseq].structure_uncompressed);
                for (int k = 0; k < l; k++)
                {
                  mask[rm] = '*';;
                  rm++;
                }
            }
          }
          else
          {
		mask[rm] = ss->substr[str].structure[i];
		rm++;
                resultant[rc] = ss->substr[str].structure[i];
                rc++;
          }
        }
        resultant[rc] = '\0';
	mask[rm] = '\0';
	return resultant;
}

void Nsga2::RNAStudy(Substructure* ss, char* target)
{
        int total = strlen(target);
        int stotal = 0;
        char s[1000];
        strcpy(s,"");
        cout << target << endl;
        for (int i=0; i < ss->n; i++)
        {
          stotal += (ss->substr[i].n * strlen(ss->substr[i].structure));
          if (ss->substr[i].n > 1)
            sprintf(s,"%s %d(%d x%d)",s,(ss->substr[i].n * strlen(ss->substr[i].structure)), strlen(ss->substr[i].structure),ss->substr[i].n);
          else
            sprintf(s,"%s %d",s,(ss->substr[i].n * strlen(ss->substr[i].structure)));
          cout << "SS#" << i <<":\t" << ss->substr[i].structure << "\t{ ";
	  for (int j=0; j < ss->substr[i].substructureID.size(); j++) cout << ss->substr[i].substructureID.at(j) << " ";
	  cout << "}\tPairs: ";
          for (int j=0; j < ss->substr[i].n; j++)
          {
            cout << "(" << ss->substr[i].pair[j].i << ", " << ss->substr[i].pair[j].j << ") ";
          }
	  cout << endl;
	  if (ss->substr[i].sequences.size() > 0) cout << " \t#SEQUENCES_FOUND: " << ss->substr[i].sequences.size() << endl;
        }
        cout << "Target length: " << total << "\t#Subsequences: " << ss->n - 1 << " (" << s << " )\tResto: " <<total -stotal << endl;
        cout << endl << endl;
}

*/
