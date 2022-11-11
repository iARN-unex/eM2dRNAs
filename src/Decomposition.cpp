#include "Decomposition.h"

using namespace std;

Decomposition::Decomposition(string targetRNA_, int randomSeed_, int popSize_, double stopCrit_, bool fastMode_, int model_, bool solvingFlag_, bool decompositionFlag_)
{
  targetRNA = targetRNA_;
  randomSeed = randomSeed_;
  popSize = popSize_;
  stopCrit = stopCrit_;
  fastMode = fastMode_;
  solvingFlag = solvingFlag_;
  decompositionFlag = decompositionFlag_;
  model = model_;
  if (model == MODEL_TURNER1999) vrna_params_load_RNA_Turner1999();

  iPrint = 0;
  stepTime = -1;
  if (stopCrit < 0)
  {
    stepTime = -stopCrit/50.0; // We study the solutions found by the algorithm every Runtime/50 seconds (number of steps)
  }
}

void Decomposition::updateParticipation()
{
  totalParticipation = 0;
  for (int i=0; i<problems.size(); i++)
    if (problems.at(i)->isAttempted() == false)
      totalParticipation += problems.at(i)->getTargetRNA().size();

  for (int i=0; i<problems.size(); i++)
//    if (problems.at(i)->getID() == 0)
//      problems.at(i)->setParticipation((problems.at(i)->getTargetRNA().size()*2.0)/totalParticipation);
//    else
    if (problems.at(i)->isAttempted() == false)
      problems.at(i)->setParticipation((problems.at(i)->getTargetRNA().size()*1.0)/totalParticipation);

  double sum = 0;
  for (int i=0; i<problems.size(); i++)
    sum += problems.at(i)->getParticipation();

//  cout << "Total Participation = " << sum << endl;
}

void Decomposition::decompose()
{
  if (!decompositionFlag)
  {
    cout << "*****************************************************************************************************************************************************" << endl;
    cout << "*                                                                                                                                                   *" << endl;
    cout << "* SOLVING THE PROBLEM WITHOUT DECOMPOSING THE TARGET STRUCTURE                                                                                      *" << endl;
    cout << "*                                                                                                                                                   *" << endl;
    cout << "*****************************************************************************************************************************************************" << endl;
    // Solve the problem without decompose at first, 0.1 of time/#evals
    //global = new Problem;
    //global->setTargetRNA(targetRNA, 1);
    //global->setID(0);
//    global->solve(randomSeed, popSize, stopCrit, fastMode, solvingFlag);
/*    cout << "TARGET:       " << targetRNA << endl;
    cout << "SIZE:         " << targetRNA.size() << endl;
    if (stopCrit > 0)
      cout << "STOPCRIT:     " <<  stopCrit << " local evaluations (" << stopCrit << " global evaluations)" << endl;
    else
      cout << "STOPCRIT:     " << -stopCrit << " seconds" << endl;
    cout << "*********************************************************************************************************************************************************" << endl;
    cout << endl;// << endl << endl << endl;
    global->solve(randomSeed, popSize, stopCrit, false, false);
*/
    Problem* p = new Problem;
    p->setTargetRNA(targetRNA, 1);
    this->addProblem(p);
    totalParticipation += targetRNA.size();
    updateParticipation();
  }
  else
  {
    if (solvingFlag)
    {
      cout << "*****************************************************************************************************************************************************" << endl;
      cout << "*                                                                                                                                                   *" << endl;
      cout << "* STEP 0: TRY TO SOLVE WITHOUT DECOMPOSING (SOLVING FLAG ACTIVATED)                                                                                 *" << endl;
      cout << "*                                                                                                                                                   *" << endl;
      cout << "*****************************************************************************************************************************************************" << endl;
      Problem* p = new Problem;
      p->setTargetRNA(targetRNA, 1);
      p->setID(0);
      cout << endl;
      p->solve(randomSeed, popSize, stopCrit*0.01, solvingFlag);
      cout << stopCrit << endl;
      stopCrit -= p->getStopCrit();
      cout << stopCrit << endl;
    }
    cout << "*****************************************************************************************************************************************************" << endl;
    cout << "*                                                                                                                                                   *" << endl;
    cout << "* STEP 1: DECOMPOSE THE TARGET STRUCTURE INTO SMALL SUB-STRUCTURES                                                                                  *" << endl;
    cout << "*                                                                                                                                                   *" << endl;
    cout << "*****************************************************************************************************************************************************" << endl;
    problems.clear();
    totalParticipation = 0;
    recursive(targetRNA);
    updateParticipation();

    int numberOfSubproblems = problems.size();
    bool* visited = new bool[numberOfSubproblems]; // Some subproblems may be deleted in a previous step
    for (int i=0; i < numberOfSubproblems; i++) visited[i] = false;

    vector<int> id;
    for (int i=0; i < problems.size(); i++)
      id.push_back(problems.at(i)->getID());

    // Delete all problems with length < 20 (insignificant) 
    for (int i=0; i < id.size(); i++)
    {
      Problem* p = getProblemByID(id.at(i));//problems.at(i);
      visited[p->getID()] = false;

      //if ( p->getTargetRNA().size() < 20 && p->getID() !=0 )
      if ( (p->getTargetRNA().size() * p->getParentsSize() ) < 20 && p->getID() != 0)
      {
	cout << "DELETING PROBLEM WITH ID: " << p->getID() << " | LENGTH < 20"<< endl;
        visited[p->getID()] = true;
        // Delete problem and re-organize structure
        p->remove();
        deleteProblemByID(p->getID(), true);
        free(p);
      }
    }

    id.clear();
    for (int i=0; i < problems.size(); i++)
      id.push_back(problems.at(i)->getID());

    // Delete all problems with ONLY ONE CHILD (except ID=0, global)
    for (int i=0; i < id.size(); i++)
    {
      Problem* p = getProblemByID(id.at(i)); //problems.at(i);
      if (p->getChildrenSize() == 1 && p->getID() != 0)
      {
	cout << "DELETING PROBLEM WITH ID: " << p->getID() << " | ONLY CHILD"<< endl;
        visited[p->getID()] = true;
        // Delete problem and re-organize structure
        p->remove();
        deleteProblemByID(p->getID(), true);
        free(p);
      }
    }

    // Finally, if problem with ID=0 only contain one children, dismiss it!
    Problem* p = getProblemByID(0);
    if (p->getChildrenSize() == 1)
    {
      Problem* child = p->getChildrenPairAt(0).second;
      cout << "DELETING PROBLEM WITH ID: " << child->getID() << " | IT IS ONLY CHILD IN GLOBAL"<< endl;
      visited[child->getID()] = true;
      child->remove();
      deleteProblemByID(child->getID(), true);
      free(child);
    }

    for(int i=0; i<numberOfSubproblems; i++)
      cout << "Subproblem ID " << i <<": " << visited[i] << endl;

    updateParticipation();
    topologicalSort(visited, numberOfSubproblems);
  }
}

Problem* Decomposition::recursive(string s)
{
  Problem* p = new Problem;
  p->setTargetRNA(s, ((1.0*s.size()) / targetRNA.size()));
  this->addProblem(p);
//  if (p->getID() == 0)
//    totalParticipation += 2*s.size(); // Double participation for GLOBAL
//  else
    totalParticipation += s.size();

  int i_prev = -1;
  int j_prev = s.size();

  for (int i=0; i<s.size(); i++)
  {
    if (s[i] == '(')
    {
      string s_;
      int c=0;
      int j;
      bool condition = true;
      bool cBracket = false; // CONDITION: The structure DOES NOT contain substructures
      for (j=i; j<s.size(); j++)
      {
        if (s[j] == ')'){ c--; cBracket=true; }
        if (s[j] == '('){ c++; if(cBracket) condition=false; }
        if (c == 0)
        {
          s_ = s.substr(i,j-i+1);
	  if ( s.compare(s_) != 0 && (i-1!=i_prev || j+1!=j_prev) )
	  {
	    Problem *p1 = this->checkIfExists(s_);
            if (p1 == NULL )
            {
	      if (!condition)
	      {
	        p1 = recursive(s_);
	      }
	      else
	      {
	        p1 = new Problem;
  		totalParticipation += s_.size();
	        p1->setTargetRNA(s_, ((1.0*s_.size()) / targetRNA.size()));
	        this->addProblem(p1);
	      }
            }
	    p1->addParent(i, p);
	    p->addChild(i, p1);
            i=j;
          }
	  i_prev = i;
          j_prev = j;
          break;
        }
      }
    }
  }
  return p;
}

Problem* Decomposition::checkIfExists(string s)
{
  for (int i=0; i < problems.size(); i++)
  {
    if (problems.at(i)->getTargetRNA().compare(s) == 0) return problems.at(i);
  }
  return NULL;
}

void Decomposition::topologicalSortUtil(int v, bool visited[], stack<Problem*>& Stack)
{
    Problem* p = getProblemByID(v);
    //cout << "-> P" << v << " (" << p->getParentsSize() << ")" << endl;
    visited[v] = true;

    for (int i=0; i < p->getParentsSize(); i++)
    {
      //cout << "    " << p->getParentPairAt(i).second->getID() << endl;
      if (!visited[p->getParentPairAt(i).second->getID()])
        topologicalSortUtil(p->getParentPairAt(i).second->getID(), visited, Stack);
    }

    Stack.push(p);
}

void Decomposition::topologicalSort(bool* visited, int numberOfSubproblems)
{
    stack<Problem*> Stack;

    for (int i = 0; i < numberOfSubproblems; i++)
        if (visited[i] == false)
            topologicalSortUtil(i, visited, Stack);

    problems.clear();
    while (Stack.empty() == false) {
        problems.push_back(Stack.top());
        cout << Stack.top()->getID() << " ";
        Stack.pop();
    }
}

vector<Problem*> Decomposition::getProblems()
{
  return problems;
}

int Decomposition::getProblemSize()
{
  return problems.size();
}

Problem* Decomposition::getProblemAt(int i)
{
  return problems.at(i);
}

Problem* Decomposition::getProblemByID(int id)
{
  for (int i=0; i<problems.size(); i++)
    if (problems.at(i)->getID() == id) return problems.at(i);
}

void Decomposition::addProblem(Problem* p)
{
  p->setID(problems.size());
  problems.push_back(p);
}

void Decomposition::deleteProblem(int i, bool update)
{
  if (update) totalParticipation -= problems.at(i)->getTargetRNA().size();
  problems.erase(problems.begin()+i);
  if (update) updateParticipation();

  for (int i=0; i < problems.size(); i++)
    problems.at(i)->sortChildrenParents();
}

void Decomposition::deleteProblemByID(int id, bool update)
{
  int i;
  for (i=0; i<problems.size(); i++)
    if (problems.at(i)->getID() == id) break;

  if (update) totalParticipation -= problems.at(i)->getTargetRNA().size();
  problems.erase(problems.begin()+i);
  if (update) updateParticipation();

  for (int i=0; i < problems.size(); i++)
    problems.at(i)->sortChildrenParents();
}

void Decomposition::print()
{
  for (int i=0; i < problems.size(); i++)
  {
    cout << " " << string(100,'-') << endl;
    problems.at(i)->print();
    if (stopCrit > 0)
      cout << "   A total of " << stopCrit * problems.at(i)->getParticipation() << ", that is to say, " << (stopCrit * problems.at(i)->getParticipation() ) * (pow(targetRNA.size(), index) / pow(problems.at(i)->getTargetRNA().size(), index) ) << " local evaluations" << endl;
    else
      cout << "   A total of " << stopCrit * problems.at(i)->getParticipation() << " seconds to evolve" << endl;
    cout << endl;
  }
}

void Decomposition::printGraphviz()
{
//  Problem* p = getProblemByID(3);
//  p->remove();
//  deleteProblemByID(3, true);
//  free(p);

/*  vector<int> problemID;
  for (int n=0; n < problems.size(); n++)
    problemID.push_back(problems.at(n)->getID());

  for (int n=0; n < problemID.size(); n++)
  {
*/
    cout << "digraph G {" << endl;
    cout << "   " << "node [shape=none];" << endl;
    for (int i=0; i < problems.size(); i++)
    {
      cout << "   " << problems.at(i)->getID() << " [fontname=\"Courier New\" label=< " << endl;
      cout << "     <table border=\"1\" cellspacing=\"0\" color=\"black\"> " << endl;
      cout << "       <tr>" << endl;
      cout << "         <td port=\"ID\" border=\"1\" bgcolor=\"black\" color=\"black\"><B><font color=\"white\">" << problems.at(i)->getID() << "</font></B></td>" << endl;
      cout << "         <td bgcolor=\"white\">" << problems.at(i)->getTargetRNASubseq()  << "</td>" << endl;
      cout << "       </tr> " << endl;
      cout << "     </table> >" << endl;
      cout << "   ];" << endl;
    }
    for (int i=0; i < problems.size(); i++)
    {
      problems.at(i)->printGraphviz();
    }
    cout << "}" << endl;
/*    Problem* p = getProblemByID(problemID.at(n));
    p->remove();
    deleteProblemByID(problemID.at(n), true);
    free(p);
  }*/
}

void Decomposition::solve()
{
  int iterations = 1;
  double remainingStopCrit = 0.0;

//  if (solvingFlag)
//    iterations = 10;
  vector<int> order;
  for (int i=0; i < problems.size(); i++)
    order.push_back(problems.at(i)->getID());

  for (int n=0; n < order.size(); n++)
  {
    cout << "*********************************************************************************************************************************************************" << endl;
    Problem* p = getProblemByID(order.at(n));
    cout << "SOLVING PROBLEM ID " << p->getID() << endl;
    cout << "*********************************************************************************************************************************************************" << endl;
    cout << endl;
    int ssLen = p->getTargetRNA().size();
    //if ( ( (ssLen * p->getParentsSize()) / (targetRNA.size()*1.0) >= 0.1 && (ssLen * p->getParentsSize()) / (targetRNA.size()*1.0) <= 0.9 ) || i==0 )
    double localStopCrit = 0.0;

    if (stopCrit > 0)
      localStopCrit = (stopCrit * p->getParticipation() ) * (pow(targetRNA.size(), index) / (1.0*pow(p->getTargetRNA().size(), index)) );
      //localStopCrit = (((stopCrit*1.0)/iterations) * problems.at(i)->getParticipation() ) * (pow(targetRNA.size(), index) / (1.0*pow(problems.at(i)->getTargetRNA().size(), index)) );
    else
      localStopCrit = stopCrit * p->getParticipation(); //((stopCrit*1.0)/iterations) * problems.at(i)->getParticipation();

//    if (stopCrit < 0 && p->getID() != 0 && -localStopCrit > p->getTargetRNA().size()*1.0)
//    {
//      remainingStopCrit -= localStopCrit + (p->getTargetRNA().size()*1.0);
//      localStopCrit = p->getTargetRNA().size() * (-1.0);
//    }

    if (p->getID() == 0)
      localStopCrit -= remainingStopCrit;

    cout << "TARGET:       " << p->getTargetRNA() << endl;
    cout << "SIZE:         " << ssLen << endl;
    if (stopCrit > 0)
      cout << "STOPCRIT:     " <<  localStopCrit << " local evaluations (" << (stopCrit * p->getParticipation() ) << " global evaluations)" << endl;
    else
      cout << "STOPCRIT:     " << -localStopCrit << " seconds" << endl;

    if (p->getChildrenSize() > 0)
    {
      cout << "SUBSEQUENCES: " << endl;
      for (int i=0; i<p->getChildrenSize(); i++)
      {
        cout << "  PROBLEM ID " << p->getChildrenPairAt(i).second->getID() << " (Position: "<< p->getChildrenPairAt(i).first <<"):\t" << p->getChildrenPairAt(i).second->getTargetRNA() << endl;
        cout << "    -> A total of " << p->getChildrenPairAt(i).second->getSequencesSize() << " different sequences found" << endl;
        //for (int j=0; j < p->getChildrenPairAt(i).second->getSequencesSize(); j++)
        //  cout << "        "<< p->getChildrenPairAt(i).second->getSequenceAt(j) << endl;
      }
      cout << endl;
    }

    cout << endl;
    p->solve(randomSeed, popSize, localStopCrit, solvingFlag);

    cout << p->getStopCrit() << "  " << endl;
    cout << stopCrit << endl;

    if (stopCrit > 0)
      stopCrit -= ( p->getStopCrit() /  (pow(targetRNA.size(), index) / (1.0*pow(p->getTargetRNA().size(), index)) ) );
    else
      stopCrit -= p->getStopCrit() ;

    cout << stopCrit << endl; //exit(1);

    updateParticipation();

    if (p->getSequencesSize() < 1 && p->getID() != 0)
    {
      // Delete problem and re-organize structure
      p->remove();
      deleteProblemByID(p->getID(), false);
      free(p);
    }
    p->printSequences();
    cout << "*********************************************************************************************************************************************************" << endl;
    cout << endl;// << endl << endl << endl;
  }
}

void Decomposition::printSolutions()
{
  getProblemByID(0)->printSequences();
}

Decomposition::~Decomposition()
{
  for(int i=0; i < problems.size(); i++)
    free(problems.at(i));
}
