#include "Problem.h"
#include "Nsga2.h"

using namespace std;

Problem::Problem()
{
  stopCrit = 0;
  attempted = false;
}

int Problem::getID()
{
  return ID;
}

void Problem::setID(int id)
{
  ID = id;
}

double Problem::getStopCrit()
{
  return stopCrit;
}

void Problem::setStopCrit(int stop)
{
  stopCrit = stop;
}

bool Problem::isAttempted()
{
  return attempted;
}

string Problem::getTargetRNA()
{
  return targetRNA;
}

string Problem::getMask()
{
  return mask;
}

void Problem::setMask(string m)
{
  mask = m;
}

void Problem::setTargetRNA(string s, double occupation_)
{
  targetRNA = s;
  mask = s;
  occupation = occupation_;
  attempted = false;
  children.clear();
  parents.clear();
}

double Problem::getOccupation()
{
  return occupation;
}

double Problem::getParticipation()
{
  return participation;
}

void Problem::setParticipation(double participation_)
{
  participation = participation_;
}

vector< pair <int, Problem*> > Problem::getParents()
{
  return parents;
}

void Problem::addParent(int pos, Problem* p)
{
  parents.push_back( make_pair(pos,p) );
}

void Problem::deleteParent(int id)
{
  for (int i=0; i < parents.size(); i++)
  {
    if (parents.at(i).second->getID() == id)
      parents.erase(parents.begin()+i);
  }
}

Problem* Problem::getParentAt(int i)
{
  return parents.at(i).second;
}

pair<int, Problem*> Problem::getParentPairAt(int i)
{
  return parents.at(i);
}

int Problem::getParentsSize()
{
  return parents.size();
}

vector< pair<int, Problem*> > Problem::getChildren()
{
  return children;
}

int Problem::getChildrenSize()
{
  return children.size();
}

Problem* Problem::getChild(int id)
{
  for(int i=0; i < children.size(); i++)
    if (children.at(i).second->getID() == id)
      return children.at(i).second;
}

pair<int, Problem*> Problem::getChildrenPairAt(int i)
{
  return children.at(i);
}

void Problem::addChild(int pos, Problem* p)
{
  children.push_back( make_pair(pos,p) );
  string asterisk (p->getTargetRNA().size(),'*');
  string newMask = mask.substr(0,pos) + asterisk + mask.substr(pos+asterisk.size(), mask.size() - pos - asterisk.size());
  mask = newMask;
}

string Problem::getCompressedTargetRNASubseq()
{
  string targetRNASubseq = targetRNA;
  int rem = 0;
  for(int i=0; i < children.size(); i++)
  {
     int pos = children.at(i).first - rem;
     Problem* p = children.at(i).second;
     string newStr = targetRNASubseq;
     string sID = "$" + std::to_string(p->getID());
     newStr.replace(pos, p->getTargetRNA().size(), sID);
     rem += p->getTargetRNA().size() - sID.size();
     targetRNASubseq = newStr;
  }
  return targetRNASubseq;
}

string Problem::getTargetRNASubseq()
{
  string targetRNASubseq = targetRNA;
  int rem = 0;
  for(int i=0; i < children.size(); i++)
  {
     int pos = children.at(i).first - rem;
     Problem* p = children.at(i).second;
     string newStr = targetRNASubseq;
     string sID = "</td>\n         <td port=\"P" + std::to_string(i) + "\" bgcolor=\"gray\"><B><font color=\"black\">" + std::to_string(p->getID()) + "</font></B></td>\n         <td bgcolor=\"white\">";
     newStr.replace(pos, p->getTargetRNA().size(), sID);
     rem += p->getTargetRNA().size() - sID.size();
     targetRNASubseq = newStr;

/*     string sID = "<B>#" + std::to_string(p->getID())+"</B>";
     string newStr = targetRNASubseq.substr(0, pos) + sID + targetRNASubseq.substr(pos+p->getTargetRNA().size(), targetRNASubseq.size() - pos - p->getTargetRNA().size() );
     rem += p->getTargetRNA().size() - sID.size();
     targetRNASubseq = newStr;*/
  }
  return targetRNASubseq;
}

void Problem::deleteChild(int id)
{
  for (int i=0; i < children.size(); i++)
  {
    if (children.at(i).second->getID() == id)
      children.erase(children.begin()+i);
  }
}

vector<string>* Problem::getSequences()
{
  return &sequences;
}

int Problem::getSequencesSize()
{
  return sequences.size();
}

string Problem::getSequenceAt(int i)
{
  return sequences.at(i);
}

void Problem::addSequence(string s)
{
  sequences.push_back(s);
}

void Problem::deleteSequence(int i)
{
  sequences.erase(sequences.begin()+i);
}

void Problem::print()
{
  cout << "   ID:         " << ID << endl;

  cout << "   Size:       " << targetRNA.size() << endl;

  cout << "   Problem:    " << targetRNA << endl;

  cout << "   Mask:       " << mask << endl;

  cout << "   Compressed: " << getCompressedTargetRNASubseq() << endl;

  cout << "   Occupation: " << occupation*100 << "%"<<endl;

  cout << "   Participation: " << participation*100 << "%"<<endl;

  cout << "   Children:   " << endl;
  if (children.size() == 0)
    cout << "   \t{empty}" << endl;
  for (int i=0; i < children.size(); i++)
    cout << "   \t[Pos: " << children.at(i).first << "]\t" << children.at(i).second->getTargetRNA() << "\t(#" << children.at(i).second->getID() << ")" << endl;

  cout << "   Parents:    " << endl;
  if (parents.size() == 0)
    cout << "   \t{empty}" << endl;
  for (int i=0; i < parents.size(); i++)
    cout << "   \t[Pos: " << parents.at(i).first << "]\t(#" << parents.at(i).second->getID() << ")" << endl;

  cout << "   Sequences:  " << endl;
  if (sequences.size() == 0)
    cout << "   \t{empty}" << endl;
  for (int i=0; i < sequences.size(); i++)
    cout << "\t" << sequences.at(i) << endl;
}

void Problem::printGraphviz()
{
  for (int i=0; i < children.size(); i++)
  {
//    cout << "   " << ID << ":ID->" << children.at(i).second->getID() << ":P" << i << endl;
    cout << "   " << ID << ":P" << i << "->" << children.at(i).second->getID() << ":ID" << endl;
  }
}

void Problem::printSequences()
{
  sortAndUnique();
  for (int i=0; i < sequences.size(); i++)
    cout << "#" << i << "#\t" << sequences.at(i) << endl;
}

void Problem::solve(int randomSeed, int popSize, double stopCrit_, bool solvingFlag)
{
  attempted = true;

  Nsga2* algorithm = NULL;

  stopCrit = stopCrit_;

  algorithm = new Nsga2(randomSeed, popSize, stopCrit, solvingFlag, this);
  clock_t t, t2;
  t = clock();
  sequences = algorithm->run();
  t2 = clock() - t;

  if (stopCrit < 0)
    stopCrit = -((double) t2 / CLOCKS_PER_SEC);
  else
    stopCrit = algorithm->getNEvals();

cout << stopCrit << endl;
//  cout << "** NEVALS ** " << algorithm->getNEvals() << endl;
//  exit(1);

  sortAndUnique();
//  limitSequences(popSize);

  delete algorithm;
}

int Problem::getEditDistance(std::string first, std::string second)
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
 
    return T[m][n];
}

double Problem::findStringSimilarity(std::string first, std::string second) {
    double max_length = std::max(first.length(), second.length());
    if (max_length > 0) {
        return (max_length - getEditDistance(first, second)) / max_length;
    }
    return 1.0;
}

/*
void Problem::limitSequences(int n)
{
  vector<string> aux;
  if (sequences.size() > n)
  {
    int step = sequences.size() / n;
    for (int i=0; i<sequences.size(); i+=step)
    {
	aux.push_back(sequences.at(i));
    }
    sequences.clear();
    for (int i=0; i<aux.size(); i++)
      sequences.push_back(aux.at(i));
  }
}
*/

void Problem::limitSequences (int n)
{
  if (sequences.size() > n)
  {
    vector<pair<double,string>> seq_sim;
    for (int i=0; i < sequences.size(); i++)
    {
      pair<double,string> p;
      p.first = sequences.at(i).size();
      p.second = sequences.at(i);
      seq_sim.push_back(p);
    }
    for (int i=0; i < seq_sim.size(); i++)
    {
      for (int j=i+1; j < seq_sim.size(); j++)
      {
        int v = 0;
        for (int k=0; k < seq_sim[i].second.size(); k++)
        {
          if (seq_sim[i].second[k] == seq_sim[j].second[k])
            v++;
        }
	seq_sim[i].first += v;
	seq_sim[j].first += v;
//	if (v < seq_sim[i].first) seq_sim[i].first = v;
//	if (v < seq_sim[j].first) seq_sim[j].first = v;
      }
      seq_sim[i].first /= ( (seq_sim.size() - 1) * 1.0 );
    }
    sort( seq_sim.begin(), seq_sim.end() );
    seq_sim.erase( unique( seq_sim.begin(), seq_sim.end() ), seq_sim.end() );

    sequences.clear();
    for (int i=0; i < n; i++)
    {
      cout << "Selected: " << seq_sim[i].second << "   " << seq_sim[i].first << endl;
      sequences.push_back(seq_sim[i].second);
    }
  }
}

/*void Problem::limitSequences(int n)
{
  if (sequences.size() > n)
  {
    vector<pair<double,string>> seq_sim;
    for (int i=0; i<sequences.size(); i++)
    {
      pair<double,string> p;
      p.first = 0;
      p.second = sequences.at(i);
      seq_sim.push_back(p);
    }
    for (int i=0; i < seq_sim.size(); i++)
    {
      for (int j=i+1; j < seq_sim.size(); j++)
      {
	  double v =findStringSimilarity(seq_sim.at(i).second, seq_sim.at(j).second);
	  seq_sim.at(i).first += v;
	  seq_sim.at(j).first += v;
      }
      seq_sim.at(i).first /= (seq_sim.size()-1);
    }
    sort( seq_sim.begin(), seq_sim.end() );
    seq_sim.erase( unique( seq_sim.begin(), seq_sim.end() ), seq_sim.end() );

    sequences.clear();
    for (int i=0; i < n; i++)
      sequences.push_back(seq_sim[i].second);
  }
}*/

void Problem::remove()
{
  for (int i=0; i < parents.size(); i++)
  {
    // Modify targetRNA of its parents
    int pos = parents.at(i).first;
    Problem* parent = parents.at(i).second;
    parent->deleteChild(ID);
    for (int j=0; j<children.size(); j++)
    {
      parent->addChild(pos + children.at(j).first, children.at(j).second);
      children.at(j).second->deleteParent(ID);
      children.at(j).second->addParent(pos+children.at(j).first, parent);
    }
  }
}

void Problem::sortAndUnique()
{
  sort( sequences.begin(), sequences.end() );
  sequences.erase( unique( sequences.begin(), sequences.end() ), sequences.end() );
}

void Problem::sortChildrenParents()
{
  sort( children.begin(), children.end() );
  sort( parents.begin(), parents.end() );
}

Problem::~Problem()
{
}
