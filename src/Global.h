#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <string>
#include <bits/stdc++.h>

extern "C"
{
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/equilibrium_probs.h>
#include <ViennaRNA/params/io.h>
}

#define MODEL_TURNER1999	0
#define MODEL_TURNER2004	1

using namespace std;

//static long rnd_uni_init;
extern bool found;
extern clock_t runtime_start;
extern int iPrint;
extern double stepTime; // Used for studying converge

#endif /* GLOBAL_H_ */
