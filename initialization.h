//
//  initialization.h
//  
//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
//

#ifndef _INITIALIZATION_H_
#define _INITIALIZATION_H_


#include <algorithm>
#include "binary.h"
#include "solar.h"
#include "particles.h"
#include "sysdyn.h"

using namespace std;

void Initialization(char* parafile, fem &alpha, fem &beta, poly &polya, poly &polyb, DEBRIS &d, ODE_OPTION &odeopt, SOLAR &ss, KEPORB &kep, BPHASE &x, CLOUD &cl);
void LoadParas(char* parafile, string &polyfileA, string &polyfileB, string &polyfileA1, string &polyfileB1, string &debrisfile, string &iniconfile, ODE_OPTION &odeopt, SOLAR &ss);
void LoadInicon(string iniconfile, DEBRIS &d, KEPORB &kep, BPHASE &x, CLOUD &cl); 
void OriginCorrection(fem &alpha, fem &beta, BPHASE &x); 

#endif /* */
