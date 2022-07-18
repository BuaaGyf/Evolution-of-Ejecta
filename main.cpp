//
//  main.cpp
//  
//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
//

#include <iostream>
#include "body.h"
#include "initialization.h"
#include "kinematics.h"
#include "binary.h"
#include "odeset.h"
#include "sysdyn.h"
#include "particles.h"
#include "solar.h"

//

using namespace std;
//
int main(int argc,char *argv[])
{
	FEM alpha, beta;
    POLY polya, polyb;
    DEBRIS d;
    BPHASE x;
    ODE_OPTION odeopt;
    KEPORB kep;
    SOLAR ss;
    CLOUD cl;
    
    setbuf(stdout,(char *)NULL);
    
    if (argc!=2)
    {
        (void) fprintf(stderr,"Usage: %s file\n",argv[0]);
        exit(1);
    }
    
    Initialization(argv[1], alpha, beta, polya, polyb, d, odeopt, ss, kep, x, cl);
    
    clock_t systim1;
    systim1 = clock();

    SolvingEqns(alpha, beta, polya, polyb, d, kep, odeopt, ss, x, cl);
    
    systim1 = clock() - systim1;
    cout<<"Simulation time = "<<systim1<<"/"<<CLOCKS_PER_SEC<<" sec"<<endl;
            
    return 0;
    
}


