//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
//

#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include "stdlib.h"
#include <sstream>
#include "body.h"
//#include <omp.h>

using namespace std;

typedef struct cloud {
    //
    int *status; // -2: on Alpha, -1: on Beta, 0: on Orbit, 1: escaped to heliocentric orbit (out of the Hill sphere)
    double **DebrisPos, **DebrisVel;
    
} CLOUD;

void CloudAlloc(DEBRIS &d, CLOUD &cl);
void CloudCopy(DEBRIS &d, CLOUD &cl1, CLOUD &cl);

#endif
