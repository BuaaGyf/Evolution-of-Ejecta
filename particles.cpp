#include "particles.h"

void CloudAlloc(DEBRIS &d, CLOUD &cl)
{
    int i;
    
    cl.status = new int[d.NumDebris];
    
    cl.DebrisPos = new double*[d.NumDebris];
    for (i=0; i<d.NumDebris; i++) cl.DebrisPos[i] = new double[3];
    
    cl.DebrisVel = new double*[d.NumDebris];
    for (i=0; i<d.NumDebris; i++) cl.DebrisVel[i] = new double[3];
}

void CloudCopy(DEBRIS &d, CLOUD &cl1, CLOUD &cl)
{
    int i,j;
    
    for (i=0; i<d.NumDebris; i++)
    {
        cl.status[i] = cl1.status[i];
        for (j=0; j<3; j++)
        {
            cl.DebrisPos[i][j] = cl1.DebrisPos[i][j];
            cl.DebrisVel[i][j] = cl1.DebrisVel[i][j];
        }
    }
}


