//
//  binary.h
//  
//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
//

#ifndef _BINARY_H_
#define _BINARY_H_

#include "stdlib.h"
#include "solar.h"
#include "kinematics.h"
#include "sysdyn.h"
#include "body.h"

#include <fstream>
#include <iomanip>
#include <omp.h>



void BinaryForceTorque(Vector &GravForceA, Vector &GravTorqueA, Vector &GravTorqueB, fem &alpha, fem &beta, MASCON &varpool, KEPORB &kep, SOLAR &ss);
void phaseCopy(BPHASE &x1, BPHASE &x2);
void gravCopy(GRAV &gt1, GRAV &gt0);
void phaseExtract(BPHASE &x, double* y);
void phaseCompress(const double* y, BPHASE &x);
//
void getMomentumMoment(fem &alpha, fem &beta, BPHASE &x, Vector MM);
void getMomentum(fem &alpha, fem &beta, BPHASE &x, Vector MT);
void getMechanicEnergy(fem &alpha, fem &beta, BPHASE& x, MASCON &varpool, double &E);

#endif /* defined(____binary__) */
