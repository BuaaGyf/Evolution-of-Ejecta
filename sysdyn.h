//
//  ode.h
//  
//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
//
#ifndef _SYSDYN_H_
#define _SYSDYN_H_

#include "kinematics.h"
#include <sstream>
//#include <iostream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include "binary.h"
#include "particles.h"
#include "body.h"



void SolvingEqns(fem& alpha, fem& beta, poly& PolyA, poly& PolyB, DEBRIS& d, KEPORB& kep, ODE_OPTION& odeopt, SOLAR& ss, BPHASE& x, CLOUD& cl);
void BinarySIStep(double t, BPHASE &x, BPHASE &x1, GRAV &gt0, GRAV &gt1, ODE_OPTION &odeopt, fem &alpha, fem& beta, KEPORB &kep, SOLAR &ss, MASCON &varpool);
void OXYZ2SXYZ(KEPORB &kep, Vector Pos, Vector Vel);
void DebBinarySIStep(double t, Vector &POS0, Vector &VEL0, Vector &POS1, Vector &VEL1, ODE_OPTION &odeopt, BPHASE &x, fem &alpha, fem &beta, DEBRIS &d, int n, KEPORB &kep, SOLAR &ss, MASCON &varpool);
void getBinaryForce(Vector Force, MASCON &varpool, BPHASE &x, fem &alpha, fem &beta, Vector Pos, DEBRIS &d, int n, KEPORB &kep, SOLAR &ss);
void DetectCollisionState(int& ind, fem &alpha, fem &beta, poly &polya, poly &polyb, ODE_OPTION &odeopt, Vector POS1, Vector VEL1, BPHASE &x1, Vector POSout, Vector VELout);
void CollisionModel(poly &poly, ODE_OPTION &odeopt, int cfid2, Vector POSin, Vector VELin, Vector POSout, Vector VELout);
void CheckCollision(poly &p, Vector r, int& ind, int& cfid);
int chkOCC(Vector ROrb, Vector RBody, Vector Rpt, Vector Ellip, Matrix TranMX);
void getHeliosForce(Vector Force, double t, Vector Pos, Vector Vel, fem &alpha, fem &beta, DEBRIS &d, int n, KEPORB &kep, SOLAR &ss);
void DebHeliosSIStep(double t, Vector &POS0, Vector &VEL0, Vector &POS1, Vector &VEL1, ODE_OPTION &odeopt, fem &alpha, fem &beta, DEBRIS &d, int n, KEPORB &kep, SOLAR &ss);
void InitialChecking(DEBRIS &d, BPHASE &x, CLOUD &cl, fem &alpha, fem &beta, poly &polya, poly &polyb, ODE_OPTION &odeopt);
void RotationX(Matrix &m, const double theta);
void RotationY(Matrix &m, const double theta);
void RotationZ(Matrix &m, const double theta);

using namespace std;

#endif /* defined(____binary__) */
