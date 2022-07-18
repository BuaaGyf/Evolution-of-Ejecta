//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
//

#ifndef _BODY_H_
#define _BODY_H_

#include <iostream>
#include <fstream>
#include <string>
#include "mat3d.h"

//#include "kinematics.h"



using namespace std;

typedef double Quaternion[4];

typedef struct fem {

	double TotalMass;
	Vector MassCenter, InertVec;
	Vector Ellipsoid;
	Matrix Inertia;
	int NumVerts, NumFaces, NodeNum, ElemNum;
	double** Vertices;
	int** Faces;
	double* NodeDensities;
	double* NodeWeights;
	double** Nodes;
	int** Elements;
	double* ElemJS;
} FEM;

typedef struct poly{

	double DeltaAngle;
	int MeshNum;
	int VerticeNum;
	int FaceNum;
	int HeadIndexNum;
	int BodyLeftIndexRowNum;
	int BodyLeftIndexColumnNum;
	double** Vertices;
	int** Faces;
	int* HeadIndex;
	int** BodyLeftIndex;
	int** BodyRightIndex;
	int* TailIndex;
	double EnvelopeRadiusSqr;
	double **FaceNormVec;
	int* FaceCollisionFlag;

} POLY;



typedef struct debris {
    int NumDebris;
    int *Id;
    double *Mass, *Radius;

} DEBRIS;

typedef struct bphase {
	//
	Vector AlphaPos; // OXYZ
	Vector AlphaVel; // OXYZ
	Vector AlphaMomentum; // OXYZ
	Matrix DCM_A; // AXYZ->AXaYaZa
	Vector AlphaAngVel; //AXaYaZa
	Vector AlphaAngularmomentum; //AXaYaZa
	//
	Vector BetaPos; // OXYZ
	Vector BetaVel; // OXYZ
	Vector BetaMomentum; // OXYZ
	Matrix DCM_B;; // BXYZ->BXbYbZb
	Vector BetaAngVel; // BXbYbZb
	Vector BetaAngularmomentum;  //BXaYaZa
	//
} BPHASE;

typedef struct grav {
	Vector GravForceA;
	Vector GravForceB;
	Vector GravTorqueA;
	Vector GravTorqueB;
} GRAV;


typedef struct mascon {

	double** alpha, ** beta;

} MASCON; // ode solver options




void LoadFEM(string filename, FEM& fe);
void AllocFEM(FEM& fe);
void LoadPoly(string filename,POLY& poly);
void PolyAlloc(POLY &poly);
void LoadDebris(string debrisfile, DEBRIS &d);
void DebrisAlloc(DEBRIS &d);

#endif
