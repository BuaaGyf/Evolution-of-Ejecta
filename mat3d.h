
#ifndef _MAT3D_H
#define _MAT3D_H

#include <math.h>
#include <assert.h>

/*
 The three-dimensional vector and matrix are specially defined here by "Vector" and "Matrix" for their generalities in our physics.
 */

typedef double Vector[3];
typedef Vector Matrix[3];

void vectorCopy(const Vector u,Vector v);

void vectorScale(const Vector u,const double s,Vector v);

void vectorAdd(const Vector v1,const Vector v2,Vector v);

void vectorSub(const Vector v1,const Vector v2,Vector v);

void vectorCross(const Vector v1,const Vector v2,Vector v);

void vectorMultiply(const Vector v1,const Vector v2,Vector v);

void vectorDivision(const Vector v1,const Vector v2,Vector v);

double vectorDot(const Vector v1,const Vector v2);

double vectorMagSq(const Vector v);

double vectorMag(const Vector v);

void vectorNorm(Vector v);

void vectorSet(Vector v,const double x,const double y,const double z);

void vectorZero(Vector v);

void vectorTransform(const Matrix m,const Vector u,Vector v);

void matrixCopy(const Matrix a,Matrix b);

void matrixZero(Matrix m);

void matrixIdentity(Matrix m);

void matrixDiagonal(const Vector v,Matrix m);

void matrixScale(const Matrix a,const double s,Matrix b);

void matrixSub(const Matrix a,const Matrix b,Matrix c);

void matrixMultiply(const Matrix a,const Matrix b,Matrix c);

void matrixTranspose(Matrix a);

double matrixDet(const Matrix a);

double matrixTrace(const Matrix a);

void matrixExtend(const Matrix m, double* a);

void matrixContract(const double* a, Matrix m);

void matrixInverse(const Matrix m, Matrix a);

void matrixNorm(const Matrix m, Matrix a);

void Matrixset(Matrix a, const double b1, const double b2, const double b3, const double b4, const double b5, const double b6, const double b7, const double b8, const double b9);

void xyz2tpr(const Vector v1, Vector v);

#endif
