//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
//

#ifndef _KINEMATICS_H_
#define _KINEMATICS_H_

#include "constant.h"
#include "gravity.h"
#include "body.h"
#include <iomanip>



void QuatSet(Quaternion q, double a0, double a1, double a2, double a3);
void QuatCopy(const Quaternion q1, Quaternion q2);
void QuatNorm(Quaternion q);
void QuatMultiply(const Quaternion q1, const Quaternion q2, Quaternion q);
void Quat2DCM(const Quaternion q, Matrix m);
void Quat2IvDCM(const Quaternion q, Matrix m);
//void DCM2Quat(const Matrix m, Quaternion q);
//void Euler2Quat(const double a[3], Quaternion q);

#endif
