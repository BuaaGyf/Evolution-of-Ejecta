//
//  binary.cpp
//  
//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
// In the simulation of the bianry motion, we only include the perturbation solar tide, the SRP is not reffered here because for the massive objects of asteroids,
// the SRP acceleration is too small and neglectable.
//
//

#include "binary.h"

void BinaryForceTorque(Vector &GravForceA, Vector &GravTorqueA, Vector &GravTorqueB, fem &alpha, fem &beta, MASCON &varpool, KEPORB &kep, SOLAR &ss)
{
    // vectorZero(GravForceA); // attraction on Alpha
    // vectorZero(GravTorqueA); // torque on Alpha
    // vectorZero(GravTorqueB); // torque on Beta

    
	
    for (int ii=0; ii<alpha.NodeNum; ii++)
    {
        double dd;
        Vector VV1,VV2,VV3,VV4;
        Vector GFAtmp,GTAtmp,GTBtmp;
        
        vectorZero(GFAtmp); //
        vectorZero(GTAtmp); //
        vectorZero(GTBtmp); //
        vectorSet(VV2,varpool.alpha[ii][3],varpool.alpha[ii][4],varpool.alpha[ii][5]);
        for (int jj=0; jj<beta.NodeNum; jj++)
        {
            vectorSet(VV3,varpool.beta[jj][0]-varpool.alpha[ii][0],varpool.beta[jj][1]-varpool.alpha[ii][1], \
                      varpool.beta[jj][2]-varpool.alpha[ii][2]);
            dd = vectorMag(VV3);
            dd = G*alpha.NodeWeights[ii]*beta.NodeWeights[jj]/(dd*dd*dd);
            vectorScale(VV3,dd,VV1);
            vectorAdd(GFAtmp,VV1,GFAtmp);
            vectorCross(VV2,VV1,VV4);
            vectorAdd(GTAtmp,VV4,GTAtmp);
            vectorSet(VV3,varpool.beta[jj][3],varpool.beta[jj][4],varpool.beta[jj][5]);
            vectorCross(VV1,VV3,VV4);
            vectorAdd(GTBtmp,VV4,GTBtmp);
        }
        
 
		vectorAdd(GravForceA,GFAtmp,GravForceA);
		vectorAdd(GravTorqueA,GTAtmp,GravTorqueA);
		vectorAdd(GravTorqueB,GTBtmp,GravTorqueB);
    }


}



void getMomentumMoment(fem& alpha , fem& beta, BPHASE& x, Vector MM)
{
    Vector MMT_A, MMR_A, MMT_B, MMR_B;
    Vector tmpV1;

    vectorCross(x.AlphaPos, x.AlphaMomentum, MMT_A);
    vectorCross(x.BetaPos, x.BetaMomentum, MMT_B);

    vectorScale(MMT_A, -1, MMT_A);
    vectorScale(MMT_B, -1, MMT_B);

    vectorTransform(x.DCM_A, x.AlphaAngularmomentum, MMR_A);
    vectorTransform(x.DCM_B, x.BetaAngularmomentum, MMR_B);
    //
    vectorAdd(MMT_A, MMT_B, tmpV1);
    vectorAdd(MMR_A, MMR_B, MM);
    //
    vectorAdd(tmpV1, MM, MM); // in OXYZ
}

void getMomentum(fem &alpha, fem &beta, BPHASE &x, Vector MT)
{

    vectorAdd(x.AlphaMomentum, x.BetaMomentum, MT); // in OXYZ
}

void getMechanicEnergy(fem &alpha, fem &beta, BPHASE &x,  MASCON &varpool, double &E)
{
    int i, j;
    double KT_A, KR_A, KT_B, KR_B, U;
    double tmpd1;
	Vector tmpV1, tmpV2;
	//Vector AccA, AccB;
	Matrix DCM_A1, IvDCM_A1, DCM_B1, IvDCM_B1;
	BPHASE x1;

    matrixCopy(x.DCM_A, DCM_A1);
    matrixCopy(x.DCM_B, DCM_B1);

    matrixInverse(DCM_A1, IvDCM_A1);
    matrixInverse(DCM_B1, IvDCM_B1);


    KT_A = 0.5 * (x.AlphaMomentum[0]*x.AlphaMomentum[0]+x.AlphaMomentum[1]*x.AlphaMomentum[1]+x.AlphaMomentum[2]*x.AlphaMomentum[2])/alpha.TotalMass;
    KT_B = 0.5 * (x.BetaMomentum[0]*x.BetaMomentum[0]+x.BetaMomentum[1]*x.BetaMomentum[1]+x.BetaMomentum[2]*x.BetaMomentum[2])/beta.TotalMass;

    KR_A = 0.5 * (x.AlphaAngularmomentum[0] * x.AlphaAngularmomentum[0] / alpha.InertVec[0] + x.AlphaAngularmomentum[1] * x.AlphaAngularmomentum[1] / alpha.InertVec[1] +x.AlphaAngularmomentum[2] * x.AlphaAngularmomentum[2] / alpha.InertVec[2]);
    KR_B = 0.5 * (x.BetaAngularmomentum[0] * x.BetaAngularmomentum[0] / beta.InertVec[0] + x.BetaAngularmomentum[1] * x.BetaAngularmomentum[1] / beta.InertVec[1] + x.BetaAngularmomentum[2] * x.BetaAngularmomentum[2] / beta.InertVec[2]);
    //
    U = 0.0;
	//
	phaseCopy(x, x1);
	//

	//
	//
	for (i = 0; i < alpha.NodeNum; i++)
	{
		vectorSet(tmpV1, alpha.Nodes[i][0], alpha.Nodes[i][1], alpha.Nodes[i][2]);
		vectorTransform(DCM_A1, tmpV1, tmpV2);
		vectorAdd(x.AlphaPos, tmpV2, tmpV1);
		for (j = 0; j < 3; j++)
		{
			varpool.alpha[i][j] = tmpV1[j];
			varpool.alpha[i][j + 3] = tmpV2[j];
		}
	}
	//
	for (i = 0; i < beta.NodeNum; i++)
	{
		vectorSet(tmpV1, beta.Nodes[i][0], beta.Nodes[i][1], beta.Nodes[i][2]);
		vectorTransform(DCM_B1, tmpV1, tmpV2);
		vectorAdd(x.BetaPos, tmpV2, tmpV1);
		for (j = 0; j < 3; j++)
		{
			varpool.beta[i][j] = tmpV1[j];
			varpool.beta[i][j + 3] = tmpV2[j];
		}
	}
	//
	
	for (int ii = 0; ii < alpha.NodeNum; ii++)
	{
		double dd;
		Vector VV2, VV3;

		vectorSet(VV2, varpool.alpha[ii][3], varpool.alpha[ii][4], varpool.alpha[ii][5]);
		for (int jj = 0; jj < beta.NodeNum; jj++)
		{
			vectorSet(VV3, varpool.beta[jj][0] - varpool.alpha[ii][0], varpool.beta[jj][1] - varpool.alpha[ii][1], \
				varpool.beta[jj][2] - varpool.alpha[ii][2]);
			dd = vectorMag(VV3);
			dd = G * alpha.NodeWeights[ii] * beta.NodeWeights[jj] / dd;
			tmpd1 = sqrt(dd * dd);
			U = U + tmpd1;
		}
	}
	

	//
    E = KT_A + KR_A + KT_B + KR_B - U;
	double flagx = 0;
}




