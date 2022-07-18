// Created by Yang Yu and Yunfeng Gao on 2022-5-1.


#include "sysdyn.h"
#include "mat3d.h"
#include "math.h"

void SolvingEqns(fem& alpha, fem& beta, poly& polya, poly& polyb, DEBRIS& d, KEPORB& kep, ODE_OPTION& odeopt, SOLAR& ss, BPHASE& x, CLOUD& cl)
{
	int i, j, count; //ID
	double t;
	//
	ofstream resfile;
	BPHASE x1;
    MASCON varpool;
    GRAV gt0,gt1;

    Matrix DCM_A1, IvDCM_A1, DCM_B1, IvDCM_B1;

    Vector tmpV1, tmpV2;
    

    varpool.alpha = new double* [alpha.NodeNum];
	for (i = 0; i < alpha.NodeNum; i++) varpool.alpha[i] = new double[6];
	varpool.beta = new double* [beta.NodeNum];
	for (i = 0; i < beta.NodeNum; i++) varpool.beta[i] = new double[6];


    matrixCopy(x.DCM_A, DCM_A1);
    matrixCopy(x.DCM_B, DCM_B1);

    matrixInverse(DCM_A1, IvDCM_A1);
    matrixInverse(DCM_B1, IvDCM_B1);


    for (i=0; i<alpha.NodeNum; i++)
    {
        vectorSet(tmpV1,alpha.Nodes[i][0],alpha.Nodes[i][1],alpha.Nodes[i][2]);
        vectorTransform(DCM_A1, tmpV1, tmpV2);
        vectorAdd(x.AlphaPos,tmpV2,tmpV1);
        for (j=0; j<3; j++)
        {
            varpool.alpha[i][j] = tmpV1[j];  // OXYZ
            varpool.alpha[i][j+3] = tmpV2[j]; // AXYZ
        }
    }
    //
    for (i=0; i<beta.NodeNum; i++)
    {
        vectorSet(tmpV1,beta.Nodes[i][0],beta.Nodes[i][1],beta.Nodes[i][2]);
        vectorTransform(DCM_B1, tmpV1, tmpV2);
        vectorAdd(x.BetaPos,tmpV2,tmpV1);
        for (j=0; j<3; j++)
        {
            varpool.beta[i][j] = tmpV1[j]; // OXYZ
            varpool.beta[i][j+3] = tmpV2[j]; //BXYZ
        }
    }

    vectorZero(gt0.GravForceA); //
    vectorZero(gt0.GravForceB); //
    vectorZero(gt0.GravTorqueA); //
    vectorZero(gt0.GravTorqueB); //


    BinaryForceTorque(gt0.GravForceA, gt0.GravTorqueA, gt0.GravTorqueB, alpha, beta, varpool, kep, ss);

    vectorScale(gt0.GravForceA, -1.0, gt0.GravForceB); // attraction on Beta in OXYZ
    // torque on Beta transfered to AXaYaZa
    vectorTransform(IvDCM_A1, gt0.GravTorqueA, gt0.GravTorqueB);
    // torque on Beta transfered to BXbYbZb
    vectorTransform(IvDCM_B1, gt0.GravTorqueA, gt0.GravTorqueB);

    vectorScale(gt0.GravTorqueA, -1.0, gt0.GravTorqueA);
    vectorScale(gt0.GravTorqueB, -1.0, gt0.GravTorqueB);



	//
	// if FunOption = 0, run the binary conservation test module, and SolarTide perturbation will be forced off
	if (odeopt.FunOption == 0)
	{
		Vector MM, MT;
		double E;

        
        ostringstream convert;
		string Filename;

		ss.SolarTide = 0;

        count = 0;
		//
		// resfile.open("BinaryConservations.bt");
		// resfile << setiosflags(ios::scientific) << setprecision(PrecDouble); // output precision and format
		// resfile << "           Time         |                      Moment of momentum (OXYZ)                     |" \
		// 	<< "                           Momentum (OXYZ)                          |    Mechanical Energy" << endl;
		//
		for (i = odeopt.StartStep; i < odeopt.EndStep; i++)
		{
			t = i * odeopt.StepSize;
			//
			getMomentumMoment(alpha, beta, x, MM);
			getMomentum(alpha, beta, x, MT);

			// phaseExtract(x, y3);
			getMechanicEnergy(alpha, beta, x, varpool, E);
			// resfile << setw(WidthDouble) << t << "  |" << setw(WidthDouble) << MM[0] << setw(WidthDouble) << MM[1] << setw(WidthDouble) << MM[2] << "  |" \
			// 	<< setw(WidthDouble) << MT[0] << setw(WidthDouble) << MT[1] << setw(WidthDouble) << MT[2] << "  |" << setw(WidthDouble) << E << endl;
			//
			BinarySIStep(t, x, x1, gt0, gt1, odeopt, alpha, beta, kep, ss, varpool);
			//

            gravCopy(gt1, gt0);
			phaseCopy(x1, x);

            count++;
				//
				//cout<<kep.Position[0]<<"  "<<kep.Position[1]<<"  "<<kep.Position[2]<<endl;
				//
				// output module
            printf("count: %d\n", count);
            if (count == odeopt.OutputInterval)
            {
                //output (x, kep, cl)
                Filename = "res.";
                convert<<setw(10)<<setfill('0')<<i+1;
                Filename.append(convert.str());
                Filename.append(".bt");
                //
                const char *outfile = Filename.c_str();
                //
                resfile.open(outfile);
                resfile<<setiosflags(ios::scientific)<<setprecision(PrecDouble); // output precision and format
                resfile<<setw(WidthDouble)<<kep.SemiMajorAxis<<setw(WidthDouble)<<kep.Eccentricity<<setw(WidthDouble)<<kep.LongAscendNode \
                <<setw(WidthDouble)<<kep.Inclination<<setw(WidthDouble)<<kep.ArgPeriapsis<<setw(WidthDouble)<<kep.MeanAnomaly<<endl;
                resfile<<setw(WidthInt)<<ID_Alpha<<setw(WidthDouble)<<x.AlphaPos[0]<<setw(WidthDouble)<<x.AlphaPos[1]<<setw(WidthDouble)<<x.AlphaPos[2] \
                <<setw(WidthDouble)<<x.AlphaVel[0]<<setw(WidthDouble)<<x.AlphaVel[1]<<setw(WidthDouble)<<x.AlphaVel[2] \
                <<setw(WidthDouble)<<x.DCM_A[0][0]<<setw(WidthDouble)<<x.DCM_A[0][1]<<setw(WidthDouble)<<x.DCM_A[0][2] \
                <<setw(WidthDouble)<<x.DCM_A[1][0]<<setw(WidthDouble)<<x.DCM_A[1][1]<<setw(WidthDouble)<<x.DCM_A[1][2] \
                <<setw(WidthDouble)<<x.DCM_A[2][0]<<setw(WidthDouble)<<x.DCM_A[2][1]<<setw(WidthDouble)<<x.DCM_A[2][2] \
                <<setw(WidthDouble)<<x.AlphaAngVel[0]<<setw(WidthDouble)<<x.AlphaAngVel[1]<<setw(WidthDouble)<<x.AlphaAngVel[2]<<endl;
                resfile<<setw(WidthInt)<<ID_Beta<<setw(WidthDouble)<<x.BetaPos[0]<<setw(WidthDouble)<<x.BetaPos[1]<<setw(WidthDouble)<<x.BetaPos[2] \
                <<setw(WidthDouble)<<x.BetaVel[0]<<setw(WidthDouble)<<x.BetaVel[1]<<setw(WidthDouble)<<x.BetaVel[2] \
                <<setw(WidthDouble)<<x.DCM_B[0][0]<<setw(WidthDouble)<<x.DCM_B[0][1]<<setw(WidthDouble)<<x.DCM_B[0][2] \
                <<setw(WidthDouble)<<x.DCM_B[1][0]<<setw(WidthDouble)<<x.DCM_B[1][1]<<setw(WidthDouble)<<x.DCM_B[1][2] \
                <<setw(WidthDouble)<<x.DCM_B[2][0]<<setw(WidthDouble)<<x.DCM_B[2][1]<<setw(WidthDouble)<<x.DCM_B[2][2] \
                <<setw(WidthDouble)<<x.BetaAngVel[0]<<setw(WidthDouble)<<x.BetaAngVel[1]<<setw(WidthDouble)<<x.BetaAngVel[2]<<endl;

                resfile<<setw(WidthDouble)<<E<<setw(WidthDouble)<<MT[0]<<setw(WidthDouble)<<MT[1]<<setw(WidthDouble)<<MT[2] \
                <<setw(WidthDouble)<<MM[0]<<setw(WidthDouble)<<MM[1]<<setw(WidthDouble)<<MM[2]<<endl;
                //
                //output cl // resfile<<endl;
                resfile.close();
                //
                convert.str("");
                count = 0;
                //
            }
			//
		}
		//
		resfile.close();
	}
	else
	{
		//cout<<"1"<<endl;
		KEPORB kep1;
		CLOUD cl1;
		//
		ostringstream convert;
		string Filename;

        
		//
		CloudAlloc(d, cl1);
		count = 0;
        int outputnum = 0;
        int daynum = 0;
        int numview = 0;
		//
		//cout<<"2"<<endl;
		for (i = odeopt.StartStep; i < odeopt.EndStep; i++)
		{
			t = i * odeopt.StepSize;

            
            if (numview == 100)
            {
                printf("The step is: %d\n",i);
                numview = 0;
            }

            numview++;


            
			//
			//-----------omp-for-loop-1-start----------------------
			int j, k;
            int orbitnum = 0;

            int ii, jj;

            GRAV gt2;        


            MASCON varpool3;

            BinarySIStep(t, x, x1, gt0, gt1, odeopt, alpha, beta, kep, ss, varpool);

            //
            KepCopy(kep, kep1);
            KepStepping(odeopt, kep1);

            gravCopy(gt1,gt2);

            matrixCopy(x1.DCM_B,DCM_B1); // DCM_B : Oxyz to OXYZ
            matrixCopy(x1.DCM_A,DCM_A1);
            matrixInverse(DCM_B1, IvDCM_B1);

            varpool3.alpha = new double* [alpha.NodeNum];
            for (ii = 0; ii < alpha.NodeNum; ii++) varpool3.alpha[ii] = new double[6];
            varpool3.beta = new double* [beta.NodeNum];
            for (ii = 0; ii < beta.NodeNum; ii++) varpool3.beta[ii] = new double[6];

            
            for (ii=0; ii<alpha.NodeNum; ii++)
                {
                    vectorSet(tmpV1,alpha.Nodes[ii][0],alpha.Nodes[ii][1],alpha.Nodes[ii][2]);
                    vectorTransform(DCM_A1, tmpV1, tmpV2);
                    vectorAdd(x1.AlphaPos,tmpV2,tmpV1);
                    for (jj=0; jj<3; jj++)
                    {
                        varpool3.alpha[ii][jj] = tmpV1[jj];  // OXYZ
                        varpool3.alpha[ii][jj+3] = tmpV2[jj]; // AXYZ
                    }
                }
                //
            for (ii=0; ii<beta.NodeNum; ii++)
                {
                    vectorSet(tmpV1,beta.Nodes[ii][0],beta.Nodes[ii][1],beta.Nodes[ii][2]);
                    vectorTransform(DCM_B1, tmpV1, tmpV2);
                    vectorAdd(x1.BetaPos,tmpV2,tmpV1);
                    for (jj=0; jj<3; jj++)
                    {
                        varpool3.beta[ii][jj] = tmpV1[jj]; // OXYZ
                        varpool3.beta[ii][jj+3] = tmpV2[jj]; //BXYZ
                    }
                }

            // vectorZero(gt2.GravForceA);
            // vectorZero(gt2.GravForceB);
            // vectorZero(gt2.GravTorqueA);
            // vectorZero(gt2.GravTorqueB);

            // BinaryForceTorque(gt2.GravForceA, gt2.GravTorqueA, gt2.GravTorqueB, alpha, beta, varpool3, kep, ss);
            // vectorScale(gt2.GravForceA, -1, gt2.GravForceB);       

            Vector domg_beta;

            domg_beta[0] = gt2.GravTorqueA[0]/alpha.InertVec[0]*(-1);
            domg_beta[1] = gt2.GravTorqueA[1]/alpha.InertVec[1]*(-1);
            domg_beta[2] = gt2.GravTorqueA[2]/alpha.InertVec[2]*(-1);






            # pragma omp parallel shared(x1,kep1,cl1) private(k)
            {
            # pragma omp for schedule (dynamic) //private ( ID )
           

			for (j = -1; j < d.NumDebris; j++)
			{
                Vector POS0, VEL0, POS1, VEL1;
                

				if (j == -1)
				{
	        		// BinarySIStep(t, x, x1, gt0, gt1, odeopt, alpha, beta, kep, ss, varpool);
      
					// //
					// KepCopy(kep, kep1);
					// KepStepping(odeopt, kep1);
                    double yyyy = 1.1111111;

                    // printf("Tasknum is %d, threadnum =  %d\n", j, omp_get_thread_num());
                    // printf("Which is the first?");
				}
				else
				{
                    // printf("Tasknum is %d, threadnum =  %d\n", j, omp_get_thread_num());
					vectorSet(POS0, cl.DebrisPos[j][0], cl.DebrisPos[j][1], cl.DebrisPos[j][2]);
					vectorSet(VEL0, cl.DebrisVel[j][0], cl.DebrisVel[j][1], cl.DebrisVel[j][2]);

					//

                    
					if (cl.status[j] == 0)
					{
                        //POS: pos of debris (OXYZ)
                        //VEL: vel of debris (OXYZ)
                        
						DebBinarySIStep(t, POS0, VEL0, POS1, VEL1, odeopt, x, alpha, beta, d, j, kep, ss, varpool); // OXYZ

                        orbitnum = orbitnum + 1;
                        
						cl1.status[j] = 2;
					}
					else if (cl.status[j] == 1)
					{
						//DebHeliosSIStep(t, POS0, VEL0, POS1, VEL1, odeopt, alpha, beta, d, j, kep, ss); //SXYZ
                        vectorCopy(POS0, POS1);
                        vectorCopy(VEL0, VEL1);
						cl1.status[j] = 1;
					}
					else if (cl.status[j] == -1)
					{
						cl1.status[j] = -1; //Fixed frames, Bxayaza
						vectorCopy(VEL0, VEL1);

                        Vector PosDBtmp;

                        Matrix DCM_B6, IvDCM_B6;
                        double vtmp;
                        //
                        matrixCopy(x.DCM_B,DCM_B6);
                        matrixInverse(DCM_B6, IvDCM_B6);

                        vectorSub(POS0,x.BetaPos,PosDBtmp);
                        vectorTransform(IvDCM_B6,PosDBtmp,PosDBtmp); 

                        matrixCopy(x1.DCM_B,DCM_B6);
                        
                        vectorTransform(DCM_B6,PosDBtmp,PosDBtmp);
                        vectorAdd(PosDBtmp,x1.BetaPos,POS1);
					}

                    else if (cl.status[j] == -2)
					{
						cl1.status[j] = cl.status[j]; //Fixed frames, Axayaza

                        int ind2, cfid2;
                        Matrix DCM_A5,DCM_B5,IvDCM_A5,IvDCM_B5;
                        Vector aa, ae, ar, ac;
                        Vector PosDBtmp, POSin1, nvec, Force;
                        
                        ind2 = 0;
                        cfid2 = 0;

                        double vtmp;
                        //
                        matrixCopy(x.DCM_A,DCM_B5); // DCM_B : Oxyz to OXYZ
                        matrixCopy(x.DCM_B,DCM_A5);
                        matrixInverse(DCM_B5, IvDCM_B5);

                        vectorSub(POS0,x.AlphaPos,PosDBtmp);
                        vectorTransform(IvDCM_B5,PosDBtmp,PosDBtmp); 

                        matrixCopy(x1.DCM_A,DCM_B5);
                        
                        vectorTransform(DCM_B5,PosDBtmp,PosDBtmp);
                        vectorAdd(PosDBtmp,x1.AlphaPos,POS1);

                        // %%%%%%%%%%%%%%%%

                        // if (j == 4819)
                        //     {
                        //         // cout << "x.BetaPos[0]:  " << x.BetaPos[0] << "   |  x.BetaPos[1]:   " << x.BetaPos[1] << "   |  x.BetaPos[2]:   " << x.BetaPos[2] << endl;
                        //         // cout << "x.AlphaPos[0]:  " << x.AlphaPos[0] << "   |  x.AlphaPos[1]:   " << x.AlphaPos[1] << "   |  x.AlphaPos[2]:   " << x.AlphaPos[2] << endl;
                        //         double xxxxx = 1;
                        //     }

                        getBinaryForce(Force,varpool3,x1,alpha,beta,POS1,d,j,kep,ss);

                        //绝对加速度

                        double tmpmass;
                        tmpmass = 1/d.Mass[j];

                        vectorScale(Force,tmpmass,aa);

                        // if (j == 0)
                        //     {
                        //         cout << "aa[0]:  " << aa[0] << "   |  aa[1]:   " << aa[1] << "   |  aa[2]:   " << aa[2] << endl;
                        //     }

                        //牵连加速度

                        if (j == 111)
                            {
                                double xxxxx = 1.1111;
                            }

                        Vector Omegatmp, Fc, Angvel_Bxyz, Angvel_BXYZ, Rtmp, acc_b, rho, RelVel; //Omegatmp: 角速度单位向量
                        double tmp1;

                        Angvel_Bxyz[0] = x.AlphaAngularmomentum[0]/alpha.InertVec[0]*(-1);
                        Angvel_Bxyz[1] = x.AlphaAngularmomentum[1]/alpha.InertVec[1]*(-1);
                        Angvel_Bxyz[2] = x.AlphaAngularmomentum[2]/alpha.InertVec[2]*(-1);


                        vectorSub(POS1,x1.AlphaPos,PosDBtmp);  
                        matrixCopy(x1.DCM_A,DCM_B5);
                        matrixInverse(DCM_B5, IvDCM_B5);
                        vectorTransform(IvDCM_B5,PosDBtmp,PosDBtmp); // 当前时刻颗粒在主星表面的位置
                        vectorCross(Angvel_Bxyz,PosDBtmp,VEL1);// 当前时刻颗粒在主星表面的速度， 固定坐标系下
                        vectorCopy(VEL1,RelVel);

                        
                        vectorTransform(DCM_B5,VEL1,VEL1); //当前时刻颗粒在空间中的速度


                        vectorTransform(DCM_B5,Angvel_Bxyz,Angvel_BXYZ);
                        vectorTransform(DCM_B5,RelVel,RelVel);

                        tmp1 = 1/alpha.TotalMass;

                        vectorScale(gt2.GravForceA,tmp1,acc_b);      

                        if (j == 0)
                            {
                                // cout << "acc_b[0]:  " << acc_b[0] << "   |  acc_b[1]:   " << acc_b[1] << "   |  acc_b[2]:   " << acc_b[2] << endl;
                                double xxxxxx = 1.1111;
                            }

                        // ae = acc_b + Qw × (Qw × rho) + [(Qw_f)w+Qdw × rho]

                        Vector tmpV5, tmpV6, tmpV7, tmpV8, tmpV9;

                        vectorSub(POS1,x.AlphaPos,rho);
                        
                        Matrix W, tmpM1;
                        Matrixset(W, 0, -Angvel_Bxyz[2], Angvel_Bxyz[1], Angvel_Bxyz[2], 0, -Angvel_Bxyz[0], -Angvel_Bxyz[1], Angvel_Bxyz[0], 0);

                        // 第一项 acc_B
                        // 第二项 

                        vectorCross(Angvel_BXYZ,rho,tmpV5);
                        vectorCross(Angvel_BXYZ,tmpV5,tmpV5);

                        //第三项

                        matrixMultiply(DCM_B5,W,tmpM1);
                        vectorTransform(tmpM1,Angvel_Bxyz,tmpV6);
                        vectorTransform(DCM_B5,domg_beta,tmpV7);
                        vectorAdd(tmpV6,tmpV7,tmpV8);
                        vectorCross(tmpV8,rho,tmpV8);

                        // 加和
                        vectorAdd(acc_b,tmpV5,tmpV9);
                        vectorAdd(tmpV9,tmpV8,ae);

                        //科氏加速度
                        vectorCross(Angvel_BXYZ,RelVel,ac);

                        //相对加速度
                        vectorSub(aa,ae,ar);
                        vectorSub(ar,ac,ar);

                        if (j == 0)
                            {
                                // cout << "ae[0]:  " << ae[0] << "   |  ae[1]:   " << ae[1] << "   |  ae[2]:   " << ae[2] << endl;
                                double xxxxxxxx = 1.1111;
                            }
                        

                        //相对加速度投影到从星上
                        vectorTransform(IvDCM_B5,ar,ar);



                        // 寻找当地外法向量
                        double tmp4,flag1;
                        Vector Nvec;

                        vectorSub(POS1, x1.AlphaPos, POSin1); // BXYZ
                        vectorTransform(IvDCM_B5, POSin1, POSin1); // BXbYbZb

                        // cout << "POSin[0]:  " << POSin[0] << "   |  POSin[1]:   " << POSin[1] << "   |  POSin[2]:   " << POSin[2] << endl;


                        CheckCollision(polya, POSin1, ind2, cfid2);

                        // if (i = 144000)
                        // {
                        //         cout << "cfid2 == " << cfid2 << "   |  ind2 ==   " << ind2 << endl;
                        // }

						

                        vectorSet(nvec, polya.FaceNormVec[cfid2][0], polya.FaceNormVec[cfid2][1], polya.FaceNormVec[cfid2][2]);


                        tmp4 = vectorMag(nvec);

                        tmp4 = 1/tmp4;

                        vectorScale(nvec,tmp4,Nvec);

                        // 相对加速度和当地法向量做点积

                        flag1 = vectorDot(ar,Nvec);

                        if (flag1 > 0)
                        {
                            cl1.status[j] = 2;
                            vectorScale(POS1,1.0025,POS1);
                            // vectorCopy(POS1,POS1);
                            // vectorCopy(VEL0, VEL1);

                        }




					}

                    else
					{
						cl1.status[j] = cl.status[j]; //Fixed frames, Axayaza or Bxbybzb
						vectorCopy(POS0, POS1);
						vectorCopy(VEL0, VEL1);
					}
					//
					for (k = 0; k < 3; k++)
					{
						cl1.DebrisPos[j][k] = POS1[k];
						cl1.DebrisVel[j][k] = VEL1[k];
					}

                      
					//
				}

                

            }

            //printf("Number of particles in the orbit: %d\n", orbitnum);

                # pragma omp for schedule (dynamic) //schedule (dynamic,1) shared(p,c,ss,cl1)//private ( ID )


				for (j = 0; j < d.NumDebris; j++)
				{
					int ColliInd = 0;
					Vector POS11, VEL11, POSout, VELout;

					if (cl1.status[j] == 2) //  status = 2 : just calculated in binary orbits
					{
						vectorSet(POS11, cl1.DebrisPos[j][0], cl1.DebrisPos[j][1], cl1.DebrisPos[j][2]);
						vectorSet(VEL11, cl1.DebrisVel[j][0], cl1.DebrisVel[j][1], cl1.DebrisVel[j][2]);

						if (vectorMag(POS11) > ss.SOI)
						{
							cl1.status[j] = 1; //  status = 1 : has entered the orbit of the solar system
							//cout<<POS1[0]<<",  "<<POS1[1]<<",  "<<POS1[2]<<endl;
							// OXYZ2SXYZ(kep1, POS1, VEL1);
							//cout<<POS1[0]<<",  "<<POS1[1]<<",  "<<POS1[2]<<endl;
							for (k = 0; k < 3; k++)
							{
								cl1.DebrisPos[j][k] = POS11[k];
								cl1.DebrisVel[j][k] = VEL11[k];
							}
						}
						else
						{
							DetectCollisionState(ColliInd, alpha, beta, polya, polyb, odeopt, POS11, VEL11, x1, POSout, VELout);
							if (ColliInd == 0)
							{
								//cl1.status[j] = ColliInd;
								cl1.status[j] = 0;
								for (k = 0; k < 3; k++)
								{
									cl1.DebrisPos[j][k] = POSout[k];
									cl1.DebrisVel[j][k] = VELout[k];
								}
							}
							else if (ColliInd == -2)
							{
                                cl1.status[j] = -2;
								for (k = 0; k < 3; k++)
								{
									cl1.DebrisPos[j][k] = POSout[k];
									cl1.DebrisVel[j][k] = VELout[k];
								}
							}
                            else 
                            {
                                cl1.status[j] = -1;
								for (k = 0; k < 3; k++)
								{
									cl1.DebrisPos[j][k] = POSout[k];
									cl1.DebrisVel[j][k] = VELout[k];
								}                           
                            }
						}
					}
				}
				//      
                }      
                gravCopy(gt1, gt0);
				phaseCopy(x1, x);
				KepCopy(kep1, kep);
				CloudCopy(d, cl1, cl);
                
                
				//
				count++;
				//
				//cout<<kep.Position[0]<<"  "<<kep.Position[1]<<"  "<<kep.Position[2]<<endl;
				//
				// output module
                //printf("count: %d\n", count);
            outputnum++;
            if (outputnum == 86400)
            {
                daynum++;
                printf("The program has run to : %d d/1825d\n", daynum);
                outputnum = 0;
            }

			if (count == odeopt.OutputInterval)
            {
                //output (x, kep, cl)
                Filename = "res.";
                convert<<setw(10)<<setfill('0')<<i+1;
                Filename.append(convert.str());
                Filename.append(".bt");
                //
                const char *outfile = Filename.c_str();
                //
                resfile.open(outfile);
                resfile<<setiosflags(ios::scientific)<<setprecision(PrecDouble); // output precision and format
                resfile<<setw(WidthDouble)<<kep.SemiMajorAxis<<setw(WidthDouble)<<kep.Eccentricity<<setw(WidthDouble)<<kep.LongAscendNode \
                <<setw(WidthDouble)<<kep.Inclination<<setw(WidthDouble)<<kep.ArgPeriapsis<<setw(WidthDouble)<<kep.MeanAnomaly<<endl;
                resfile<<setw(WidthInt)<<ID_Alpha<<setw(WidthDouble)<<x.AlphaPos[0]<<setw(WidthDouble)<<x.AlphaPos[1]<<setw(WidthDouble)<<x.AlphaPos[2] \
                <<setw(WidthDouble)<<x.AlphaVel[0]<<setw(WidthDouble)<<x.AlphaVel[1]<<setw(WidthDouble)<<x.AlphaVel[2] \
                <<setw(WidthDouble)<<x.DCM_A[0][0]<<setw(WidthDouble)<<x.DCM_A[0][1]<<setw(WidthDouble)<<x.DCM_A[0][2] \
                <<setw(WidthDouble)<<x.DCM_A[1][0]<<setw(WidthDouble)<<x.DCM_A[1][1]<<setw(WidthDouble)<<x.DCM_A[1][2] \
                <<setw(WidthDouble)<<x.DCM_A[2][0]<<setw(WidthDouble)<<x.DCM_A[2][1]<<setw(WidthDouble)<<x.DCM_A[2][2] \
                <<setw(WidthDouble)<<x.AlphaAngVel[0]<<setw(WidthDouble)<<x.AlphaAngVel[1]<<setw(WidthDouble)<<x.AlphaAngVel[2]<<endl;
                resfile<<setw(WidthInt)<<ID_Beta<<setw(WidthDouble)<<x.BetaPos[0]<<setw(WidthDouble)<<x.BetaPos[1]<<setw(WidthDouble)<<x.BetaPos[2] \
                <<setw(WidthDouble)<<x.BetaVel[0]<<setw(WidthDouble)<<x.BetaVel[1]<<setw(WidthDouble)<<x.BetaVel[2] \
                <<setw(WidthDouble)<<x.DCM_B[0][0]<<setw(WidthDouble)<<x.DCM_B[0][1]<<setw(WidthDouble)<<x.DCM_B[0][2] \
                <<setw(WidthDouble)<<x.DCM_B[1][0]<<setw(WidthDouble)<<x.DCM_B[1][1]<<setw(WidthDouble)<<x.DCM_B[1][2] \
                <<setw(WidthDouble)<<x.DCM_B[2][0]<<setw(WidthDouble)<<x.DCM_B[2][1]<<setw(WidthDouble)<<x.DCM_B[2][2] \
                <<setw(WidthDouble)<<x.BetaAngVel[0]<<setw(WidthDouble)<<x.BetaAngVel[1]<<setw(WidthDouble)<<x.BetaAngVel[2]<<endl;
                //
                for (int j=0; j<d.NumDebris; j++)
                {
                    resfile<<setw(WidthInt)<<d.Id[j]<<setw(WidthInt)<<cl.status[j]<<setw(WidthDouble)<<cl.DebrisPos[j][0]<<setw(WidthDouble)<<cl.DebrisPos[j][1] \
                    <<setw(WidthDouble)<<cl.DebrisPos[j][2]<<setw(WidthDouble)<<cl.DebrisVel[j][0]<<setw(WidthDouble)<<cl.DebrisVel[j][1]<<setw(WidthDouble) \
                    <<cl.DebrisVel[j][2]<<endl;
                }
                //output cl // resfile<<endl;
                resfile.close();
                //
                convert.str("");
                count = 0;
                //
            }

            for (ii = 0; ii < alpha.NodeNum; ii++) delete[] varpool3.alpha[ii];
            delete[] varpool3.alpha;

            for (ii = 0; ii < beta.NodeNum; ii++) delete[] varpool3.beta[ii];
            delete[] varpool3.beta;
			
			//
		}
		//
	}
}

void BinarySIStep(double t, BPHASE &x, BPHASE &x1, GRAV &gt0, GRAV &gt1, ODE_OPTION &odeopt, fem &alpha, fem &beta, KEPORB &kep, SOLAR &ss, MASCON &varpool)
{

	int i, j;
    Matrix DCM_A1, IvDCM_A1, DCM_B1, IvDCM_B1;

    Vector tmpV1, tmpV2;
    Vector AlphaAngularmomentumtmp, BetaAngularmomentumtmp,AlphaMomentumtmp,BetaMomentumtmp,Alphapostmp,Betapostmp;
    Vector AlphaAngularmomentumtmp1, BetaAngularmomentumtmp1;
    Matrix R1a,R2a,R3a,R4a,R5a,R1b,R2b,R3b,R4b,R5b;
    double tmp1,tmp2,tmp3,theta;

    matrixCopy(x.DCM_A, DCM_A1);
    matrixCopy(x.DCM_B, DCM_B1);

	matrixNorm(DCM_A1,DCM_A1);
    matrixNorm(DCM_B1,DCM_B1);

    matrixInverse(DCM_A1, IvDCM_A1);
    matrixInverse(DCM_B1, IvDCM_B1);

    tmp1 = 0.5*odeopt.StepSize;

    vectorScale(gt0.GravTorqueA,tmp1,AlphaAngularmomentumtmp);
    vectorAdd(x.AlphaAngularmomentum,AlphaAngularmomentumtmp,AlphaAngularmomentumtmp);
    vectorScale(gt0.GravTorqueB,tmp1,BetaAngularmomentumtmp);
    vectorAdd(x.BetaAngularmomentum,BetaAngularmomentumtmp,BetaAngularmomentumtmp);

    vectorScale(gt0.GravForceA,tmp1,AlphaMomentumtmp);
    vectorAdd(x.AlphaMomentum,AlphaMomentumtmp,AlphaMomentumtmp);
    vectorScale(gt0.GravForceB,tmp1,BetaMomentumtmp);
    vectorAdd(x.BetaMomentum,BetaMomentumtmp,BetaMomentumtmp);

    tmp2 = odeopt.StepSize/alpha.TotalMass;
    vectorScale(AlphaMomentumtmp,tmp2,Alphapostmp);
    vectorAdd(x.AlphaPos,Alphapostmp,Alphapostmp);
    vectorCopy(Alphapostmp,x1.AlphaPos);


    tmp3 = odeopt.StepSize/beta.TotalMass;
    vectorScale(BetaMomentumtmp,tmp3,Betapostmp);
    vectorAdd(x.BetaPos,Betapostmp,Betapostmp);
    vectorCopy(Betapostmp,x1.BetaPos);

    Vector DCM_A11, DCM_A12, DCM_A13;
    Vector DCM_B11, DCM_B12, DCM_B13;

    vectorSet(DCM_A11,DCM_A1[0][0],DCM_A1[0][1],DCM_A1[0][2]);
    vectorSet(DCM_A12,DCM_A1[1][0],DCM_A1[1][1],DCM_A1[1][2]);
    vectorSet(DCM_A13,DCM_A1[2][0],DCM_A1[2][1],DCM_A1[2][2]);

    vectorSet(DCM_B11,DCM_B1[0][0],DCM_B1[0][1],DCM_B1[0][2]);
    vectorSet(DCM_B12,DCM_B1[1][0],DCM_B1[1][1],DCM_B1[1][2]);
    vectorSet(DCM_B13,DCM_B1[2][0],DCM_B1[2][1],DCM_B1[2][2]);



    theta = 0.5*odeopt.StepSize*AlphaAngularmomentumtmp[0]/alpha.InertVec[0];
    RotationX(R1a,theta);   
    vectorTransform(R1a, AlphaAngularmomentumtmp, AlphaAngularmomentumtmp1);
    vectorTransform(R1a,DCM_A11,DCM_A11);
    vectorTransform(R1a,DCM_A12,DCM_A12);
    vectorTransform(R1a,DCM_A13,DCM_A13);



    theta = 0.5*odeopt.StepSize*AlphaAngularmomentumtmp1[1]/alpha.InertVec[1];
    RotationY(R2a,theta);
    vectorTransform(R2a, AlphaAngularmomentumtmp1, AlphaAngularmomentumtmp1);   
    vectorTransform(R2a,DCM_A11,DCM_A11);
    vectorTransform(R2a,DCM_A12,DCM_A12);
    vectorTransform(R2a,DCM_A13,DCM_A13);

    theta = odeopt.StepSize*AlphaAngularmomentumtmp1[2]/alpha.InertVec[2];
    RotationZ(R3a,theta);
    vectorTransform(R3a, AlphaAngularmomentumtmp1, AlphaAngularmomentumtmp1);
    vectorTransform(R3a,DCM_A11,DCM_A11);
    vectorTransform(R3a,DCM_A12,DCM_A12);
    vectorTransform(R3a,DCM_A13,DCM_A13);

    theta = 0.5*odeopt.StepSize*AlphaAngularmomentumtmp1[1]/alpha.InertVec[1];
    RotationY(R4a,theta);
    vectorTransform(R4a, AlphaAngularmomentumtmp1, AlphaAngularmomentumtmp1);
    vectorTransform(R4a,DCM_A11,DCM_A11);
    vectorTransform(R4a,DCM_A12,DCM_A12);
    vectorTransform(R4a,DCM_A13,DCM_A13);

    theta = 0.5*odeopt.StepSize*AlphaAngularmomentumtmp1[0]/alpha.InertVec[0];
    RotationX(R5a,theta);
    vectorTransform(R5a, AlphaAngularmomentumtmp1, AlphaAngularmomentumtmp1);
    vectorTransform(R5a,DCM_A11,DCM_A11);
    vectorTransform(R5a,DCM_A12,DCM_A12);
    vectorTransform(R5a,DCM_A13,DCM_A13);

	theta = 0.5 * odeopt.StepSize * BetaAngularmomentumtmp[0] / beta.InertVec[0];
	RotationX(R1b, theta);
	vectorTransform(R1b, BetaAngularmomentumtmp, BetaAngularmomentumtmp1);

    vectorTransform(R1b,DCM_B11,DCM_B11);
    vectorTransform(R1b,DCM_B12,DCM_B12);
    vectorTransform(R1b,DCM_B13,DCM_B13);

	theta = 0.5 * odeopt.StepSize * BetaAngularmomentumtmp1[1] / beta.InertVec[1];
	RotationY(R2b, theta);
	vectorTransform(R2b, BetaAngularmomentumtmp1, BetaAngularmomentumtmp1);
    
    vectorTransform(R2b,DCM_B11,DCM_B11);
    vectorTransform(R2b,DCM_B12,DCM_B12);
    vectorTransform(R2b,DCM_B13,DCM_B13);



	theta = odeopt.StepSize * BetaAngularmomentumtmp1[2] / beta.InertVec[2];
	RotationZ(R3b, theta);
	vectorTransform(R3b, BetaAngularmomentumtmp1, BetaAngularmomentumtmp1);
    
    vectorTransform(R3b,DCM_B11,DCM_B11);
    vectorTransform(R3b,DCM_B12,DCM_B12);
    vectorTransform(R3b,DCM_B13,DCM_B13);


	theta = 0.5 * odeopt.StepSize * BetaAngularmomentumtmp1[1] / beta.InertVec[1];
	RotationY(R4b, theta);
	vectorTransform(R4b, BetaAngularmomentumtmp1, BetaAngularmomentumtmp1);
    
    vectorTransform(R4b,DCM_B11,DCM_B11);
    vectorTransform(R4b,DCM_B12,DCM_B12);
    vectorTransform(R4b,DCM_B13,DCM_B13);



	theta = 0.5 * odeopt.StepSize * BetaAngularmomentumtmp1[0] / beta.InertVec[0];
	RotationX(R5b, theta);
	vectorTransform(R5b, BetaAngularmomentumtmp1, BetaAngularmomentumtmp1);
    
    vectorTransform(R5b,DCM_B11,DCM_B11);
    vectorTransform(R5b,DCM_B12,DCM_B12);
    vectorTransform(R5b,DCM_B13,DCM_B13);

    Matrixset(x1.DCM_A, DCM_A11[0], DCM_A11[1], DCM_A11[2], DCM_A12[0], DCM_A12[1], DCM_A12[2], DCM_A13[0], DCM_A13[1], DCM_A13[2]);
    Matrixset(x1.DCM_B, DCM_B11[0], DCM_B11[1], DCM_B11[2], DCM_B12[0], DCM_B12[1], DCM_B12[2], DCM_B13[0], DCM_B13[1], DCM_B13[2]);


    
    matrixInverse(x1.DCM_A, IvDCM_A1);
    matrixInverse(x1.DCM_B, IvDCM_B1);

    for (i=0; i<alpha.NodeNum; i++)
    {
        vectorSet(tmpV1,alpha.Nodes[i][0],alpha.Nodes[i][1],alpha.Nodes[i][2]);
        vectorTransform(x1.DCM_A, tmpV1, tmpV2);
        vectorAdd(x1.AlphaPos,tmpV2,tmpV1);
        for (j=0; j<3; j++)
        {
            varpool.alpha[i][j] = tmpV1[j];  // OXYZ
            varpool.alpha[i][j+3] = tmpV2[j]; // AXYZ
        }
    }
    //
    for (i=0; i<beta.NodeNum; i++)
    {
        vectorSet(tmpV1,beta.Nodes[i][0],beta.Nodes[i][1],beta.Nodes[i][2]);
        vectorTransform(x1.DCM_B, tmpV1, tmpV2);
        vectorAdd(x1.BetaPos,tmpV2,tmpV1);
        for (j=0; j<3; j++)
        {
            varpool.beta[i][j] = tmpV1[j]; // OXYZ
            varpool.beta[i][j+3] = tmpV2[j]; //BXYZ
        }
    }

    vectorZero(gt1.GravForceA); //
    vectorZero(gt1.GravForceB); //
    vectorZero(gt1.GravTorqueA); //
    vectorZero(gt1.GravTorqueB); //


    BinaryForceTorque(gt1.GravForceA, gt1.GravTorqueA, gt1.GravTorqueB, alpha, beta, varpool, kep, ss);

    vectorScale(gt1.GravForceA, -1.0, gt1.GravForceB); // attraction on Beta in OXYZ
    // torque on Beta transfered to AXaYaZa
    vectorTransform(IvDCM_A1, gt1.GravTorqueA, gt1.GravTorqueA);
    // torque on Beta transfered to BXbYbZb
    vectorTransform(IvDCM_B1, gt1.GravTorqueB, gt1.GravTorqueB);


    if (ss.SolarTide == 1)
        {
            Vector tmpV4, tmpV5, tmpV1, tmpV2, tmpV3;
            double tmpd1;
            // The attractions from the sun should be added, gravitational torques from the sun are neglectable
            //
            vectorAdd(kep.Position, x1.AlphaPos, tmpV1);    // Gravitational acceleration of the primary from the Sun, SXYZ / OXYZ
            tmpd1 = vectorMag(tmpV1);

            // cout << tmpd1 << "  |Distance|"  << endl;
            tmpd1 = -G*SM/(tmpd1*tmpd1*tmpd1);
            vectorScale(tmpV1, tmpd1, tmpV1);
            //
            vectorAdd(kep.Position, x1.BetaPos, tmpV2);     // Gravitational acceleration of the secondary from the Sun, SXYZ / OXYZ
            tmpd1 = vectorMag(tmpV2);
            tmpd1 = -G*SM/(tmpd1*tmpd1*tmpd1);
            vectorScale(tmpV2, tmpd1, tmpV2);
            //
            tmpd1 = vectorMag(kep.Position);               // Inertial acceleration, OXYZ
            tmpd1 = G*SM/(tmpd1*tmpd1*tmpd1);
            vectorScale(kep.Position, tmpd1, tmpV3);
            //
            vectorAdd(tmpV1, tmpV3, tmpV4);     // The solar tide Acc of the primary
            vectorAdd(tmpV2, tmpV3, tmpV5);     // The solar tide Acc of the secondary
            vectorScale(tmpV4,alpha.TotalMass,tmpV1);
            vectorScale(tmpV5,beta.TotalMass,tmpV2);

            
            // cout << tmpV1[0] << "  Sun" << tmpV1[1] << "  |" << tmpV1[2] << endl;
            // cout << gt1.GravForceA[0] << "  Gta" << gt1.GravForceA[1] << "  |" << gt1.GravForceA[2] << endl;

            vectorAdd(gt1.GravForceA,tmpV1,gt1.GravForceA);
            vectorAdd(gt1.GravForceB,tmpV2,gt1.GravForceB);

            // vectorAdd(tmpV4, AccA, tmpV1);
            // vectorAdd(tmpV5, AccB, tmpV2);
            // vectorCopy(tmpV1, AccA);
            // vectorCopy(tmpV2, AccB);
        }



    vectorScale(gt1.GravTorqueA, -1.0, gt1.GravTorqueA);
    vectorScale(gt1.GravTorqueB, -1.0, gt1.GravTorqueB);

    tmp1 = 0.5*odeopt.StepSize;

    vectorScale(gt1.GravTorqueA, tmp1,x1.AlphaAngularmomentum);
    vectorAdd(AlphaAngularmomentumtmp1,x1.AlphaAngularmomentum,x1.AlphaAngularmomentum);
    vectorScale(gt1.GravTorqueB, tmp1,x1.BetaAngularmomentum);
    vectorAdd(BetaAngularmomentumtmp1,x1.BetaAngularmomentum,x1.BetaAngularmomentum);

    vectorScale(gt1.GravForceA,tmp1,x1.AlphaMomentum);
    vectorAdd(AlphaMomentumtmp,x1.AlphaMomentum,x1.AlphaMomentum);
    vectorScale(gt1.GravForceB,tmp1,x1.BetaMomentum);
    vectorAdd(BetaMomentumtmp,x1.BetaMomentum,x1.BetaMomentum);

    tmp2 = 1/alpha.TotalMass;
    vectorScale(x1.AlphaMomentum,tmp2,x1.AlphaVel);
    tmp3 = 1/beta.TotalMass;
    vectorScale(x1.BetaMomentum,tmp3,x1.BetaVel);

    vectorDivision(x1.AlphaAngularmomentum,alpha.InertVec,x1.AlphaAngVel);
    vectorDivision(x1.BetaAngularmomentum,beta.InertVec,x1.BetaAngVel);
    
}

//
void OXYZ2SXYZ(KEPORB &kep, Vector Pos, Vector Vel)
{
    vectorAdd(kep.Position, Pos, Pos);
    vectorAdd(kep.Velocity, Vel, Vel);
}

//
void DebBinarySIStep(double t, Vector &POS0, Vector &VEL0, Vector &POS1, Vector &VEL1, ODE_OPTION &odeopt, BPHASE &x, fem &alpha, fem &beta, DEBRIS &d, int n, KEPORB &kep, SOLAR &ss, MASCON &varpool)
{
    
    int i, j;

    double tmp1, tmp2;
    Vector tmpV1, tmpV2, Pos, Vel, Momentum, Force, Momentumtmp;
    Vector Postmp;
    Matrix DCM_A1, DCM_B1;

    //POS: pos of debris (OXYZ)
    //VEL: vel of debris (OXYZ)

    vectorCopy(POS0,Pos);
    vectorCopy(VEL0,Vel);

    matrixCopy(x.DCM_A, DCM_A1);
    matrixCopy(x.DCM_B, DCM_B1);

    for (i=0; i<alpha.NodeNum; i++)
    {
        vectorSet(tmpV1,alpha.Nodes[i][0],alpha.Nodes[i][1],alpha.Nodes[i][2]);
        vectorTransform(DCM_A1, tmpV1, tmpV2);
        vectorAdd(x.AlphaPos,tmpV2,tmpV1);
        for (j=0; j<3; j++)
        {
            varpool.alpha[i][j] = tmpV1[j];  // OXYZ
            varpool.alpha[i][j+3] = tmpV2[j]; // AXYZ
        }
    }
    //
    for (i=0; i<beta.NodeNum; i++)
    {
        vectorSet(tmpV1,beta.Nodes[i][0],beta.Nodes[i][1],beta.Nodes[i][2]);
        vectorTransform(DCM_B1, tmpV1, tmpV2);
        vectorAdd(x.BetaPos,tmpV2,tmpV1);
        for (j=0; j<3; j++)
        {
            varpool.beta[i][j] = tmpV1[j]; // OXYZ
            varpool.beta[i][j+3] = tmpV2[j]; //BXYZ
        }
    }

    vectorZero(Force);
    getBinaryForce(Force,varpool,x,alpha,beta,Pos,d,n,kep,ss);

    vectorScale(Vel,d.Mass[n],Momentum);

    tmp1 = 0.5*odeopt.StepSize;
    vectorScale(Force,tmp1,Momentumtmp);
    vectorAdd(Momentum,Momentumtmp,Momentumtmp);

    tmp2 = odeopt.StepSize/d.Mass[n];
    vectorScale(Momentumtmp,tmp2,Postmp);
    vectorAdd(Pos,Postmp,Pos);

    vectorZero(Force);
    getBinaryForce(Force,varpool,x,alpha,beta,Pos,d,n,kep,ss);

    
    vectorScale(Force, tmp1, Momentum);
    vectorAdd(Momentumtmp, Momentum, Momentum);
    tmp1 = 1/d.Mass[n];
    vectorScale(Momentum, tmp1, Vel);

    vectorCopy(Pos,POS1);
    vectorCopy(Vel,VEL1);

}

//
void getBinaryForce(Vector Force, MASCON &varpool, BPHASE &x, fem &alpha, fem &beta, Vector Pos, DEBRIS &d, int n, KEPORB &kep, SOLAR &ss)
{
    int i, j;

    double tmpd, Coef, ddtest, GFAtest;
    Vector tmpV1, tmpV2, tmpV3;
    Vector GravForceA, GravForceB;
    Matrix DCM_A1, DCM_B1;

    // cout << "Apos"  << x.AlphaPos[0] << "  |" << x.AlphaPos[1] << "  |" << x.AlphaPos[2] << endl;
    // cout << "Bpos"  << x.BetaPos[0] << "  |" << x.BetaPos[1] << "  |" << x.BetaPos[2] << endl;
    // cout << "Dpos"  << Pos[0] << "  |" << Pos[1] << "  |" << Pos[2] << endl;

	vectorZero(GravForceA);
    for (i = 0; i < alpha.NodeNum; i++)
	{
		double dd, dd1;
		Vector VV1, VV2;

		vectorSet(VV2, Pos[0] - varpool.alpha[i][0], Pos[1] - varpool.alpha[i][1], \
			Pos[2] - varpool.alpha[i][2]);
		dd = vectorMag(VV2);
		dd1 = -1 * G * d.Mass[n] * alpha.NodeWeights[i] / (dd * dd * dd);
		vectorScale(VV2, dd1, VV1);
		// ddtest = vectorMag(VV1);
		if (dd < 2)
		{
			//cout << ddtest << endl;
			vectorZero(VV1);
		}
		
		vectorAdd(GravForceA, VV1, GravForceA);
		GFAtest = vectorMag(GravForceA);
	}

    // cout << "GravForceA"  << GravForceA[0] << "  |" << GravForceA[1] << "  |" << GravForceA[2] << endl;

	vectorZero(GravForceB);
	for (j = 0; j < beta.NodeNum; j++)
	{
		double dd, dd1;
		Vector VV3, VV4;

		vectorSet(VV4, Pos[0] - varpool.beta[j][0], Pos[1] - varpool.beta[j][1], \
			Pos[2] - varpool.beta[j][2]);
		dd = vectorMag(VV4);
		dd1 = -1 * G * d.Mass[n] * beta.NodeWeights[j] / (dd * dd * dd);
		vectorScale(VV4, dd1, VV3);

        if (dd < 2)
		{
			//cout << ddtest << endl;
			vectorZero(VV3);
		}

		vectorAdd(GravForceB, VV3, GravForceB);
	}

    // cout << "GravForceB"  << GravForceB[0] << "  |" << GravForceB[1] << "  |" << GravForceB[2] << endl;
	//cout << Pos[0] << "  |" << Pos[1] << "  |" << Pos[2] << endl;

	//cout << GravForceA[0] << "  |" << GravForceA[1] << "  |" << GravForceA[2] << "  |" << GravForceB[0] << "  |" << GravForceB[1] << "  |" << GravForceB[2] << endl;

    
    vectorAdd(GravForceA, GravForceB, Force);


    //  cout << "Binary"  << Force[0] << "  |" << Force[1] << "  |" << Force[2] << endl;
    //tmpd = G * d.Mass[n];
    //vectorScale(Force, tmpd, Force);
    //
    // Current position of the binary mass center OXYZ
    if (ss.SolarTide == 1 || ss.SolarPressure == 1)
    {
        vectorAdd(kep.Position, Pos, tmpV1);
        tmpd = vectorMag(tmpV1);
        tmpd = 1.0/(tmpd*tmpd*tmpd);
        vectorScale(tmpV1, tmpd, tmpV3);
        //
        if (ss.SolarTide == 1)
        {
            Coef = -G*SM*d.Mass[n];
            // Attraction from the Sun OXYZ
            vectorScale(tmpV3, Coef, tmpV1);
            // Inertial force of OXYZ
            tmpd = vectorMag(kep.Position);
            tmpd = -Coef/(tmpd*tmpd*tmpd);
            vectorScale(kep.Position, tmpd, tmpV2);
            // add to force
            vectorAdd(tmpV1, tmpV2, tmpV1);
            vectorAdd(Force, tmpV1, Force);

            //  cout << "SolarTide"  << tmpV1[0] << "  |" << tmpV1[1] << "  |" << tmpV1[2] << endl;
        }
        if (ss.SolarPressure == 1)
        {
            if (!chkOCC(kep.Position, x.AlphaPos, Pos, alpha.Ellipsoid, DCM_A1))
            {
                if (!chkOCC(kep.Position, x.BetaPos, Pos, beta.Ellipsoid, DCM_B1))
                {
                    // Solar radiation pressure from the Sun OXYZ
                    Coef = ss.Reflection*PI*d.Radius[n]*d.Radius[n]*SRF*AU*AU/LS;
                    vectorScale(tmpV3, Coef, tmpV1);
                    // add to force
                    vectorAdd(Force, tmpV1, Force);
                    //  cout << "SolarPressure"  << tmpV1[0] << "  |" << tmpV1[1] << "  |" << tmpV1[2] << endl;

					//  double xxxx;
					//  xxxx = 1;
                }
            }
        }
    }


}

//
void DetectCollisionState(int& ind, fem &alpha, fem &beta, poly &polya, poly &polyb, ODE_OPTION &odeopt, Vector POS1, Vector VEL1, BPHASE &x1, Vector POSout, Vector VELout)
{
    int ind1, ind2, cfid2, cfid1;
	ind1 = 0;
	ind2 = 0;
	cfid1 = 0;
	cfid2 = 0;
    Vector POSin, VELin, VELtmp, AngVel;
    Matrix DCM_A1, DCM_B1, IvDCM_A1, IvDCM_B1;
    double vtmp;
    //
    matrixCopy(x1.DCM_A,DCM_A1);
    matrixCopy(x1.DCM_B,DCM_B1);

    matrixInverse(DCM_A1, IvDCM_A1);
    matrixInverse(DCM_B1, IvDCM_B1);


    vectorSub(POS1, x1.BetaPos, POSin); // BXYZ
    vectorTransform(IvDCM_B1, POSin, POSin); // BXbYbZb
    CheckCollision(polyb, POSin, ind2, cfid2);
    //
    if (ind2)  
    {
        ind = ID_Beta;

        vectorSub(VEL1, x1.BetaVel, VELin); // BXYZ  Relative vel of B
        vectorTransform(IvDCM_B1, VELin, VELin); //BXbYbZb

        vectorDivision(x1.BetaAngularmomentum, beta.InertVec, AngVel);

        vectorCross(AngVel,POSin,VELtmp); //
        vectorSub(VELin,VELtmp,VELin);
        // vectorAdd(VELin, VELtmp, VELin);

        CollisionModel(polyb, odeopt, cfid2, POSin, VELin, POSout, VELout);

		vtmp = vectorMag(VELout);

		if (vtmp < 0.000002)
		{
            vectorAdd(VELout,VELtmp,VELout);
            // vectorSub(VELout, VELtmp, VELout);

            vectorTransform(DCM_B1, POSout, POSout); //BXYZ
            vectorTransform(DCM_B1, VELout, VELout); //BXYZ
            vectorAdd(POSout, x1.BetaPos, POSout); //OXYZ
            vectorAdd(VELout, x1.BetaVel, VELout); //OXYZ
			ind = -1;
			return;
		}

        vectorAdd(VELout,VELtmp,VELout);
        // vectorSub(VELout, VELtmp, VELout);

        vectorTransform(DCM_B1, POSout, POSout); //BXYZ
        vectorTransform(DCM_B1, VELout, VELout); //BXYZ
        vectorAdd(POSout, x1.BetaPos, POSout); //OXYZ
        vectorAdd(VELout, x1.BetaVel, VELout); //OXYZ
    }
    else
    {
        // DCM_A: the direction cosine matrix (transfer from AXYZ to AXaYaZa)
        vectorSub(POS1, x1.AlphaPos, POSin); // AXYZ
        vectorTransform(IvDCM_A1, POSin, POSin); //AXaYaZa
        CheckCollision(polya, POSin, ind1, cfid1);
        //
        if (ind1)
        {
            ind = ID_Alpha;
            vectorSub(VEL1, x1.AlphaVel, VELin); // AXYZ  Relative vel of A
            vectorTransform(IvDCM_A1, VELin, VELin); //AXaYaZa

            vectorDivision(x1.AlphaAngularmomentum, alpha.InertVec, AngVel);


			// cout << POSin[0] << "  |Posin|" << POSin[1] << "  |" << POSin[2] << endl;
			// cout << AngVel[0] << "  |AngVel|" << AngVel[1] << "  |" << AngVel[2] << endl;

            vectorCross(AngVel,POSin,VELtmp); //??????????????????????????????????
			//vectorCross(POSin, AngVel, VELtmp);

			//cout << VELtmp[0] << "  |1|" << VELtmp[1] << "  |" << VELtmp[2]  << endl;

            vectorSub(VELin,VELtmp,VELin);
			// vectorAdd(VELin, VELtmp, VELin);// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

			//cout << VELin[0] << "  |2|" << VELin[1] << "  |" << VELin[2] << endl;

            CollisionModel(polya, odeopt, cfid1, POSin, VELin, POSout, VELout);

			//cout << VELout[0] << "  |3|" << VELout[1] << "  |" << VELout[2] << endl;

            vtmp = vectorMag(VELout);

			if (vtmp < 0.00002)
			{
                vectorAdd(VELout,VELtmp,VELout);
                vectorTransform(DCM_A1, POSout, POSout); //BXYZ
                vectorTransform(DCM_A1, VELout, VELout); //BXYZ
                vectorAdd(POSout, x1.AlphaPos, POSout); //OXYZ
                vectorAdd(VELout, x1.AlphaVel, VELout); //OXYZ
				ind = -2;
				return;
			}
            vectorAdd(VELout,VELtmp,VELout);
			// vectorSub(VELout, VELtmp, VELout);//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

			// cout << VELout[0] << "  |4|" << VELout[1] << "  |" << VELout[2] << endl;



            vectorTransform(DCM_A1, POSout, POSout); //BXYZ
            vectorTransform(DCM_A1, VELout, VELout); //BXYZ
            vectorAdd(POSout, x1.AlphaPos, POSout); //OXYZ
            vectorAdd(VELout, x1.AlphaVel, VELout); //OXYZ
            //
        }
        else
        {
            ind = 0;
            vectorCopy(POS1, POSout);
            vectorCopy(VEL1, VELout);
        }
    }
}




void CollisionModel(poly &poly, ODE_OPTION &odeopt, int cfid2, Vector POSin, Vector VELin, Vector POSout, Vector VELout)
{
    int N;
    Vector nvec, vt_vec, vn_vec, vw_vec, VELtmp;
    double vni, vti, vns, vts, nvi, VNi, VTi, r, theta, ud, vd, wd, vnout, vtout, vwout, tmp;
    double mu, alphaN, alphaT, alphaN1, alphaT1, seed1, seed2, seed3, seed4;
    mu = 0.5;
    alphaN = 0.5;
    alphaT = 0.5;
    N = 999;

	//cout <<POSin[0] << "  |" << POSin[1] << "  |" << POSin[2] << "  |"<< VELin[0] << "  |" << VELin[1] << "  |" << VELin[2] << endl;

    vectorSet(nvec, poly.FaceNormVec[cfid2][0], poly.FaceNormVec[cfid2][1], poly.FaceNormVec[cfid2][2]);
    vectorNorm(nvec);
    vni = vectorDot(VELin,nvec);
    
    // vntmp = -1*vn;
    vectorScale(nvec,vni,vn_vec);//入射速度向量在平面法线上的投影
   

	vt_vec[0] = VELin[0] - vn_vec[0];
	vt_vec[1] = VELin[1] - vn_vec[1];
	vt_vec[2] = VELin[2] - vn_vec[2]; //入射速度向量在切线上的投影
	

	//cout << VELin[0] << "  |" << VELin[1] << "  |" << VELin[2] << "  |" << vt_vec[0] << "  |" << vt_vec[1] << "  |" << vt_vec[2] << endl;
    
    vti = vectorMag(vt_vec);

    
    vns = -alphaN*vni;
    vts = alphaT*vti;

	double vvtmp = sqrt(vns * vns + vts * vts);

	if (vvtmp < 0.0001)
	{
		if (poly.FaceCollisionFlag[cfid2] = 0)
		{
			vectorZero(VELout);
			vectorCopy(POSin, POSout);
			return;
		}
		
	}

	/*cout << vns << "  |" << vts <<"  |" << vvtmp << endl;*/

    // vni = -vni;

    alphaN1 = alphaN;
    alphaT1 = alphaT*(2-alphaT);


    nvi = sqrt(vni*vni + vti*vti);
    VNi = vni/nvi;
    VTi = vti/nvi;

    srand(time(NULL));

    seed1 = rand() % (N + 1) / (float)(N + 1);
    seed2 = rand() % (N + 1) / (float)(N + 1);
    seed3 = rand() % (N + 1) / (float)(N + 1);
    seed4 = rand() % (N + 1) / (float)(N + 1);

    r = sqrt(-alphaN1*log(seed1));
    theta = 2*PI*seed2;
    ud = sqrt(r*r+(1-alphaN1)*VNi*VNi+2*r*sqrt(1-alphaN1)*VNi*cos(theta))*nvi;
    r = sqrt(-alphaT1*log(seed3));
    theta = 2*PI*seed4;
    vd = (sqrt(1-alphaT1)*VTi+r*cos(theta))*nvi;
    wd = r*sin(theta)*nvi;

    vnout = mu*vns + (1-mu)*ud;
    vtout = mu*vts + (1-mu)*vd;
    vwout = (1-mu)*wd;

    vectorNorm(vn_vec);
	vectorScale(vn_vec, -1, vn_vec);
    vectorNorm(vt_vec);

    vectorCross(vt_vec,vn_vec,vw_vec);

	vectorNorm(vw_vec);

    vectorScale(vn_vec, vnout, vn_vec);
    vectorScale(vt_vec, vtout, vt_vec);
    vectorScale(vw_vec, vwout, vw_vec);

    vectorZero(VELout);
    vectorAdd(VELout,vn_vec,VELout);
    vectorAdd(VELout,vt_vec,VELout);
    vectorAdd(VELout,vw_vec,VELout);

    vectorScale(VELin, odeopt.StepSize, VELtmp);
    vectorSub(POSin, VELtmp, POSout);

	// cout << POSout[0] << "  |" << POSout[1] << "  |" << POSout[2] << "  |" << VELout[0] << "  |" << VELin[1] << "  |" << VELin[2] << endl;
    
}




//
void CheckCollision(poly &p, Vector r, int &ind, int &cfid)
{
    int tmp2;
    double tmp1,tmp3;
    Vector tpr_vec;
	double thetaN, phiN, flagg;
    int ti, pj, fid, flag;
    
    tmp1 = vectorDot(r,r);
    double xxx = 10000^4;

    if (tmp1>250000) 
    {
        ind = 0;
        cfid = 0;
        return;
    }
    else
    {
		//cout << r[0] << "  |" << r[1] << "  |" << r[2] << endl;
        xyz2tpr(r,tpr_vec);
        thetaN = tpr_vec[0]/p.DeltaAngle;
        phiN = tpr_vec[1]/p.DeltaAngle;

        ti = floor(thetaN);
        pj = floor(phiN);

        if (ti == 0)
        {
            fid = p.HeadIndex[pj];
        }
        else if(ti<p.MeshNum-1)
        {
            flagg = (thetaN - phiN) - (ti - pj);
            if (flagg>0)
            {
                fid = p.BodyLeftIndex[ti-1][pj];
            }
            else
            {
                fid = p.BodyRightIndex[ti-1][pj];
            }
        }
        else
        {
            fid = p.TailIndex[pj];
        }

        fid = fid - 1;

        Vector VV1,VV2;
        vectorSet(VV1, p.FaceNormVec[fid][0],p.FaceNormVec[fid][1], p.FaceNormVec[fid][2]);
        tmp2 = p.Faces[fid][0] - 1;
        vectorSet(VV2,p.Vertices[tmp2][0],p.Vertices[tmp2][1],p.Vertices[tmp2][2]);
        vectorSub(r,VV2,VV2);

        tmp3 = vectorDot(VV1,VV2);

        if (tmp3>= 0)
        {
            flag = 0;
            cfid = fid;
        }
        else
        {
            flag = 1;
            cfid = fid;
        }
        ind = flag;
    }
    
}

//

//
int chkOCC(Vector ROrb, Vector RBody, Vector Rpt, Vector Ellip, Matrix TranMX) //kep.Position, x.AlphaPos, Pos, p.Ellipsoid, DCM_A
{
    int Flag;
    double lambda, dis;
    Vector P0, tmpV, P01;
    
    vectorAdd(ROrb, RBody, P0); // P0: kep.Pos + x.alphaPos    alphaPos in SXYZ
    vectorScale(P0, -1.0, P0);
    
    vectorSub(Rpt, RBody, tmpV); // tmpV: deb��� alphaλ��
    
    vectorTransform(TranMX, P0, P0); // P0: sunPos in AXYZ
    vectorTransform(TranMX, tmpV, tmpV);// debPos in AXYZ
    
    vectorSub(tmpV, P0, P01); // P01: debPos - sunPos in AXYZ
    vectorSet(tmpV, 1.0/Ellip[0], 1.0/Ellip[1], 1.0/Ellip[2]);
    
    vectorMultiply(P0, tmpV, P0);
    vectorMultiply(P01, tmpV, P01);
    
    lambda = -vectorDot(P0,P01)/vectorMagSq(P01);
    if (lambda < 0.0) lambda = 0.0;
    if (lambda > 1.0) lambda = 1.0;
    
    vectorScale(P01, lambda, tmpV);
    vectorAdd(P0, tmpV, tmpV);
    dis = vectorMag(tmpV);
    
    if (dis<1.0) Flag = 1; // occulated
    else Flag = 0;
    
    return Flag;
}

//
void getHeliosForce(Vector Force, double t, Vector Pos, Vector Vel, fem &alpha, fem &beta, DEBRIS &d, int n, KEPORB &kep, SOLAR &ss)
{
    double BinaryMass, tmpd;
    Vector tmpV1, tmpV2;
    
    BinaryMass = alpha.TotalMass + beta.TotalMass;
    vectorScale(Pos, -1.0, tmpV1);
    vectorSub(kep.Position, Pos, tmpV2);
    
    tmpd = vectorMag(tmpV1);
    tmpd = G*SM*d.Mass[n]/(tmpd*tmpd*tmpd);
    vectorScale(tmpV1, tmpd, Force);
    tmpd = vectorMag(tmpV2);
    tmpd = G*BinaryMass*d.Mass[n]/(tmpd*tmpd*tmpd);
    vectorScale(tmpV2, tmpd, tmpV2);
    vectorAdd(Force, tmpV2, Force);
    
    if (ss.SolarPressure == 1)
    {
        tmpd = vectorMag(tmpV1);
        tmpd = -1.0/(tmpd*tmpd*tmpd);
        vectorScale(tmpV1, tmpd, tmpV2);
        //
        tmpd = ss.Reflection*PI*d.Radius[n]*d.Radius[n]*SRF*AU*AU/LS;
        vectorScale(tmpV2, tmpd, tmpV1);
        vectorAdd(Force, tmpV1, Force);
    }
    
    //
}

//
void DebHeliosSIStep(double t, Vector &POS0, Vector &VEL0, Vector &POS1, Vector &VEL1, ODE_OPTION &odeopt, fem &alpha, fem &beta, DEBRIS &d, int n, KEPORB &kep, SOLAR &ss)
{
        
    double tmp1, tmp2;
    Vector Pos, Vel, Momentum, Force, Momentumtmp;
    Vector Postmp;
    
    vectorCopy(POS0,Pos);
    vectorCopy(VEL0,Vel);

    vectorZero(Force);
    getHeliosForce(Force, t, Pos, Vel, alpha, beta, d, n, kep, ss);

    vectorScale(Vel,d.Mass[n],Momentum);

    tmp1 = 0.5*odeopt.StepSize;
    vectorScale(Force,tmp1,Momentumtmp);
    vectorAdd(Momentum,Momentumtmp,Momentumtmp);

    tmp2 = odeopt.StepSize/alpha.TotalMass;
    vectorScale(Momentumtmp,tmp2,Postmp);
    vectorAdd(Pos,Postmp,Pos);

    vectorZero(Force);
    getHeliosForce(Force, t, Pos, Vel, alpha, beta, d, n, kep, ss);

    vectorScale(Force, tmp1, Momentum);
    vectorAdd(Momentumtmp, Momentum, Momentum);
    tmp1 = 1/d.Mass[n];
    vectorScale(Momentum, tmp1, Vel);

    vectorCopy(Pos,POS1);
    vectorCopy(Vel,VEL1);
}

//
void InitialChecking(DEBRIS &d, BPHASE &x, CLOUD &cl, fem &alpha, fem &beta, poly &polya, poly &polyb, ODE_OPTION &odeopt)
{
    int i, ind;
	ind = 0;
    Vector Pos, Vel, POSout, VELout;
    //Matrix DCM_A, IvDCM_B;
    
    for (i=0; i<d.NumDebris; i++)
    {

        if (cl.status[i]==0)
        {
            vectorSet(Pos, cl.DebrisPos[i][0], cl.DebrisPos[i][1], cl.DebrisPos[i][2]); // OXYZ
            vectorSet(Vel, cl.DebrisVel[i][0], cl.DebrisVel[i][1], cl.DebrisVel[i][2]); // OXYZ

            DetectCollisionState(ind, alpha, beta, polya, polyb, odeopt, Pos, Vel, x, POSout, VELout);
            if (ind == ID_Alpha)
            {
                cout<<"Incorrect initialization: the on-orbit particle of No."<<d.Id[i]<<" overlapping with the primary!"<<endl;
                exit(1);
            }
            if (ind == ID_Beta)
            {
                cout<<"Incorrect initialization: the on-orbit particle of No."<<d.Id[i]<<" overlapping with the secondary!"<<endl;
                exit(1);
            }
        }
        //
    }
    
    cout<<"Initial checking correct! "<<endl;
}


void phaseCopy(BPHASE &x1, BPHASE &x2)
{
    vectorCopy(x1.AlphaPos, x2.AlphaPos);
    vectorCopy(x1.AlphaVel, x2.AlphaVel);
    vectorCopy(x1.AlphaMomentum, x2.AlphaMomentum);
    vectorCopy(x1.AlphaAngVel, x2.AlphaAngVel);
    matrixCopy(x1.DCM_A, x2.DCM_A);
    vectorCopy(x1.AlphaAngularmomentum,x2.AlphaAngularmomentum);
    //
    vectorCopy(x1.BetaPos, x2.BetaPos);
    vectorCopy(x1.BetaVel, x2.BetaVel);
    vectorCopy(x1.BetaMomentum, x2.BetaMomentum);
    vectorCopy(x1.BetaAngVel, x2.BetaAngVel);
    matrixCopy(x1.DCM_B, x2.DCM_B);
    vectorCopy(x1.BetaAngularmomentum,x2.BetaAngularmomentum);
}

void gravCopy(GRAV &gt1, GRAV &gt0)
{
    vectorCopy(gt1.GravForceA, gt0.GravForceA);
    vectorCopy(gt1.GravForceB, gt0.GravForceB);
    vectorCopy(gt1.GravTorqueA, gt0.GravTorqueA);
    vectorCopy(gt1.GravTorqueB, gt0.GravTorqueB);
}

void phaseExtract(BPHASE &x, double* y)
{
    y[0] = x.AlphaPos[0];
    y[1] = x.AlphaPos[1];
    y[2] = x.AlphaPos[2];
    y[3] = x.AlphaVel[0];
    y[4] = x.AlphaVel[1];
    y[5] = x.AlphaVel[2];
    y[6] = x.AlphaMomentum[0];
    y[7] = x.AlphaMomentum[1];
    y[8] = x.AlphaMomentum[2];
    y[9] = x.DCM_A[0][0];
    y[10] = x.DCM_A[0][1];
    y[11] = x.DCM_A[0][2];
    y[12] = x.DCM_A[1][0];
    y[13] = x.DCM_A[1][1];
    y[14] = x.DCM_A[1][2];
    y[15] = x.DCM_A[2][0];
    y[16] = x.DCM_A[2][1];
    y[17] = x.DCM_A[2][2];
    y[18] = x.AlphaAngVel[0];
    y[19] = x.AlphaAngVel[1];
    y[20] = x.AlphaAngVel[2];
    y[21] = x.AlphaAngularmomentum[0];
    y[22] = x.AlphaAngularmomentum[1];
    y[23] = x.AlphaAngularmomentum[2];
//
    y[24] = x.BetaPos[0];
    y[25] = x.BetaPos[1];
    y[26] = x.BetaPos[2];
    y[27] = x.BetaVel[0];
    y[28] = x.BetaVel[1];
    y[29] = x.BetaVel[2];
    y[30] = x.BetaMomentum[0];
    y[31] = x.BetaMomentum[1];
    y[32] = x.BetaMomentum[2];
    y[33] = x.DCM_B[0][0];
    y[34] = x.DCM_B[0][1];
    y[35] = x.DCM_B[0][2];
    y[36] = x.DCM_B[1][0];
    y[37] = x.DCM_B[1][1];
    y[38] = x.DCM_B[1][2];
    y[39] = x.DCM_B[2][0];
    y[40] = x.DCM_B[2][1];
    y[41] = x.DCM_B[2][2];
    y[42] = x.BetaAngVel[0];
    y[43] = x.BetaAngVel[1];
    y[44] = x.BetaAngVel[2];
    y[45] = x.BetaAngularmomentum[0];
    y[46] = x.BetaAngularmomentum[1];
    y[47] = x.BetaAngularmomentum[2];
}

void phaseCompress(const double* y, BPHASE &x)
{
    x.AlphaPos[0] = y[0];
    x.AlphaPos[1] = y[1];
    x.AlphaPos[2] = y[2];
    x.AlphaVel[0] = y[3];
    x.AlphaVel[1] = y[4];
    x.AlphaVel[2] = y[5];
    x.AlphaMomentum[0] = y[6];
    x.AlphaMomentum[1] = y[7];
    x.AlphaMomentum[2] = y[8];
    x.DCM_A[0][0] = y[9];
    x.DCM_A[0][1] = y[10];
    x.DCM_A[0][2] = y[11];
    x.DCM_A[1][0] = y[12];
    x.DCM_A[1][1] = y[13];
    x.DCM_A[1][2] = y[14];
    x.DCM_A[2][0] = y[15];
    x.DCM_A[2][1] = y[16];
    x.DCM_A[2][2] = y[17];
    x.AlphaAngVel[0] = y[18];
    x.AlphaAngVel[1] = y[19];
    x.AlphaAngVel[2] = y[20];
    x.AlphaAngularmomentum[0] = y[21];
    x.AlphaAngularmomentum[1] = y[22];
    x.AlphaAngularmomentum[2] = y[23];
//
    x.BetaPos[0] = y[24];
    x.BetaPos[1] = y[25];
    x.BetaPos[2] = y[26];
    x.BetaVel[0] = y[27];
    x.BetaVel[1] = y[28];
    x.BetaVel[2] = y[29];
    x.BetaMomentum[0] = y[30];
    x.BetaMomentum[1] = y[31];
    x.BetaMomentum[2] = y[32];
    x.DCM_B[0][0] = y[33];
    x.DCM_B[0][1] = y[34];
    x.DCM_B[0][2] = y[35];
    x.DCM_B[1][0] = y[36];
    x.DCM_B[1][1] = y[37];
    x.DCM_B[1][2] = y[38];
    x.DCM_B[2][0] = y[39];
    x.DCM_B[2][1] = y[40];
    x.DCM_B[2][2] = y[41];
    x.BetaAngVel[0] = y[42];
    x.BetaAngVel[1] = y[43];
    x.BetaAngVel[2] = y[44];
    x.BetaAngularmomentum[0] = y[45];
    x.BetaAngularmomentum[1] = y[46];
    x.BetaAngularmomentum[2] = y[47];
}

// void RotationX(Matrix &m, const double theta)
// {

// 	m[0][0] = 1;
// 	m[0][1] = 0;
// 	m[0][2] = 0;
// 	m[1][0] = 0;
// 	m[1][1] = (1-theta*theta/4)/(1+theta*theta/4);
// 	m[1][2] = -1*theta/(1+theta*theta/4);
// 	m[2][0] = 0;
// 	m[2][1] = theta/(1+theta*theta/4);
// 	m[2][2] = (1-theta*theta/4)/(1+theta*theta/4);

// }

// void RotationY(Matrix &m, const double theta)
// {

// 	m[0][0] = (1-theta*theta/4)/(1+theta*theta/4);
// 	m[0][1] = 0;
// 	m[0][2] = theta/(1+theta*theta/4);
// 	m[1][0] = 0;
// 	m[1][1] = 1;
// 	m[1][2] = 0;
// 	m[2][0] = -1*theta/(1+theta*theta/4);
// 	m[2][1] = 0;
// 	m[2][2] = (1-theta*theta/4)/(1+theta*theta/4);

// }

// void RotationZ(Matrix &m, const double theta)
// {

// 	m[0][0] = (1-theta*theta/4)/(1+theta*theta/4);
// 	m[0][1] = -1*theta/(1+theta*theta/4);
// 	m[0][2] = 0;
// 	m[1][0] = theta/(1+theta*theta/4);
// 	m[1][1] = (1-theta*theta/4)/(1+theta*theta/4);
// 	m[1][2] = 0;
// 	m[2][0] = 0;
// 	m[2][1] = 0;
// 	m[2][2] = 1;

// }

void RotationX(Matrix &m, const double theta)
{

	m[0][0] = 1;
	m[0][1] = 0;
	m[0][2] = 0;
	m[1][0] = 0;
	m[1][1] = cos(theta);
	m[1][2] = -1*sin(theta);
	m[2][0] = 0;
	m[2][1] = sin(theta);
	m[2][2] = cos(theta);

}

void RotationY(Matrix &m, const double theta)
{

	m[0][0] = cos(theta);
	m[0][1] = 0;
	m[0][2] = sin(theta);
	m[1][0] = 0;
	m[1][1] = 1;
	m[1][2] = 0;
	m[2][0] = -1*sin(theta);
	m[2][1] = 0;
	m[2][2] = cos(theta);

}

void RotationZ(Matrix &m, const double theta)
{

	m[0][0] = cos(theta);
	m[0][1] = -1*sin(theta);
	m[0][2] = 0;
	m[1][0] = sin(theta);
	m[1][1] = cos(theta);
	m[1][2] = 0;
	m[2][0] = 0;
	m[2][1] = 0;
	m[2][2] = 1;

}


