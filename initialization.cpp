//
//  initialization.cpp
//  
//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
//

#include "initialization.h"
#include "body.h"

void Initialization(char* parafile, fem &alpha, fem &beta, poly &polya, poly &polyb,  DEBRIS &d, ODE_OPTION &odeopt, SOLAR &ss, KEPORB &kep, BPHASE &x, CLOUD &cl)
{
    //
    double massratio;
    string alphafile, betafile, polyfileA1, polyfileB1, debrisfile, iniconfile;
    
    LoadParas(parafile, alphafile,  betafile, polyfileA1, polyfileB1, debrisfile, iniconfile, odeopt, ss);
    
	LoadFEM(alphafile, alpha);

	LoadFEM(betafile, beta);

    LoadPoly(polyfileA1, polya);

    LoadPoly(polyfileB1, polyb);
    
    LoadDebris(debrisfile, d);
    
    if (d.NumDebris==0)
    {
        cout<<"No debris gonna be traced! Program exits."<<endl;
        exit(1);
    }
    
    LoadInicon(iniconfile, d, kep, x, cl);
    //
    massratio = (alpha.TotalMass + beta.TotalMass)/SM;
    ss.SOI = kep.SemiMajorAxis*pow(massratio,0.4);
    //
    KEP2RV(kep);
    //
    // if (odeopt.StartStep == 0) // if this is the original step
    // {
    //     OriginCorrection(alpha, beta, x);
    //     InitialChecking(d, x, cl, alpha, beta, polya, polyb, odeopt);       
    // }
    
    //check debris num ~=0
    //check debris in binary ////collision & time horizon
    //check
    
    
}

//

void LoadParas(char* parafile, string &polyfileA, string &polyfileB, string &polyfileA1, string &polyfileB1, string &debrisfile, string &iniconfile, ODE_OPTION &odeopt, SOLAR &ss)
{
    int numchar, id, n;
    int flag;

    ifstream infile;
    string line, core, name, value;
    
    flag = 0;

    infile.open(parafile);
    
    while(getline(infile,line))
    {
        numchar = line.length();
        for (n=0; n<numchar; n++) if (line[n] == '#') break;
        core.assign(line,0,n);
        //
        core.erase(std::remove(core.begin(), core.end(), '\t'), core.end());
        core.erase(std::remove(core.begin(), core.end(), ' '), core.end());
        //
        if (! core.empty())
        {
            //
            numchar = core.length();
            id = core.find('=');
            assert((id>0) && (id+1<numchar));
            name.assign(core,0,id);
            n = numchar - id - 1;
            assert(n>0);
            value.assign(core,id+1,n);
            //
            if (name == "PolyhedronAFile")       {polyfileA = value;                              flag++;}
            else if (name == "PolyhedronBFile")     { polyfileB = value;                              flag++;}
            else if (name == "PolyA")     { polyfileA1 = value;                              flag++;}
            else if (name == "PolyB")     { polyfileB1 = value;                              flag++;}
            else if (name == "DebrisFile")      {debrisfile = value;                            flag++;}
            else if (name == "IniConFile")      {iniconfile = value;                            flag++;}
            //
            else if (name == "FunOption")       {odeopt.FunOption = atof(value.c_str());        flag++;}
            else if (name == "StepSize")        {odeopt.StepSize = atof(value.c_str());         flag++;}
            else if (name == "EndStep")         {odeopt.EndStep = atof(value.c_str());         flag++;}
            else if (name == "StartStep")       {odeopt.StartStep = atof(value.c_str());        flag++;}
            else if (name == "OutputInterval")  {odeopt.OutputInterval = atof(value.c_str());   flag++;}
            //
            else if (name == "SolarTide")       {ss.SolarTide = atof(value.c_str());            flag++;}
            else if (name == "SolarPressure")   {ss.SolarPressure = atof(value.c_str());        flag++;}
            else if (name == "ReflectionRate")  {ss.Reflection = atof(value.c_str());           flag++;}
            //
            else {cout<< "Unrecoganized item "<<name<<"!"<<endl; exit(1);}
            //
        }
    }
    infile.close();
    
    if (flag != NumParameter)
    {
        cout<<flag<<endl;
        cout<<"Important parameter setting missing!"<<endl;
        exit(1);
    }
}

void LoadInicon(string iniconfile, DEBRIS &d, KEPORB &kep, BPHASE &x, CLOUD &cl)
{
    int i, tmpi;
    const char* file;
    FILE * infile;
    
    CloudAlloc(d, cl);
    
    file = iniconfile.c_str();
    
    infile = fopen(file, "r");
    //
    fscanf(infile,"%lf%lf%lf%lf%lf%lf",&kep.SemiMajorAxis, &kep.Eccentricity, &kep.LongAscendNode, &kep.Inclination, &kep.ArgPeriapsis, &kep.MeanAnomaly);
    
    fscanf(infile,"%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&tmpi, &x.AlphaPos[0], &x.AlphaPos[1], &x.AlphaPos[2], &x.AlphaVel[0], \
                                                              &x.AlphaVel[1], &x.AlphaVel[2], &x.AlphaMomentum[0], &x.AlphaMomentum[1], &x.AlphaMomentum[2], \
                                                              &x.DCM_A[0][0], &x.DCM_A[0][1], &x.DCM_A[0][2], &x.DCM_A[1][0], \
                                                              &x.DCM_A[1][1], &x.DCM_A[1][2], &x.DCM_A[2][0], &x.DCM_A[2][1], \
                                                              &x.DCM_A[2][2], &x.AlphaAngVel[0], &x.AlphaAngVel[1], &x.AlphaAngVel[2], \
                                                              &x.AlphaAngularmomentum[0], &x.AlphaAngularmomentum[1], &x.AlphaAngularmomentum[2]);
    
    fscanf(infile,"%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&tmpi, &x.BetaPos[0], &x.BetaPos[1], &x.BetaPos[2], &x.BetaVel[0], \
                                                              &x.BetaVel[1], &x.BetaVel[2], &x.BetaMomentum[0], &x.BetaMomentum[1], &x.BetaMomentum[2], \
                                                              &x.DCM_B[0][0], &x.DCM_B[0][1], &x.DCM_B[0][2], &x.DCM_B[1][0], \
                                                              &x.DCM_B[1][1], &x.DCM_B[1][2], &x.DCM_B[2][0], &x.DCM_B[2][1], \
                                                              &x.DCM_B[2][2], &x.BetaAngVel[0], &x.BetaAngVel[1], &x.BetaAngVel[2], \
                                                              &x.BetaAngularmomentum[0], &x.BetaAngularmomentum[1], &x.BetaAngularmomentum[2]);
    
    for(i=0; i<d.NumDebris; i++) fscanf (infile,"%d%d%lf%lf%lf%lf%lf%lf",&tmpi, &cl.status[i], &cl.DebrisPos[i][0], &cl.DebrisPos[i][1], \
                                                                         &cl.DebrisPos[i][2], &cl.DebrisVel[i][0], &cl.DebrisVel[i][1], &cl.DebrisVel[i][2]);
    
    fclose (infile);
   
    // QuatNorm(x.AlphaOrien);
    // QuatNorm(x.BetaOrien);
}

void OriginCorrection(fem &alpha, fem &beta, BPHASE &x)
{
    double tmpd;
    Vector tmpV1, tmpV2, tmpV3;
    
    vectorScale(x.AlphaPos, alpha.TotalMass, tmpV1);
    vectorScale(x.BetaPos, beta.TotalMass, tmpV2);
    vectorAdd(tmpV1, tmpV2, tmpV3);
    //
    tmpd = 1.0/(alpha.TotalMass + beta.TotalMass);
    vectorScale(tmpV3, tmpd, tmpV3);
    //
    vectorSub(x.AlphaPos, tmpV3, x.AlphaPos);
    vectorSub(x.BetaPos, tmpV3, x.BetaPos);
    //
    vectorScale(x.AlphaVel, alpha.TotalMass, tmpV1);
    vectorScale(x.BetaVel, beta.TotalMass, tmpV2);
    vectorAdd(tmpV1, tmpV2, tmpV3);
    //
    vectorScale(tmpV3, tmpd, tmpV3);
    //
    vectorSub(x.AlphaVel, tmpV3, x.AlphaVel);
    vectorSub(x.BetaVel, tmpV3, x.BetaVel);
}


