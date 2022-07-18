#include "body.h"

void LoadFEM(string filename, FEM& fe)
{
	int i, j;
	double tmpd, tmpp;

	const char* file;
	FILE* infile;

	file = filename.c_str();
	//
	infile = fopen(file, "r");
	//
	fscanf(infile, "%lf", &fe.TotalMass);
	fscanf(infile, "%lf%lf%lf", &fe.MassCenter[0], &fe.MassCenter[1], &fe.MassCenter[2]);
	for (i = 0; i < 3; i++) fscanf(infile, "%lf%lf%lf", &fe.Inertia[i][0], &fe.Inertia[i][1], &fe.Inertia[i][2]);
	//
	fscanf(infile, "%d%d%d%d", &fe.NumVerts, &fe.NumFaces, &fe.NodeNum, &fe.ElemNum);
	//
	AllocFEM(fe);
	//
	for (i = 0; i < fe.NumVerts; i++) fscanf(infile, "%lf%lf%lf", &fe.Vertices[i][0], &fe.Vertices[i][1], &fe.Vertices[i][2]);
	//
	for (i = 0; i < fe.NumFaces; i++) fscanf(infile, "%d%d%d", &fe.Faces[i][0], &fe.Faces[i][1], &fe.Faces[i][2]);
	//
	for (i = 0; i < fe.NodeNum; i++) fscanf(infile, "%lf%lf%lf%lf%lf", &fe.NodeDensities[i], &fe.NodeWeights[i], &fe.Nodes[i][0], &fe.Nodes[i][1], &fe.Nodes[i][2]);
	//
	for (i = 0; i < fe.ElemNum; i++) fscanf(infile, "%d%d%d%d%lf", &fe.Elements[i][0], &fe.Elements[i][1], &fe.Elements[i][2], &fe.Elements[i][3], &fe.ElemJS[i]);
	//
	fclose(infile);
	//
	for (i = 0; i < 3; i++) fe.InertVec[i] = fe.Inertia[i][i];

	for (i = 0; i < 3; i++)
	{
		tmpd = 0.0;
		for (j = 0; j < fe.NumVerts; j++)
		{
			tmpp = fabs(fe.Vertices[j][i]);
			if (tmpp > tmpd) tmpd = tmpp;
		}
		fe.Ellipsoid[i] = tmpd;
	}

	//
}

void AllocFEM(FEM& fe)
{
	int i;

	fe.Vertices = new double* [fe.NumVerts];
	for (i = 0; i < fe.NumVerts; i++) fe.Vertices[i] = new double[3];

	fe.Faces = new int* [fe.NumFaces];
	for (i = 0; i < fe.NumFaces; i++) fe.Faces[i] = new int[3];

	fe.NodeDensities = new double[fe.NodeNum];
	fe.NodeWeights = new double[fe.NodeNum];

	fe.Nodes = new double* [fe.NodeNum];
	for (i = 0; i < fe.NodeNum; i++) fe.Nodes[i] = new double[3];

	fe.Elements = new int* [fe.ElemNum];
	for (i = 0; i < fe.ElemNum; i++) fe.Elements[i] = new int[4];

	fe.ElemJS = new double[fe.ElemNum];

	//cout<<"8"<<endl;

}


void LoadPoly(string filename, POLY& poly)
{
	int i, j;

	const char* file;
	FILE* infile;

	file = filename.c_str();
	//
	infile = fopen(file, "r");
	//
	fscanf(infile, "%lf%lf%d%d%d%d%d%d", &poly.DeltaAngle, &poly.EnvelopeRadiusSqr, &poly.MeshNum, &poly.VerticeNum, &poly.FaceNum, &poly.HeadIndexNum, &poly.BodyLeftIndexRowNum, &poly.BodyLeftIndexColumnNum);
	//
	PolyAlloc(poly);
	//
	for (i = 0; i < poly.VerticeNum; i++) fscanf(infile, "%lf%lf%lf", &poly.Vertices[i][0], &poly.Vertices[i][1], &poly.Vertices[i][2]);
	//
	for (i = 0; i < poly.FaceNum; i++) fscanf(infile, "%d%d%d", &poly.Faces[i][0], &poly.Faces[i][1], &poly.Faces[i][2]);
	//
	for (i = 0; i < poly.HeadIndexNum; i++) fscanf(infile, "%d", &poly.HeadIndex[i]);
	//
	for (i = 0; i < poly.HeadIndexNum; i++) fscanf(infile, "%d", &poly.TailIndex[i]);
	//
	for (i = 0; i < poly.BodyLeftIndexRowNum; i++)
	{
		for (j = 0; j < poly.BodyLeftIndexColumnNum; j++)
		{
			fscanf(infile, "%d", &poly.BodyLeftIndex[i][j]);
			//cout << poly.BodyLeftIndex[i][j] << endl;
		}
	}
	//
	for (i = 0; i < poly.BodyLeftIndexRowNum; i++)
	{
		for (j = 0; j < poly.BodyLeftIndexColumnNum; j++)
		{
			fscanf(infile, "%d", &poly.BodyRightIndex[i][j]);
		}
	}
	//
	for (i = 0; i < poly.FaceNum; i++) fscanf(infile, "%lf%lf%lf", &poly.FaceNormVec[i][0], &poly.FaceNormVec[i][1], &poly.FaceNormVec[i][2]);
	//
	for (i = 0; i < poly.FaceNum; i++) fscanf(infile, "%d", &poly.FaceCollisionFlag[i]);
	//
}

void PolyAlloc(POLY& poly)
{
	int i;

	poly.Vertices = new double* [poly.VerticeNum];
	for (i = 0; i < poly.VerticeNum; i++) poly.Vertices[i] = new double[3];

	poly.Faces = new int* [poly.FaceNum];
	for (i = 0; i < poly.FaceNum; i++) poly.Faces[i] = new int[3];

	poly.HeadIndex = new int[poly.HeadIndexNum];
	poly.TailIndex = new int[poly.HeadIndexNum];

	poly.BodyLeftIndex = new int* [poly.BodyLeftIndexRowNum];
	for (i = 0; i < poly.BodyLeftIndexRowNum; i++) poly.BodyLeftIndex[i] = new int[poly.BodyLeftIndexColumnNum];

	poly.BodyRightIndex = new int* [poly.BodyLeftIndexRowNum];
	for (i = 0; i < poly.BodyLeftIndexRowNum; i++) poly.BodyRightIndex[i] = new int[poly.BodyLeftIndexColumnNum];

	poly.FaceNormVec = new double* [poly.FaceNum];
	for (i = 0; i < poly.FaceNum; i++) poly.FaceNormVec[i] = new double[3];

	poly.FaceCollisionFlag = new int[poly.FaceNum];

	//cout<<"8"<<endl;

}



void LoadDebris(string debrisfile, DEBRIS &d)
{
    int i;
    const char* file;
    
    FILE * infile;
    
    file = debrisfile.c_str();
    
    infile = fopen(file, "r");
    
    fscanf(infile, "%d", &d.NumDebris);
    
    DebrisAlloc(d);
  
    for(i=0; i<d.NumDebris; i++) fscanf(infile,"%d%lf%lf",&d.Id[i],&d.Mass[i],&d.Radius[i]);
    
    fclose (infile);
}

void DebrisAlloc(DEBRIS &d)
{
    d.Id = new int[d.NumDebris];
    d.Mass = new double[d.NumDebris];
    d.Radius = new double[d.NumDebris];
}




