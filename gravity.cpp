//
//  gravity.cpp
//  
//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
// Calculat the attraction and potential on a united mass point from the polyhedron, just the geometric component (density and gravitational constant are not included in these subroutines)
//

#include "gravity.h"

//void GravAttraction(fem &alpha, Vector r, Vector F)
//{
//    int i, j;
//    double tmp1, tmp2, omegaf, le;
//    Vector tmpV1, tmpV2, FF, EE;
//    Vector tmpVS[3];
//    
//    double **RAY;
//    
//    RAY = new double*[alpha.NumVerts];
//    for (i=0; i<p.NumVerts; i++) RAY[i] = new double[4];
//    
//    for (i=0; i<p.NumVerts; i++)
//    {
//        vectorSet(tmpV1, p.Vertices[i][0], p.Vertices[i][1], p.Vertices[i][2]);
//        vectorSub(tmpV1, r, tmpV2);
//        tmp1 = vectorMag(tmpV2);
//        for (j=0; j<3; j++) RAY[i][j] = tmpV2[j];
//        RAY[i][3] = tmp1;
//    }
//    
//    vectorZero(FF);
//    for (i=0; i<p.NumFaces; i++)
//    {
//        for (j=0; j<3; j++)
//        {
//            vectorSet(tmpV1, RAY[p.Faces[i][j]][0], RAY[p.Faces[i][j]][1], RAY[p.Faces[i][j]][2]);
//            vectorScale(tmpV1, 1.0 / RAY[p.Faces[i][j]][3], tmpVS[j]);
//        }
//        vectorCross(tmpVS[1], tmpVS[2], tmpV2);
//        tmp1 = vectorDot(tmpV2, tmpVS[0]);
//        tmp2 = vectorDot(tmpVS[0], tmpVS[1]) + vectorDot(tmpVS[1], tmpVS[2]) + vectorDot(tmpVS[2], tmpVS[0]) + 1.0;
//        omegaf = 2.0 * atan2(tmp1, tmp2);
//        vectorSet(tmpV2, p.FaceNormVecs[i][0], p.FaceNormVecs[i][1], p.FaceNormVecs[i][2]);
//        tmp1 = vectorDot(tmpV1, tmpV2);
//        vectorScale(tmpV2, tmp1*omegaf, tmpV2);
//        vectorAdd(FF, tmpV2, FF);
//    }
//
//    vectorZero(EE);
//    for (i=0; i<p.NumEdges; i++)
//    {
//        tmp1 = RAY[p.Edges[i][0]][3] + RAY[p.Edges[i][1]][3];
//        le = log((tmp1 + p.EdgeLens[i]) / (tmp1 - p.EdgeLens[i]));
//        vectorSet(tmpV1, RAY[p.Edges[i][0]][0], RAY[p.Edges[i][0]][1], RAY[p.Edges[i][0]][2]);
//        for (j=0; j<2; j++)
//        {
//            vectorSet(tmpV2, p.EdgeNormVecs[i][3*j], p.EdgeNormVecs[i][3*j+1], p.EdgeNormVecs[i][3*j+2]);
//            tmp1 = vectorDot(tmpV2, tmpV1);
//            vectorSet(tmpV2,p.FaceNormVecs[p.Edges[i][j+2]][0], p.FaceNormVecs[p.Edges[i][j+2]][1], p.FaceNormVecs[p.Edges[i][j+2]][2]);
//            vectorScale(tmpV2, tmp1*le, tmpV2);
//            vectorAdd(EE, tmpV2, EE);
//        }
//    }
//
//    vectorSub(FF, EE, F);
//
//    for (i=0; i<p.NumVerts; i++) delete [] RAY[i];
//    delete [] RAY;
//}
//
//
//
//// r, F in BXbYbZb, no G and m_ejecta multiplied
//void GravBeta(CLUSTER &c, Vector r, Vector F)
//{
//    int i;
//    double tmpd;
//    Vector tmpV;
//    
//    vectorZero(F);
//    for (i=0; i<c.Number; i++)
//    {
//        vectorSet(tmpV, c.Positions[i][0], c.Positions[i][1], c.Positions[i][2]);
//        vectorSub(r, tmpV, tmpV);
//        //
//        tmpd = vectorMag(tmpV);
//        assert(tmpd>0.0);
//        tmpd = - c.Mass[i] / (tmpd*tmpd*tmpd);
//        //
//        vectorScale(tmpV, tmpd, tmpV);
//        vectorAdd(F, tmpV, F);
//    }
//}














