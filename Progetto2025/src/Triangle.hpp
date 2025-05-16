#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PoligonalLibrary {

DA CAMBIARE TUTTOOOOO
CI SONO SOLO POCHE COSE FATTE BENE

struct TriangularMesh
{
   
    
    vector<unsigned int> Cell0DsId;
    vector<unsigned int> Cell1DsId;
    vector<unsigned int> Cell2DsId;
    vector< unsigned int > Cell0DsMarker;
    vector< unsigned int > Cell1DsMarker;

    MatrixXd Cell0DsCoordinates;
    MatrixXi Cell1DsExtremes;

    vector<vector<unsigned int>> Cell2DsVertices;
    vector<vector<unsigned int>> Cell2DsEdges;

    map<Vector3d<unsigned int> , unsigned int> coordinate_punti;  

    map<unsigned int,list <unsigned int>> MarkerCell1Ds;  //mappa marker → lista di segmenti (Cell1D)
     
};
}
