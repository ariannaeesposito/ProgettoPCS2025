
#pragma once

#include <iostream>
#include <map>
#include <vector>
#include "Eigen/Eigen"


using namespace std;
using namespace Eigen;

namespace PolygonalLibrary{

    struct PolygonalMesh{
        unsigned int Dim0D;
        unsigned int Dim1D;
        unsigned int Dim2D;// lunghezza della matrice 2D (numero di righe, dunque numero di poligoni)
        
        vector<unsigned int> Cell0DsId;
        vector<unsigned int> Cell1DsId;
        vector<unsigned int> Cell2DsId;

        MatrixXd Cell0DsCoordinates;

        unsigned int classe;
        unsigned int d , q , p ;
        
        string nomefile0;
        string nomefile1;
        string nomefile2;

        map <unsigned int, list<unsigned int>> marker0D;
        map <unsigned int, list<unsigned int>> marker1D;
        map <unsigned int, list<unsigned int>> marker2D;

        Eigen::MatrixXd M0D;
        Eigen::MatrixXi M1D; 
        vector<vector< int>> M2D_vertici;

        vector<vector< int>> M2D_spigoli; //facce
    };

};
