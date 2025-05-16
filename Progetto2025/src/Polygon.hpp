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
        
        unsigned int classe;
        
        string nomefile0;
        string nomefile1;
        string nomefile2;

        map <unsigned int, list<unsigned int>> marker0D;
        map <unsigned int, list<unsigned int>> marker1D;
        map <unsigned int, list<unsigned int>> marker2D;

        Eigen::MatrixXd M0D;
        Eigen::MatrixXi M1D; 
        vector<vector< int>> M2D_vertici;
        vector<vector< int>> M2D_spigoli;
    };

};
