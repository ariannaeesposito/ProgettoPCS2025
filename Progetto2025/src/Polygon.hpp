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

        unsigned int classe;
        unsigned int d , q , p;
        double lunghezza_lato_triangolino; // numero di spigoli
        // V, E, F :  n di vertici totali (con mesh), n spigoli totali (con mesh), n facce di tutti i triangoli (con mesh)
        unsigned int V, E, F, T;   

        string nomefile0;
        string nomefile1;
        string nomefile2;

        map <unsigned int, list<unsigned int>> marker0D;
        map <unsigned int, list<unsigned int>> marker1D;
        map <unsigned int, list<unsigned int>> marker2D;

        MatrixXd M0D; //coordinate
        MatrixXi M1D; //estremi
        MatrixXi M1D_triangolini; // triangolini

        MatrixXd M2D; // triangolini 

        vector<vector<unsigned int>> M2D_vertici; // essendo un vector di vector oer accedere ai singoli elt devo utilizzare [][] e non ( , )
        vector<vector<unsigned int>> M2D_spigoli; 
       

} ;

}