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

        MatrixXd M0D; // 3 * ID punti --> coordinate 
        MatrixXi M1D; //spigoli poligono platonico
        MatrixXi M1D_triangolini; // matrice spigoli triangolini 2*num_spigolini , per ogni j-esimo spigolino ho i due punti sulle righe 
        MatrixXi M_pt_spigoli; //matrice che conterrà l'ID del j-esimo punto ( anche generato dalla mesh ) sul i-esimo spigolo.  matrice id_spigolo*num_punti(sullo spigolo) 
        MatrixXi M1D_spigoli_intermedi; //matrice di dimensione (numero_spigoli × d) (d = NUMERO DI SEGMENTI su ogni lato)
        MatrixXi M2D; // coordinate e spigoli triangolini

        vector<vector<unsigned int>> M2D_vertici; // essendo un vector di vector oer accedere ai singoli elt devo utilizzare [][] e non ( , )
        vector<vector<unsigned int>> M2D_spigoli; 
} ;

}