#pragma once

#include <iostream>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolygonalLibrary{

    struct PolygonalMesh{
     
        unsigned int p, q, b, c, d;
        // V, E, F :  n di vertici totali (con mesh), n spigoli totali (con mesh), n facce di tutti i triangoli (con mesh)
        unsigned int V, E, F;   
        
        unsigned int nodo_i, nodo_f; // vertici iniziale e finale per il percorso geodetico
		vector<unsigned int> percorso; // vector con gli id del percorso

        MatrixXd M0D; // 3 * ID punti --> coordinate 
        MatrixXi M1D; //spigoli poligono platonico
        vector<vector<unsigned int>> M2D_vertici; // essendo un vector di vector oer accedere ai singoli elt devo utilizzare [][] e non ( , )
        vector<vector<unsigned int>> M2D_spigoli; 
        vector<unsigned int> M3D ; 

        MatrixXi M1D_spigoli_intermedi; //matrice di dimensione (numero_spigoli × d) (d = NUMERO DI SEGMENTI su ogni lato)
        MatrixXi M_pt_spigoli; //matrice che conterrà l'ID del j-esimo punto ( anche generato dalla mesh ) sul i-esimo spigolo.  matrice id_spigolo*num_punti(sullo spigolo) 


} ;

}
