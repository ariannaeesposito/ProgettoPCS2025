#pragma once

#include<vector>
#include <iostream>
#include "Polygon.hpp"
#include "Triangle.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Gedim;

namespace PolygonalLibrary
{
    bool ImportMesh(PolygonalMesh& mesh);
    bool input_solido_platonico(PolygonalMesh& mesh, int dim ,char* argv[]);

    bool ImportCell0Ds(PolygonalMesh& mesh);
    bool ImportCell1Ds(PolygonalMesh& mesh);
    bool ImportCell2Ds(PolygonalMesh& mesh);

 
    bool Inizializzazione_vertici( PolygonalMesh& Pmesh);
    bool Inizializzazione_punti_interni(PolygonalMesh& Pmesh);
    bool Proiezione_sfera(PolygonalMesh& Pmesh);
	bool Duale(PolygonalMesh& Pmesh, PolygonalMesh& Dmesh);
    bool stampa_geodetico(PolygonalMesh& Pmesh);
    bool Inizializzazione_punti_interni_classe2(PolygonalMesh& Pmesh);
    bool crea_triangolo(PolygonalMesh& Pmesh, const unsigned int& id_triangolo ,const unsigned int& id_pt_1, const unsigned int& id_pt_2, const unsigned int& id_pt_3, const unsigned int& id_sp_1, const unsigned int& id_sp_2, const unsigned int& id_sp_3);



    // unsigned int accedimappa(TriangularMesh& mesh, PolygonalMesh& Pmesh, Vector3d<unsigned int> coord);
    //void triangolazione(PolygonalMesh& mesh);
    // void triangolazione1(PolygonalMesh& mesh);
    // void triangolazione2(PolygonalMesh& mesh);

}