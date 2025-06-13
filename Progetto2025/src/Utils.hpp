#pragma once
#include "Polygon.hpp"
#include "Triangle.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Gedim;

namespace PolygonalLibrary

{
    bool ImportMesh(PolygonalMesh& poligono, const string& poligono_platonico);
    bool ImportCell0Ds(PolygonalMesh& poligono, const string& poligono_platonico);
    bool ImportCell1Ds(PolygonalMesh& poligono, const string& poligono_platonico);
    bool ImportCell2Ds(PolygonalMesh& poligono, const string& poligono_platonico);

    bool input_solido_platonico(PolygonalMesh& poligono, PolygonalMesh& mesh, int argc ,char* argv[], string& poligono_platonico, unsigned int& classe);
    bool Inizializzazione_vertici( PolygonalMesh& Pmesh, const unsigned int& classe, PolygonalMesh& poligono);

    bool Triangolazione_1_classe(PolygonalMesh& Pmesh, PolygonalMesh& poligono);
    bool Triangolazione_2_classe(PolygonalMesh& Pmesh, PolygonalMesh& poligono);

    bool Proiezione_sfera(PolygonalMesh& Pmesh);

	bool Duale(PolygonalMesh& Pmesh, PolygonalMesh& Dmesh);

    bool stampa_solido(PolygonalMesh& Pmesh,const string& poligono_platonico);
   
    bool CamminoMinimo(PolygonalMesh& Pmesh);
   //bool CamminoMinimo2(PolygonalMesh& Pmesh);

    bool crea_triangolo(PolygonalMesh& Pmesh ,const unsigned int& id_pt_1, const unsigned int& id_pt_2, const unsigned int& id_pt_3, const unsigned int& id_sp_1, const unsigned int& id_sp_2, const unsigned int& id_sp_3);

}