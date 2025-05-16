#pragma once

#include<vector>
#include <iostream>
#include "PolygonalMesh.hpp"
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

    void triangolazione(PolygonalMesh& mesh);
    void triangolazione1(PolygonalMesh& mesh);
    void triangolazione2(PolygonalMesh& mesh);

    unsigned int accedimappa(TriangularMesh& mesh, PolygonalMesh& Pmesh, Vector3d<unsigned int> coord);
 
 
}
