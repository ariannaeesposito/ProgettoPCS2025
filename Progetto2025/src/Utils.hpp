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
 
}
