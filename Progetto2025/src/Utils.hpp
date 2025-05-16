#pragma once

#include<vector>
#include <iostream>
#include "PolygonalMesh.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Gedim;

namespace PolygonalLibrary
{

bool ImportMesh(PolygonalMesh& mesh,const string& nomefile0,const string& nomefile1,const string& nomefile2);
bool input_solido_platonico(PolygonalMesh& mesh, int dim ,char* argv[]);
 
}
