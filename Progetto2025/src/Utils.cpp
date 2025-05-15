#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "UCDUtilities.hpp"


using namespace Eigen;

namespace PolygonalLibrary 
{

bool ImportMesh(PolygonalMesh& mesh,const string& nomefile0,const string& nomefile1,const string& nomefile2)
{

    if(!ImportCell0Ds(mesh,nomefile0))
        return false;

    if(!ImportCell1Ds(mesh,nomefile1))
        return false;

    if(!ImportCell2Ds(mesh,nomefile2))
        return false;

    return true;

}
