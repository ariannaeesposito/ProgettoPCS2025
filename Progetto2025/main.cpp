int main()
{
	#include <iostream>
	#include "Polygon.hpp"
	#include "Utils.hpp"
	#include "UCDUtilities.hpp"
	
	using namespace std;
	using namespace Eigen;
	using namespace PolygonalLibrary;
	
	int main(int dim ,char* argv[])
	{
		PolygonalMesh mesh;
		cout << input_solido_platonico(mesh, dim, argv) << endl;
	
		if(!input_solido_platonico(mesh, dim, argv));
		{
			cerr << "Errorore" << endl;
			return 1;
		}
	
		if(!ImportMesh(mesh))
		{
			cerr << "file not found" << endl;
			return 1;
		}
		else
		{
			cout << "File imported successfully" << endl;
		}
		
		Gedim::UCDUtilities utilities;
		utilities.ExportPoints("./Cell0Ds.inp",
								mesh.M0D);
		utilities.ExportSegments("./Cell1Ds.inp",
								  mesh.M0D,
								  mesh.M1D);
	
	
	return 0;
	
	}
	
	}
