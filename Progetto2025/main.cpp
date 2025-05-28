#include <iostream>
	#include "Polygon.hpp"
	#include "Utils.hpp"
	#include "UCDUtilities.hpp"
	
	using namespace std;
	using namespace Eigen;
	using namespace PolygonalLibrary;
	
	int main(int argc ,char* argv[])
	{
		PolygonalMesh mesh;
		
		if(!input_solido_platonico(mesh, argc, argv))
		{
			cerr << "Error Input non imported successfully" << endl;
			return 1;
		}
		else
		{
			cout << "Input imported successfully" << endl;
		}

		if(!ImportMesh(mesh))
		{
			cerr << "Error File not found" << endl;
			return 1;
		}
		else
		{
			cout << "File imported successfully" << endl;
		}

		Inizializzazione_vertici(mesh);
		Inizializzazione_punti_interni(mesh);
		Proiezione_sfera(mesh);
		stampa_geodetico(mesh);
		
		Gedim::UCDUtilities utilities;
		utilities.ExportPoints("./Cell0Ds.inp",
								mesh.M0D);
		utilities.ExportSegments("./Cell1Ds.inp",
								  mesh.M0D,
								  mesh.M1D_triangolini);
	
		PolygonalMesh Dmesh;
		Duale(mesh, Dmesh);

		//Proiezione_sfera(mesh);

		utilities.ExportPoints("./Cell0DsD.inp",
								Dmesh.M0D);
		utilities.ExportSegments("./Cell1DsD.inp",
								  mesh.M0D,
								  Dmesh.M1D_triangolini);
	return 0;

	
	
	}