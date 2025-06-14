#include <iostream>
#include "Polygon.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"
	
using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;
	
	Gedim::UCDUtilities utilities;

	int main(int argc ,char* argv[])
	{
		PolygonalMesh mesh, platonico, Dmesh;
		unsigned int classe;
		string poligono_platonico;
		
		if(input_solido_platonico(platonico, mesh, argc, argv, poligono_platonico, classe))
			cout << "Input imported successfully" << endl;
		else{
			cout << "Error in input" << endl;
			return 1;
		}
		if(ImportMesh(platonico, poligono_platonico))
			cout << "File imported successfully" << endl;
		
		if(Inizializzazione_vertici(mesh, classe, platonico))
			cout << "Vertices initialized successfully" << endl;
	
		if (classe == 1){
			if (Triangolazione_1_classe(mesh, platonico))
				cout << "Triangulation of class 1 completed successfully" << endl;
		}
		else{
			if (Triangolazione_2_classe(mesh, platonico))
				cout << "Triangulation of class 2 completed successfully" << endl;			
		}

	// 	//stampa su paraview Solido triangolato

		utilities.ExportPoints("./Triangolazione0Ds.inp",
								mesh.M0D);
		utilities.ExportSegments("./Triangolazione1Ds.inp",
								  mesh.M0D,
								  mesh.M1D);
						
		//stampa su paraview Solido proiettato sulla sfera

		Proiezione_sfera(mesh);

		utilities.ExportPoints("./ProiezioneSfera0Ds.inp",
								mesh.M0D);
		utilities.ExportSegments("./ProiezioneSfera1Ds.inp",
			  					mesh.M0D,
			  					mesh.M1D);

	if ( argc == 7 ){
		CamminoMinimo(mesh);
		}
		//stampa su paraview duale

		Duale(mesh, Dmesh);

		utilities.ExportPoints("./Duale0Ds.inp",
								Dmesh.M0D);
		utilities.ExportSegments("./Duale1Ds.inp",
			  					Dmesh.M0D,
			  					Dmesh.M1D);

		//stampa_solido(mesh, poligono_platonico);
		stampa_solido(platonico, poligono_platonico);	
		//stampa_solido(Dmesh, poligono_platonico);


		if ( argc == 7 ){
			CamminoMinimo(Dmesh);
		}
	return 0;
	}