#include "Utils.hpp"
#include "UCDUtilities.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;
using namespace Eigen;

namespace PolygonalLibrary

{
bool input_solido_platonico(PolygonalMesh& mesh, int dim ,char* argv[]){
cout << "ciao"<< endl;
cout << "dim = " << dim << endl;
if (dim !=5)
    {
        cout << "Error: the number of vertices is not correct" << endl;
        return false;
    }
    int p, q, b, c;  
    istringstream convert0(argv[1]);
    convert0 >> p;
    istringstream convert1(argv[2]);
    convert1 >> q;
    istringstream convert2(argv[3]);
    convert2 >> b;
    istringstream convert3(argv[4]);
    convert3 >> c;
    cout << "p = " << p << endl;
    cout << "q = " << q << endl;
    cout << "b = " << b << endl;
    cout << "c = " << c << endl;
    int flag = 0;
if (p==3){
    switch (q)
    {
    case 3:
        mesh.nomefile0 = "../PolygonalMesh/Cell0Ds_tetraedro.csv";
        mesh.nomefile1 = "../PolygonalMesh/Cell1Ds_tetraedro.csv"; 
        mesh.nomefile2 = "../PolygonalMesh/Cell2Ds_tetraedro.csv";
        break;
    case 4:
        mesh.nomefile0 = "../PolygonalMesh/Cell0Ds_ottaedro.csv";
        mesh.nomefile1 = "../PolygonalMesh/Cell1Ds_ottaedro.csv";
        mesh.nomefile2 = "../PolygonalMesh/Cell2Ds_ottaedro.csv";
        break;
    case 5:
        mesh.nomefile0 = "../PolygonalMesh/Cell0Ds_icosaedro.csv";
        mesh.nomefile1 = "../PolygonalMesh/Cell1Ds_icosaedro.csv";
        mesh.nomefile2 = "../PolygonalMesh/Cell2Ds_icosaedro.csv";
        break;
    default:
        cout << "errore f" << endl;
        return false;
        break;
    }
    flag = 1;
}
cout << "flag = " << flag << endl;
if (q==3 && p!=3){
switch (p)
{
case 4:
    mesh.nomefile0 = "../PolygonalMesh/Cell0Ds_tetraedro.csv";
    mesh.nomefile1 = "../PolygonalMesh/Cell1Ds_tetraedro.csv";
    mesh.nomefile2 = "../PolygonalMesh/Cell2Ds_tetraedro.csv";
    break;
case 5:
    mesh.nomefile0 = "../PolygonalMesh/Cell0Ds_icosaedro.csv";
    mesh.nomefile1 = "../PolygonalMesh/Cell1Ds_icosaedro.csv";
    mesh.nomefile2 = "../PolygonalMesh/Cell2Ds_icosaedro.csv";
    break;
default:
    cout << "errore q" << endl;
    return false;
}
flag = 1;
}
if (flag == 0){
    cerr << p << "," << q << "non è un solido platonico"<< endl;
    cout << "errore flag pq" << endl;
    return false;
}
int classe;
flag = 0;
if (b==0 && c >= 1){
    mesh.classe = 1;
	mesh.d=c;
    flag = 1;
}
if (c==0 && b >= 1){
	mesh.classe=1;
	mesh.d=b;
	flag=1;
}

if ((b==c && b >= 1)){
    mesh.classe = 2;
    flag = 1;
}
if (flag == 0){
    cerr << b << "," << c << " non sono triangolazioni possibili" << endl;
    cout << "errore flag bc" << endl;
    return false;
}
cout << "arrivati qui" << endl;

return true;
}

bool ImportMesh(PolygonalMesh& mesh)
    {
    if(!ImportCell0Ds(mesh))
        return false;
    if(!ImportCell1Ds(mesh))
        return false;   
    if(!ImportCell2Ds(mesh))
        return false;

    return true;     
    }

bool ImportCell0Ds(PolygonalMesh& mesh)
{
    ifstream file_Cell0Ds(mesh.nomefile0);
    if (!file_Cell0Ds.is_open())
    {
        cerr << "Error opening file "<< mesh.nomefile0 << endl;
        return false;
    }
    string riga;
    list<string> lista_dim;
        while (getline(file_Cell0Ds, riga))
    {
        lista_dim.push_back(riga);
    }
    file_Cell0Ds.close();
    lista_dim.pop_front();//tolgo l'intestazione
    mesh.Dim0D = lista_dim.size();

    if (mesh.Dim0D == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }

    mesh.M0D = Eigen::MatrixXd::Zero(3, mesh.Dim0D); 
// voglio stampare la lista

    for (string& li : lista_dim) // itero una lista
    {
        replace(li.begin(), li.end(), ';', ' ');
        istringstream convertitore(li);
        unsigned int id,marker;
        convertitore >> id >> marker >> mesh.M0D(0,id) >> mesh.M0D(1,id) >> mesh.M0D(2,id)   ;
    if (marker != 0)
    {
        auto it = mesh.marker0D.find(marker);
        if (it != mesh.marker0D.end())//verrifico se il marker è già presente // se metto == campio if con else
              (*it).second.push_back(id);//memoriazza gli id nei marker gia presenti
        else
            mesh.marker0D.insert({marker, {id}});//memoriazza gli id nei marker
        }

    }
return true;
}

bool ImportCell1Ds(PolygonalMesh& mesh)

    {
        {
            ifstream file_Cell1Ds(mesh.nomefile1);
            if (!file_Cell1Ds.is_open())
            {
                cerr << "Error opening file "<< mesh.nomefile1 << endl;
                return false;
            }
            string riga;
            list<string> lista_dim;
            while (getline(file_Cell1Ds, riga))
            {
                lista_dim.push_back(riga);
            }

            file_Cell1Ds.close();
            lista_dim.pop_front(); //tolgo l'intestazione

            mesh.Dim1D = lista_dim.size();
            if (mesh.Dim1D == 0)
            {
                cerr << "There is no cell 1D" << endl;
                return false;
            }

            mesh.M1D = MatrixXi::Zero(2, mesh.Dim1D); ; 
            for (string& li : lista_dim) // itero una lista
            {
              replace(li.begin(), li.end(), ';', ' ');
              istringstream convertitore(li);
              unsigned int id,marker;
              convertitore >> id >> marker >> mesh.M1D(0,id)>> mesh.M1D(1,id); 
              if (marker != 0)
              {
                auto it = mesh.marker1D.find(marker);
                if (it != mesh.marker1D.end())// 
                    (*it).second.push_back(id);//memoriazza gli id nei marker
                else
                    mesh.marker1D.insert({marker, {id}});//memoriazza g
              }
            } 
        }
        return true;
    }   

bool ImportCell2Ds(PolygonalMesh& mesh)
    {
        ifstream file_Cell2Ds(mesh.nomefile2);
        if (!file_Cell2Ds.is_open())
        {
            cerr << "Error opening file "<< mesh.nomefile2 << endl;
            return false;
        }
        
        string riga;
        list<string> lista_dim;
        while (getline(file_Cell2Ds, riga))
        {
            lista_dim.push_back(riga);
        }
        file_Cell2Ds.close();
        lista_dim.pop_front();//tolgo l'intestazione
        mesh.Dim2D = lista_dim.size();
        if (mesh.Dim2D == 0)
        {
            cerr << "There is no cell 2D" << endl;
            return false;
        }

        for (string& li : lista_dim) // itero una lista
        
        {
          replace(li.begin(), li.end(), ';', ' ');
          istringstream convertitore(li);
          unsigned int id,marker,n_vertici_spigoli;
          convertitore >> id >> marker >> n_vertici_spigoli;
          if (marker != 0)
              {
                auto it = mesh.marker2D.find(marker);
                if (it != mesh.marker2D.end())// 
                    (*it).second.push_back(id);//memoriazza gli id nei marker
                else
                    mesh.marker2D.insert({marker, {id}});//memoriazza g
              }
// lavoro su matrice vertici
          mesh.M2D_vertici.reserve(mesh.Dim2D);
          vector<int> linea1;
          linea1.resize(n_vertici_spigoli);// uso resize per allocare la memoria
          for (unsigned int j = 0; j < n_vertici_spigoli ; j++)
          {
            convertitore >> linea1[j];
          }
          mesh.M2D_vertici.push_back(linea1);
          
        
          // scarto numero spigoli         
          convertitore >> n_vertici_spigoli;
          // lavoro su matrice spigoli
          mesh.M2D_spigoli.reserve(mesh.Dim2D);    
          vector<int> linea2;
          linea2.reserve(n_vertici_spigoli);
          for (unsigned int j = 0; j < n_vertici_spigoli ; j++)
          {
            convertitore >> linea2[j];
          }

          mesh.M2D_spigoli.push_back(linea2);
          
        }
        
    return true;
    }

void triangolazione(PolygonalMesh& mesh){
	if (mesh.classe ==1){
		triangolazione1(mesh);
	}
	else{
		triangolazione2(mesh);
	}
}

void triangolazione1(PolygonalMesh& mesh){
	
}

void triangolazione2(PolygonalMesh& mesh){
	
}

unsigned int accedimappa(TriangularMesh& mesh, PolygonalMesh& Pmesh, Vector3d<unsigned int> coord){
	double tol = pow(10,-15); \\definiamo tolleranza
;	for (const auto& coppia : mesh.coordinate_punti){ \\iteriamo sulla mappa
		diff = coppia.first-coord; \\calcoliamo la differnaza tra i due set di coordinate
		if (diff.norm() < tol){ \\se sono vicine allora le consideriamo le stesse
			return coppia.second \\restituiamo l'id associato a quelle coordinate
		}
	}
	\\ se non c'è nella mappa lo creiamo nuovo
	unsigned int dim = Pmesh.Dim0D; \\dim corrente
	mesh.coordinate_punti[coord] = dim+1; \\cell0ds già piena stiamo, stiamo aggiungendo alla mappa una chiave e l'id relativo
	Pmesh.dim0D += 1; \\aumentiamo il conteggio dei punti già creati
	Pmesh.Cell0DsCoordinates(0, dim+1)=coord(0); \\mettiamo il nuovo punto dentro la polygonal mesh
	Pmesh.Cell0DsCoordinates(1, dim+1)=coord(1);
	Pmesh.Cell0DsCoordinates(2, dim+1)=coord(2);
}

}
