#pragma once
#include <gtest/gtest.h>
#include "Utils.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
using namespace std;
using namespace PolygonalLibrary;

template <typename T>
void stampa_vector(const std::vector<T>& v) {
    std::cout << "[ ";
    for (const auto& elem : v) {
        std::cout << elem << " ";
    }
    std::cout << "]" << std::endl;
}

PolygonalMesh crea_mesh(char valenza, char d, unsigned int classe, string duale){
	PolygonalMesh mesh, platonico, Dmesh;
	string poligono_platonico;
	
	
	char b[2], c[2];
    if (classe == 1){
        b[0] = d; b[1] = '\0';
        c[0] = '0'; c[1] = '\0';
    }
    else{
        b[0] = d; b[1] = '\0';
        c[0] = d; c[1] = '\0';
    }
    char val[2];
    val[0] = valenza; val[1] = '\0';

    char v1[2] = {'0', '\0'};
    char v2[2] = {'3', '\0'};
	char v3[2] = {'6', '\0'};
    char v4[2] = {'7', '\0'};

    char* vec[7] = {v1, v2, val, b, c,v3,v4};
		if(!input_solido_platonico(platonico, mesh, 7, vec, poligono_platonico, classe))
		{
			cerr << "Error Input non imported successfully" << endl;
		}
		else
		{
			cout << "Input imported successfully" << endl;
		}

		if(ImportMesh(platonico, poligono_platonico))

		if(Inizializzazione_vertici(mesh, classe, platonico))

		if (classe == 1){
			if (Triangolazione_1_classe(mesh, platonico))
				cout << "Triangulation of class 1 completed successfully" << endl;
		}
		else{
			if (Triangolazione_2_classe(mesh, platonico))
				cout << "Triangulation of class 2 completed successfully" << endl;			
		}
	if (duale == "duale"){
		Duale(mesh, Dmesh);
		return Dmesh;
	}
	else{
	return mesh;
	}
}

unsigned int controlla_ordinamento(char valenza, unsigned int classe, string duale){
	char d = '6';

	PolygonalMesh mesh = crea_mesh(valenza,d, classe, duale);
	unsigned int counter = 0;
	unsigned int a;
	unsigned int b;
	unsigned int k;
	cout << "mesh pt: "<< mesh.V << "sp "<< mesh.E << "f "<< mesh.F;
	for (unsigned int i=0; i<mesh.M2D_vertici.size(); i++){
		//cout <<endl << "faccia: " << i;
		//stampa_vector(mesh.M2D_vertici[i]);
		for (unsigned int j=0; j<mesh.M2D_vertici[i].size()-1; j++){
			a = mesh.M1D(0,mesh.M2D_spigoli[i][j]);
			b = mesh.M1D(1,mesh.M2D_spigoli[i][j]);
			//cout << "spigolo: "<< j << " da " << a << " a " << b << endl;
			if (mesh.M2D_vertici[i][j] == a && mesh.M2D_vertici[i][j+1]==b){
				//cout << "yay j" << endl;
			}
			else if (mesh.M2D_vertici[i][j] == b && mesh.M2D_vertici[i][j+1]==a){
				//cout << "yay j" << endl;
			}
			else{
				//cout << "non ordinato " << i << " "<< j <<endl;
				//cout << "non ordinato " << i << " " << j << endl;
				counter++;
			}
		}
		k = mesh.M2D_vertici[i].size()-1;
		a = mesh.M1D(0,mesh.M2D_spigoli[i][k]);
		b = mesh.M1D(1,mesh.M2D_spigoli[i][k]);
		//cout << "spigolo: "<< k << " da " << a << " a " << b << endl;
		if (mesh.M2D_vertici[i][k] == a && mesh.M2D_vertici[i][0]==b){
				//cout << "yay k" << endl;
			}
			else if (mesh.M2D_vertici[i][k] == b && mesh.M2D_vertici[i][0]==a){
				//cout << "yay k" << endl;
			}
			else{
				//cout << "non ordinato " << i << " " << k << endl;
				counter++;
			}
	}
	
	return counter;
}

unsigned int controlla_conteggio(char valenza, unsigned int classe, string duale){
	char d = '3';

	PolygonalMesh mesh = crea_mesh(valenza, d, classe, duale);
	unsigned int counter =0;
	unsigned int errori = 0;
	for (unsigned int i=0; i<mesh.M1D.cols(); i++){
		for (unsigned int k=0; k<mesh.M2D_vertici.size();k++){
			for (unsigned int j =0;j<mesh.M2D_vertici[k].size(); j++){
			
		if (mesh.M2D_spigoli[k][j]==i){
			counter++;
		}
		}
		}
		if (counter == 2){
		}
		else{
			errori++;
		}
		counter = 0;
	}
	return errori;
}


TEST(PolygonalTest, TestOrdinamentoM2D){
	char v1='3';
	char v2='4';
	char v3='5';
	unsigned int counter31n = controlla_ordinamento(v1,1,"");
	unsigned int counter32n = controlla_ordinamento(v1,2,"");
	unsigned int counter31d = controlla_ordinamento(v1,1,"duale");
	unsigned int counter32d = controlla_ordinamento(v1,2,"duale");
	unsigned int counter41n = controlla_ordinamento(v2,1,"");
	unsigned int counter42n = controlla_ordinamento(v2,2,"");
	unsigned int counter41d = controlla_ordinamento(v2,1,"duale");
	unsigned int counter42d = controlla_ordinamento(v2,2,"duale");
	unsigned int counter51n = controlla_ordinamento(v3,1,"");
	unsigned int counter52n = controlla_ordinamento(v3,2,"");
	unsigned int counter51d = controlla_ordinamento(v3,1,"duale");
	unsigned int counter52d = controlla_ordinamento(v3,2,"duale");
	
	ASSERT_EQ(counter31n,0);
	ASSERT_EQ(counter32n,0);
	cout << "123" <<endl;
	ASSERT_EQ(counter31d,0);
	ASSERT_EQ(counter32d,0);
	ASSERT_EQ(counter41n,0);
	ASSERT_EQ(counter42n,0);
	cout << "456" <<endl;
	ASSERT_EQ(counter41d,0);
	ASSERT_EQ(counter42d,0);
	ASSERT_EQ(counter51n,0);
	ASSERT_EQ(counter52n,0);
	cout << "789" <<endl;
	ASSERT_EQ(counter51d,0);
	ASSERT_EQ(counter52d,0);
}
TEST(PolygonalTest, TestConteggioSpigoli){
	char v1='3';
	char v2='4';
	char v3='5';
	unsigned int counter31n = controlla_conteggio(v1,1,"");
	unsigned int counter32n = controlla_conteggio(v1,2,"");
	unsigned int counter31d = controlla_conteggio(v1,1,"duale");
	unsigned int counter32d = controlla_conteggio(v1,2,"duale");
	unsigned int counter41n = controlla_conteggio(v2,1,"");
	unsigned int counter42n = controlla_conteggio(v2,2,"");
	unsigned int counter41d = controlla_conteggio(v2,1,"duale");
	unsigned int counter42d = controlla_conteggio(v2,2,"duale");
	unsigned int counter51n = controlla_conteggio(v3,1,"");
	unsigned int counter52n = controlla_conteggio(v3,2,"");
	unsigned int counter51d = controlla_conteggio(v3,1,"duale");
	unsigned int counter52d = controlla_conteggio(v3,2,"duale");
	
	ASSERT_EQ(counter31n,0);
	ASSERT_EQ(counter32n,0);
	ASSERT_EQ(counter31d,0);
	ASSERT_EQ(counter32d,0);
	ASSERT_EQ(counter41n,0);
	ASSERT_EQ(counter42n,0);
	ASSERT_EQ(counter41d,0);
	ASSERT_EQ(counter42d,0);
	ASSERT_EQ(counter51n,0);
	ASSERT_EQ(counter52n,0);
	ASSERT_EQ(counter51d,0);
	ASSERT_EQ(counter52d,0);
}

TEST(PolygonalTest, TestTriangolazione1){
vector<vector<unsigned int>> m1vn = {
	{0,4,5},{5,4,7},{4,1,7},{5,7,2},{0,5,6},
	{6,5,8},{5,2,8},{6,8,3},{0,6,4},{4,6,9},
	{6,3,9},{4,9,1},{1,7,9},{9,7,8},{7,2,8},{9,8,3}
};

vector<vector<unsigned int >> m1sn = {
	{0,12,2},{12,14,13},{1,6,14},{13,7,3},{2,15,4},
	{15,17,16},{3,8,17},{16,9,5},{4,18,0},{18,20,19},
	{5,10,20},{19,11,1},{6,21,11},{21,23,22},{7,8,23},{22,9,10}
};

	char valenza = '3';
	char d = '2';

	PolygonalMesh mesh = crea_mesh(valenza, d, 1, "");
	
	for (unsigned int i=0; i<mesh.M2D_vertici.size();i++){
		for (unsigned int j=0; j<mesh.M2D_vertici[i].size();j++){
			ASSERT_EQ(mesh.M2D_vertici[i][j], m1vn[i][j]);
		}
	}
	for (unsigned int i=0 ; i<mesh.M2D_spigoli.size();i++){
		for (unsigned int j=0; j<mesh.M2D_spigoli[i].size();j++){
			ASSERT_EQ(mesh.M2D_spigoli[i][j],m1sn[i][j]);
		}
	}
}

TEST(PolygonalTest, TestTriangolazione2){
	vector<vector<unsigned int >> m2vn = {
	{0,4,10},{4,1,10},{0,10,5},{5,10,2},{10,2,7},{10,7,1},{0,5,11},{5,2,11},
	{0,11,6},{6,11,3},{11,3,8},{11,8,2},{0,6,12},{6,3,12},{0,12,4},{4,12,1},
	{12,1,9},{12,9,3},{1,7,13},{7,2,13},{1,13,9},{9,13,3},{13,3,8},{13,8,2}
	};

	vector<vector<unsigned int >> m2sn = {
	{0,15,14},{1,16,15},{14,12,2},{12,13,3},{13,7,17},{17,6,16},{2,21,20},{3,22,21},
	{20,18,4},{18,19,5},{19,9,23},{23,8,22},{4,27,26},{5,28,27},{26,24,0},{24,25,1},
	{25,11,29},{29,10,28},{6,33,32},{7,34,33},{32,30,11},{30,31,10},{31,9,35},{35,8,34}
	};

	char valenza = '3';
	char d = '1';

	PolygonalMesh mesh = crea_mesh(valenza, d, 2, "");
	
	for (unsigned int i=0; i<mesh.M2D_vertici.size();i++){
		for (unsigned int j=0; j<mesh.M2D_vertici[i].size();j++){
			ASSERT_EQ(mesh.M2D_vertici[i][j],m2vn[i][j]);
		}
	}
	for (unsigned int i=0 ; i<mesh.M2D_spigoli.size();i++){
		for (unsigned int j=0; j<mesh.M2D_spigoli[i].size();j++){
			ASSERT_EQ(mesh.M2D_spigoli[i][j],m2sn[i][j]);
		}
	}
}

TEST(PolygonalTest, TestTriangolazioneDuale){
	vector<vector<unsigned int>> m1vd = {
	{0,2,1},
	{0,2,3},
	{0,1,3},
	{1,2,3}
	};

	vector<vector<unsigned int>> m1sd = {
	{0,2,1},
	{0,5,3},
	{1,4,3},
	{2,5,4}
	};
	char valenza = '3';
	char d = '1';

	PolygonalMesh mesh = crea_mesh(valenza, d, 1, "duale");
	for (unsigned int i=0; i<mesh.M2D_vertici.size();i++){
		for (unsigned int j=0; j<mesh.M2D_vertici[i].size();j++){
			ASSERT_EQ(mesh.M2D_vertici[i][j],m1vd[i][j]);
		}
	}
	for (unsigned int i=0; i<mesh.M2D_spigoli.size();i++){
		for (unsigned int j=0; j<mesh.M2D_spigoli[i].size();j++){
			ASSERT_EQ(mesh.M2D_spigoli[i][j],m1sd[i][j]);
		}
	}
}


TEST(PolygonalTest, TestCamminoMinimo){
	char valenza = '3';
	char d = '2';

	cout << "prova"<< endl;
	vector<unsigned int> percorso1={7,4,6};
	vector<unsigned int> percorso2={7,5,6};
	
	PolygonalMesh mesh = crea_mesh(valenza, d, 1, "");
	CamminoMinimo(mesh);
	EXPECT_TRUE(mesh.percorso == percorso1 || mesh.percorso == percorso2);
}