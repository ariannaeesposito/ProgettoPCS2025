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
#include <queue>

using namespace std;
using namespace Eigen;

namespace PolygonalLibrary

{
    
bool input_solido_platonico(PolygonalMesh& poligono, PolygonalMesh& mesh, int argc ,char* argv[],  string& poligono_platonico, unsigned int& classe){

    if ((argc != 7) && (argc != 5)){
        cout << "Error: the number of input is not correct" << endl  ;
        return false;
    }

    unsigned int b, c, T ,flag = 0; 

    stringstream convert;
    for ( int a = 1 ; a < argc ; a++ )
        convert << argv[a] << ' ';

    convert >> mesh.p >> mesh.q >> b >> c; //numero di lati per faccia

    if ((argc == 7)){
        convert >> mesh.nodo_i >> mesh.nodo_f; // vertice iniziale per il cammino minimo e // vertice iniziale per il cammino minimo
    }

    mesh.d=max(b,c);

if (mesh.p==3){ //facce triangolari
    switch (mesh.q)
    {
    case 3: //tetraedro
        poligono_platonico = "tetraedro";
        poligono.V = 4;
        poligono.E = 6;
        poligono.F = 4;
        break;
    case 4: //ottaedro 
        poligono_platonico = "ottaedro";
        poligono.V = 6;
        poligono.E = 12;
        poligono.F = 8;
        break;
    case 5: //icosaedro
        poligono_platonico = "icosaedro";
        poligono.V = 12;
        poligono.E = 30;
        poligono.F = 20;
        break;
    default:
        cerr << " p ammissibile ma q no" << endl;
        return false;
        break;
    }

    flag = 1;
} 

if (mesh.q==3 && mesh.p!=3){ //duale di un solido a facce triangolari

    switch (mesh.p)
    {
    case 4: //cubo -> tetraedro duale
        poligono_platonico = "tetraedro";
        poligono.V = 4;
        poligono.E = 6;
        poligono.F = 4;
        break;
    case 5: //dodecaedro -> icosaedro duale
        poligono_platonico = "icosaedro";
        poligono.V = 12;
        poligono.E = 30;
        poligono.F = 20;
        break;
    default:
        cerr << "q ammissibile ma p non lo è" << endl;
        return false;
        break;
    }

flag = 1;

}

if (flag == 0){
    cerr << "q e p non ammissibili" << endl;
    return false;
}

flag = 0;

//Determinimo la classe della triangolazione
T = b*b + b*c + c*c;

if ((b==0 && c >= 1) || (c==0 && b >= 1)){
    classe = 1;
    flag = 1;
    mesh.V = poligono.F / 2 * T + 2;
    mesh.E = poligono.E * T;
    mesh.F = poligono.F * T;
}
if (b==c && b >= 1){
    classe = 2;
    flag = 1;
    mesh.V = poligono.V + poligono.E * ( 2 * b-1 ) + poligono.F * (( 3 * pow (b,2) - 3*b ) / 2 + 1) ;
    mesh.E = poligono.E * 2 * b + poligono.F * (9 * pow(b,2) + 3 * b ) / 2;
    mesh.F = poligono.E * 2 * b * (b+1);
    
}

if (flag == 0){
    cerr <<"b e c Non sono triangolazioni possibili" << endl;
    return false;
}

poligono.M3D={poligono.V, poligono.E, poligono.F}; // inizializzo M3D con il numero di vertici, spigoli e facce
mesh.M3D={mesh.V, mesh.E, mesh.F}; // inizializzo M3D con il numero di vertici, spigoli e facce
return true;
}

bool ImportMesh(PolygonalMesh& poligono, const string& poligono_platonico)
    {
    if(!ImportCell0Ds(poligono, poligono_platonico))
        return false;
    if(!ImportCell1Ds(poligono, poligono_platonico))
        return false;   
    if(!ImportCell2Ds(poligono, poligono_platonico))
        return false;
    return true;     
    }

bool ImportCell0Ds(PolygonalMesh& poligono, const string& poligono_platonico)
{
    string nome_file = "./PolygonalMesh/Cell0Ds_" + poligono_platonico + ".csv";
    ifstream file_Cell0Ds(nome_file);
    if (!file_Cell0Ds.is_open())
    {
        cerr << "Error opening file "<< endl;
        return false;
    }

    poligono.M0D = MatrixXd::Zero(3, poligono.V); // 3 righe, numero di vertici

    string riga;

    getline (file_Cell0Ds, riga); // leggo l'intestazione

    for ( unsigned int i = 0; i < poligono.V; i++)
    {
        getline (file_Cell0Ds, riga);
        replace(riga.begin(), riga.end(), ';', ' ');
        istringstream convertitore(riga);
        unsigned int id;
        convertitore >> id >> poligono.M0D(0,id) >> poligono.M0D(1,id) >> poligono.M0D(2 ,id);
    }
    
    file_Cell0Ds.close();
    return true;
}

bool ImportCell1Ds(PolygonalMesh& poligono, const string& poligono_platonico){
    
    string nome_file = "./PolygonalMesh/Cell1Ds_" + poligono_platonico + ".csv";
    ifstream file_Cell1Ds(nome_file);
    
    if (!file_Cell1Ds.is_open())
    {
        cerr << "Error opening file " << endl;
        return false;
    }

    poligono.M1D = MatrixXi::Zero(2, poligono.E); // 2 righe, numero di spigoli
    
    string riga;
    getline(file_Cell1Ds, riga); // leggo l'intestazione

    for (unsigned int i = 0; i < poligono.E; i++){
        getline(file_Cell1Ds, riga); // leggo l'intestazione
        replace(riga.begin(), riga.end(), ';', ' ');
        istringstream convertitore(riga);
        unsigned int id;
        convertitore >> id >> poligono.M1D(0,id) >> poligono.M1D(1,id);
    }
    file_Cell1Ds.close();
    return true;
}   

bool ImportCell2Ds(PolygonalMesh& poligono, const string& poligono_platonico){
    
    string nome_file = "./PolygonalMesh/Cell2Ds_" + poligono_platonico + ".csv";
    ifstream file_Cell2Ds(nome_file);
    
    if (!file_Cell2Ds.is_open()){
        cerr << "Error opening file " << endl;
        return false;
    }

    poligono.M2D_vertici.reserve(poligono.F);
    poligono.M2D_spigoli.reserve(poligono.F); // riserva spazio per le facce

    string riga;
    getline(file_Cell2Ds, riga); // leggo l'intestazione

    for (unsigned int i = 0; i < poligono.F; i++){
        
        getline(file_Cell2Ds, riga); // leggo l'intestazione
        replace(riga.begin(), riga.end(), ';', ' ');
        
        istringstream convertitore(riga);
        vector<unsigned int> vertici, spigoli; // vettore per i vertici e spigoli della faccia
        
        vertici.resize(3); // uso resize per allocare la memoria
        spigoli.resize(3); // riservo spazio per 3 vertici
        
        convertitore >> vertici[0] >> vertici[1] >> vertici[2]; // leggo i 3 vertici della faccia
        convertitore >> spigoli[0] >> spigoli[1] >> spigoli[2]; // leggo i 3 vertici della faccia

        poligono.M2D_vertici.push_back(vertici); // aggiungo la linea alla matrice dei vertici
        poligono.M2D_spigoli.push_back(spigoli); // aggiungo la linea alla matrice dei vertici

    }
    file_Cell2Ds.close();
    return true;
}

bool Inizializzazione_vertici( PolygonalMesh& Pmesh , const unsigned int& classe, PolygonalMesh& poligono)
{   
    Pmesh.M0D = MatrixXd::Zero(3, Pmesh.V); // 3 righe, numero di vertici  

    for (unsigned int i = 0; i < poligono.V; i++)
    {
        Pmesh.M0D(0,i) = poligono.M0D(0,i); // X
        Pmesh.M0D(1,i) = poligono.M0D(1,i); // Y
        Pmesh.M0D(2,i) = poligono.M0D(2,i); // Z
    }    
    unsigned int b = Pmesh.d*classe;    
    const auto s = VectorXd::LinSpaced(b+1,0.0,1.0);
    //matrice che conterrà l'ID del j-esimo punto sul i-esimo spigolo.   
    Pmesh.M_pt_spigoli= MatrixXi::Zero(poligono.E,b+1); 

    for (unsigned int i = 0; i < poligono.E ; i++) // lo faccio per ogni spigolo
    {  

        Pmesh.M_pt_spigoli(i,0) = poligono.M1D(0,i); // ID punto iniziale
        Pmesh.M_pt_spigoli(i,b) = poligono.M1D(1,i); // ID punto finale
        Vector3d Acoord = Pmesh.M0D.col(poligono.M1D(0,i)); //estrae coordinate primo estremo da M0D
        Vector3d Bcoord = Pmesh.M0D.col(poligono.M1D(1,i)); //estrae coordinate secondo estremo da M0D
 
        //Creiamo i punti intermedi lungo lo spigolo
        for (unsigned int j = 1; j < b; j++) // ci muoviamo lungo interno saltando il primo e lultimo 
        {               
            unsigned int id = poligono.V +i*(b-1)+j-1; //id punti intermedi  :poligono.V è l'offset per partire dopo i vertici iniziali, i * (b-1) tiene conto dei punti dei precedenti spigoli, j è il punto specifico su questo spigolo
            Pmesh.M_pt_spigoli(i,j) = id; //inserisce id punti intermedi nella matrice degli spigoli TOTALI
            Vector3d coord = Acoord + s(j)* (Bcoord-Acoord); //coordinate punti intermedi
            Pmesh.M0D.col(id) = coord; //inserisce coordinate punti intermedi nella matrice globale
        }

    }	

    unsigned int id_spigoli= 0;
	Pmesh.M1D = MatrixXi::Zero(2,Pmesh.E);
	Pmesh.M1D_spigoli_intermedi= MatrixXi::Zero(poligono.E, b);
	
    for (unsigned int k=0; k < poligono.E; k++)
    {
        for (unsigned int h=0; h < b; h++)
        {
            //inizialmente riempiamo M1D inserendo per ogni id_spigoli gli id dei suoi vertici
            Pmesh.M1D(0,id_spigoli) = Pmesh.M_pt_spigoli(k,h);
            Pmesh.M1D(1,id_spigoli) = Pmesh.M_pt_spigoli(k,h+1);

            //ID dello spigolino h-esimo dello spigolo originale k-esimo.
            //Per ogni spigolo originale k, teniamo traccia dei suoi b spigolini

			Pmesh.M1D_spigoli_intermedi(k,h)=id_spigoli; 
            id_spigoli ++;
        }	 

    }	
    //Dopo aver eseguito tutta la funzione:
    //M0D conterrà i vertici originali + tutti i punti interpolati sugli spigoli
    //M_pt_spigoli ti dirà quali punti sono su ciascuno spigolo (per triangolare le facce)
    return true;
}

Vector3d baricentro(const Vector3d& p0, const Vector3d& p1, const Vector3d& p2){
	Vector3d bar = (p0+p1+p2)/3.0;
	return bar;
}

bool Triangolazione_1_classe(PolygonalMesh& Pmesh, PolygonalMesh& poligono){ //Iteriamo su tutte le facce Pmesh.M2D_vertici
    unsigned int d = Pmesh.d;
    unsigned int num_facce = poligono.F; // numero di facce
    unsigned int n_nuovi_punti = num_facce * (d - 2) * (d - 1) / 2;
    unsigned int id_pt_attuale = Pmesh.V - n_nuovi_punti; 
    unsigned int id_attuale_spigoli = poligono.E*d; // numero di spigoli già creati (dalla triangolazione geodetica)
    Pmesh.M2D_spigoli.reserve(Pmesh.F); 
    Pmesh.M2D_vertici.reserve(Pmesh.F);

    // ciclo sulle facce
    for (unsigned int faccia_id =0;  faccia_id < num_facce;  faccia_id++)
    {   
        // ID dei 3 vertici della faccia triangolare.
        int A_id = poligono.M2D_vertici[faccia_id][0];
        int B_id = poligono.M2D_vertici[faccia_id][1];
        int C_id = poligono.M2D_vertici[faccia_id][2];
       // M2D_spigoli[faccia_id] è un vettore di 3 interi: gli ID degli spigoli. in questo modo trovo l'ID di ogni spigolo per faccia
        int AB_id = poligono.M2D_spigoli[faccia_id][0];
        int BC_id = poligono.M2D_spigoli[faccia_id][1];
        int CA_id = poligono.M2D_spigoli[faccia_id][2];

        //inserisco nel vettore tutti i punti sullo spigolo 
        VectorXi AB = Pmesh.M_pt_spigoli.row(AB_id);
        VectorXi BC = Pmesh.M_pt_spigoli.row(BC_id);
        VectorXi CA = Pmesh.M_pt_spigoli.row(CA_id);
		
        //inserisco nel vettore tutti gli id degli spigolini ( tra i punti nel vettore sopra ) dello spigolo 
		VectorXi AB_spigoli = Pmesh.M1D_spigoli_intermedi.row(AB_id);
		VectorXi BC_spigoli = Pmesh.M1D_spigoli_intermedi.row(BC_id);
		VectorXi CA_spigoli = Pmesh.M1D_spigoli_intermedi.row(CA_id);

        // orientiamo i vertici 

        // vuoi che AB[0] sia A, altrimenti inverti
        if (AB[0] != A_id) {
            VectorXi AB_reversed = AB.reverse().eval();
            AB = AB_reversed;
			VectorXi AB_spigoli_reversed = AB_spigoli.reverse().eval();
            AB_spigoli= AB_spigoli_reversed;
        }

        // vuoi che BC[0] sia B, altrimenti inverti
        if (BC[0] != B_id) {
            VectorXi BC_reversed = BC.reverse().eval();
            BC = BC_reversed;
			VectorXi BC_spigoli_reversed = BC_spigoli.reverse().eval();
            BC_spigoli= BC_spigoli_reversed;
        }

        // vuoi che CA[0] sia C (per costruire da C ad A)
       if (CA[0]==C_id){
            VectorXi CA_reversed = CA.reverse().eval();
            CA = CA_reversed;
			VectorXi CA_spigoli_reversed = CA_spigoli.reverse().eval();
            CA_spigoli= CA_spigoli_reversed;
        }
              
        VectorXi base = AB; // vettore di tutti i punti tra A e B
		VectorXi base_spigoli = AB_spigoli; // vettore di id di spigolini tra A e B
        VectorXi tetto;
		VectorXi tetto_spigoli;

		unsigned int flag_spigoli = 0;
        for (unsigned int i = 0; i < d-1; i ++){

            Pmesh.M1D(0,id_attuale_spigoli)=base[1];
            Pmesh.M1D(1,id_attuale_spigoli)=CA[i+1];
			id_attuale_spigoli++;
                        
            crea_triangolo(Pmesh, base[0],base[1],CA[i+1],base_spigoli[0], id_attuale_spigoli-1,CA_spigoli[i] );

            unsigned int dim = base.size();
            tetto.resize(dim-1);
            tetto[0]=CA[i+1];
			tetto[dim-2]=BC[i+1];
            tetto_spigoli.resize(dim-2);
			const auto s = VectorXd::LinSpaced(d-i,0.0,1.0);// linspaced da 0 a 1 con d-i punti

            // itero sui parallelogrammi dopo aver creato il primo triangolo 
            for (unsigned int j = 1; j < d - i; j ++) 
            {   
				
                if (j == d-i-1){ // ultimo parallelogramma 
					flag_spigoli = 1;
                }
                else{
					Vector3d pt = Pmesh.M0D.col(tetto[0])+s(j)*(Pmesh.M0D.col(tetto[dim-2])-Pmesh.M0D.col(tetto[0]));
                    Pmesh.M0D.col(id_pt_attuale) = pt;
                    tetto[j]=id_pt_attuale;
                    id_pt_attuale++;
                }
				Pmesh.M1D(0,id_attuale_spigoli)=tetto[j-1];
				Pmesh.M1D(1,id_attuale_spigoli)=tetto[j];
				tetto_spigoli[j-1] = id_attuale_spigoli;
                id_attuale_spigoli++;
				
                Pmesh.M1D(0,id_attuale_spigoli)=base[j]; //crea spoigolo obliquo interno
				Pmesh.M1D(1,id_attuale_spigoli)=tetto[j];
				
                crea_triangolo(Pmesh, tetto[j-1], base[j],tetto[j], id_attuale_spigoli-2 ,id_attuale_spigoli, id_attuale_spigoli-1 );
                
                id_attuale_spigoli++;
				
				if (flag_spigoli){
                    crea_triangolo(Pmesh, base[j],base[j+1],tetto[j], base_spigoli[j], BC_spigoli[i] , id_attuale_spigoli-1);
				}
				else{
				    Pmesh.M1D(0,id_attuale_spigoli)=base[j+1]; //crea spoigolo obliquo interno
				    Pmesh.M1D(1,id_attuale_spigoli)=tetto[j];
                    crea_triangolo(Pmesh, base[j],base[j+1],tetto[j], base_spigoli[j], id_attuale_spigoli , id_attuale_spigoli-1);
                    id_attuale_spigoli++;
				}
            }
            // lo riscrivo per poterlo riutilizzare nel ciclo successivo
			flag_spigoli= 0;
            base.resize(dim-1);
            base = tetto;
            base_spigoli.resize(dim-2);
            base_spigoli = tetto_spigoli;
        }
            crea_triangolo(Pmesh, CA[d-1], BC[d-1], CA[d], base_spigoli[0], BC_spigoli[d-1], CA_spigoli[d-1]);
    }

return true;
}

bool Proiezione_sfera(PolygonalMesh& Pmesh)
{
    //Proiezione sferica dei punti
     for (unsigned int i = 0; i < Pmesh.M0D.cols(); i++)
     {
        Vector3d coord = Pmesh.M0D.col(i);
        Pmesh.M0D.col(i) = coord / coord.norm(); // Normalizza le coordinate
    }
    return true;
}

bool stampa_solido(PolygonalMesh& Pmesh,const string& poligono_platonico)
{
    string  nome_file0 = "Cell0Ds" + poligono_platonico + ".txt";
    ofstream file_C0D(nome_file0); 
    if (!file_C0D.is_open())
    {
        cerr << "Error opening file Cell0Ds.txt" << endl;
        return false;
    }
    file_C0D << "ID;X;Y;Z" << endl;
    for (unsigned int i = 0; i < Pmesh.M0D.cols(); i++)
    {
        file_C0D << i << ";" << Pmesh.M0D(0,i) << ";" << Pmesh.M0D(1,i) << ";" << Pmesh.M0D(2,i) << endl;
    }
    file_C0D.close();

    string nome_file1 = "Cell1Ds_" + poligono_platonico + ".txt";
    ofstream file_C1D(nome_file1);
    if (!file_C1D.is_open())
    {
        cerr << "Error opening file Cell1Ds.txt" << endl;
        return false;
    } 
    file_C1D << "ID;Vertice1;Vertice2" << endl;
    for (unsigned int i = 0; i < Pmesh.M1D.cols(); i++)
    {
        file_C1D << i << ";" << Pmesh.M1D(0,i) << ";" << Pmesh.M1D(1,i) << endl;
    }
    file_C1D.close();

    string nome_file2 = "Cell2Ds_" + poligono_platonico + ".txt";
    ofstream file_C2D(nome_file2);
    if (!file_C2D.is_open())
    {
        cerr << "Error opening file Cell2Ds.txt" << endl;
        return false;
    }
    file_C2D << "ID;Vertici;Spigoli" << endl;
    for (unsigned int i = 0; i < size(Pmesh.M2D_vertici) ; i++)
    {
        file_C2D << i << ";";
        for (unsigned int j = 0; j < size(Pmesh.M2D_vertici[i]); j++)
        {
            file_C2D << Pmesh.M2D_vertici[i][j] << ";";
        }
        for (unsigned int j = 0; j < size(Pmesh.M2D_spigoli[i]); j++)
        {
            file_C2D << Pmesh.M2D_spigoli[i][j] << ";";
        }
        file_C2D.seekp(-1, std::ios::end); // vai all'ultimo carattere
        file_C2D.put(' '); // lo sovrascrivi (es. con uno spazio)
        file_C2D << endl; // vai a capo
    }
    file_C2D.close();

    string nome_file3 = "Cell3Ds_" + poligono_platonico + ".txt";
    ofstream file_C3D(nome_file3);
    if (!file_C3D.is_open())
    {
        cerr << "Error opening file Cell3Ds.txt" << endl;
        return false;
    }
    file_C3D << "ID;N^vertici;N^Spigoli;N^Facce" << endl;
    file_C3D << "0;" << Pmesh.V << ";" << Pmesh.E << ";" << Pmesh.F << endl;
    file_C3D.close();

    return true;
}

bool Duale(PolygonalMesh& Pmesh, PolygonalMesh& Dmesh)
{
	Dmesh.M0D = MatrixXd::Zero(3,Pmesh.F);

    Dmesh.M1D = -1* MatrixXi::Ones(2,Pmesh.E); //righe: flag, primo bar, secondo bar colonna id spigolo

    Dmesh.M2D_spigoli.reserve(Pmesh.V); 
    Dmesh.M2D_vertici.reserve(Pmesh.V); 

	//MatrixXi spigoli_baricentri = -1* MatrixXd::ones(2,Pmesh.E); //righe: flag, primo bar, secondo bar colonna id spigolo

	for (unsigned int i=0; i < Pmesh.F; i++){
		Dmesh.M0D.col(i) = baricentro(Pmesh.M0D.col(Pmesh.M2D_vertici[i][0]), Pmesh.M0D.col(Pmesh.M2D_vertici[i][1]), Pmesh.M0D.col(Pmesh.M2D_vertici[i][2]));

        for (unsigned int j=0; j<3; j++){
			if (Dmesh.M1D(0,Pmesh.M2D_spigoli[i][j])== -1)
			{
				Dmesh.M1D(0,Pmesh.M2D_spigoli[i][j]) = i;
			}
			else{
				Dmesh.M1D(1,Pmesh.M2D_spigoli[i][j]) = i;
			}
		}
	}
	
    MatrixXi id_spigoli = -1*MatrixXi::Ones(10,Pmesh.F);
	
	unsigned int id_vertice;
	unsigned int id_spigolo_precedente;
	unsigned int id_spigolo_successivo;
	unsigned int id_vertice_iniziale;
	
	for (unsigned int k =0; k < Pmesh.E; k++){

		unsigned int A_id = Pmesh.M1D(0,k);// k che è l'id dello spigolo (e cpincide con la i di prima)  richiama gli id dei due vertici delo spigolo
		unsigned int B_id = Pmesh.M1D(1,k);
		
		for (unsigned int j=0; j<10; j++){
			if (id_spigoli(j,A_id) == -1){
				id_spigoli(j,A_id) = k;
				break;
			}
		}
		for (unsigned int j=0; j<10; j++){
			if (id_spigoli(j,B_id) == -1){
				id_spigoli(j,B_id) = k;
				break;
			}
		}
	}

	for (unsigned int h = 0; h < Pmesh.V; h++){
        
        vector<unsigned int> vertici;
        vector<unsigned int> spigoli;

        vertici.reserve(10);
        spigoli.reserve(10);

		id_spigolo_precedente = id_spigoli(0,h);
        
        spigoli.push_back(id_spigolo_precedente);
		
		id_vertice_iniziale = Dmesh.M1D(0,id_spigolo_precedente);
        
        vertici.push_back(id_vertice_iniziale);

		id_vertice = Dmesh.M1D(1,id_spigolo_precedente);
		
		while (id_vertice != id_vertice_iniziale){
			for(unsigned int j=1; j<10; j++){
				id_spigolo_successivo = id_spigoli(j,h);
				
				if (id_spigolo_successivo != id_spigolo_precedente){
					
					unsigned int A_id = Dmesh.M1D(0,id_spigolo_successivo);
					unsigned int B_id = Dmesh.M1D(1,id_spigolo_successivo);
					
					if (A_id == id_vertice){
                        vertici.push_back(A_id);
                        spigoli.push_back( id_spigolo_successivo);
						id_vertice = B_id;
						id_spigolo_precedente = id_spigolo_successivo;
						break;
					}
					
					if (B_id == id_vertice){
                        vertici.push_back(B_id);
                        spigoli.push_back( id_spigolo_successivo);
						id_vertice = A_id;
						id_spigolo_precedente = id_spigolo_successivo;
						break;
					}
				}
			}
		}
    Dmesh.M2D_vertici.push_back(vertici);
    Dmesh.M2D_spigoli.push_back(spigoli);
  
	}
    Dmesh.V = Pmesh.F;
    Dmesh.E = Pmesh.E;
    Dmesh.F = Pmesh.V; 
    Dmesh.M3D = {Dmesh.V, Dmesh.E, Dmesh.F}; // inizializzo M3D con il numero di vertici, spigoli e facce
    Dmesh.nodo_i = Pmesh.nodo_i;
    Dmesh.nodo_f = Pmesh.nodo_f;
	return true;
}

bool crea_triangolo(PolygonalMesh& Pmesh, const unsigned int& id_pt_1, const unsigned int& id_pt_2, const unsigned int& id_pt_3, const unsigned int& id_sp_1, const unsigned int& id_sp_2, const unsigned int& id_sp_3){

    Pmesh.M2D_vertici.push_back({id_pt_1,id_pt_2,id_pt_3});
    Pmesh.M2D_spigoli.push_back({id_sp_1,id_sp_2,id_sp_3}); 
    
    return true;
}

bool Triangolazione_2_classe(PolygonalMesh& Pmesh, PolygonalMesh& poligono){
	unsigned int d = Pmesh.d;
    unsigned int num_facce = poligono.F; // numero di facce
    unsigned int n_nuovi_punti = num_facce *  (pow(d,2)+((d-2)*(d-1))/2);
    unsigned int id_pt_attuale = Pmesh.V - n_nuovi_punti;
    unsigned int id_attuale_spigolo = poligono.E*2*d; // numero di spigoli già creati (dalla triangolazione geodetica) + quelli che creeremo ora
    
    vector<int> id_baricentri_CA ;// max numero elementi necessari in un ciclo riga

    Pmesh.M2D_spigoli.reserve(Pmesh.F);
    Pmesh.M2D_vertici.reserve(Pmesh.F); 

	
	for (unsigned int faccia_id =0;  faccia_id < num_facce;  faccia_id++)
    {   
        // ID dei 3 vertici della faccia triangolare.
        int A_id = poligono.M2D_vertici[faccia_id][0];
        int B_id = poligono.M2D_vertici[faccia_id][1];
        int C_id = poligono.M2D_vertici[faccia_id][2];

       //M2D_spigoli[faccia_id] è un vettore di 3 interi: gli ID degli spigoli. in questo modo trovo l'ID di ogni spigolo per faccia
        int AB_id = poligono.M2D_spigoli[faccia_id][0];
        int BC_id = poligono.M2D_spigoli[faccia_id][1];
        int CA_id = poligono.M2D_spigoli[faccia_id][2];

        //inserisco nel vettore tutti i punti sullo spigolo 
        VectorXi AB = Pmesh.M_pt_spigoli.row(AB_id);
        VectorXi BC = Pmesh.M_pt_spigoli.row(BC_id);
        VectorXi CA = Pmesh.M_pt_spigoli.row(CA_id);
		
        //inserisco nel vettore tutti gli id degli spigolini ( tra i punti nel vettore sopra ) dello spigolo 
		VectorXi AB_spigoli = Pmesh.M1D_spigoli_intermedi.row(AB_id);
		VectorXi BC_spigoli = Pmesh.M1D_spigoli_intermedi.row(BC_id);
		VectorXi CA_spigoli = Pmesh.M1D_spigoli_intermedi.row(CA_id);

        // orientiamo i vertici 

        // vuoi che AB[0] sia A, altrimenti inverti
        if (AB[0] != A_id) {
            VectorXi AB_reversed = AB.reverse().eval();
            AB = AB_reversed;
			VectorXi AB_spigoli_reversed = AB_spigoli.reverse().eval();
            AB_spigoli= AB_spigoli_reversed;
        }

        // vuoi che BC[0] sia B, altrimenti inverti
        if (BC[0] != B_id) {
            VectorXi BC_reversed = BC.reverse().eval();
            BC = BC_reversed;
			VectorXi BC_spigoli_reversed = BC_spigoli.reverse().eval();
            BC_spigoli= BC_spigoli_reversed;
        }

        // vuoi che CA[0] sia C (per costruire da C ad A)
       if (CA[0]==C_id){
            VectorXi CA_reversed = CA.reverse().eval();
            CA = CA_reversed;
			VectorXi CA_spigoli_reversed = CA_spigoli.reverse().eval();
            CA_spigoli= CA_spigoli_reversed;
        }

        //Controlla che Pmesh.M2D_spigoli[faccia_id] abbia ID coerenti per ogni faccia        
        VectorXi base = AB; // vettore di tutti i punti tra A e B
		VectorXi base_spigoli = AB_spigoli;
        VectorXi tetto;
		VectorXi tetto_spigoli;
		VectorXd s;

        //itero sui livelli del triangolo 
		for (unsigned int i = 0; i < d; i ++){

            tetto = VectorXi::Zero(d*2-2*i-1);
            tetto_spigoli = VectorXi::Zero(tetto.size()-1);

            tetto[0]= CA[2*i+2];
            tetto[tetto.size()-1]=BC[2*i+2];
            
            //linspace sul tetto di i 
			s = VectorXd::LinSpaced(d-i, 0.0, 1.0);

            Vector3d coord_iniziale = Pmesh.M0D.col(CA[(i+1)*2]);
            Vector3d coord_finale = Pmesh.M0D.col(BC[(i+1)*2]);

			for (unsigned int l=1; l < d-i-1; l++)
            {// iteriamo sulla lunghezza se = 2 non lo fa
                Pmesh.M0D.col(id_pt_attuale)= coord_iniziale +s(l)*(coord_finale-coord_iniziale);
				tetto[2*l] = id_pt_attuale; // tetto continene id dei punti sul tetto di i 
                id_pt_attuale++;
			}

            unsigned int obliquo = CA[2*i+1]; 
            unsigned int obliquo_spigolo1 = CA_spigoli[i*2];
            unsigned int obliquo_spigolo2 = CA_spigoli[i*2+1];

            for (unsigned int k=0; k < d-i-1; k++)
             {                     
                Pmesh.M0D.col(id_pt_attuale) = baricentro(Pmesh.M0D.col(base[2*k]), Pmesh.M0D.col(base[2*k+2]), Pmesh.M0D.col(tetto[k*2]));
                //unsigned int baricentro1 = id_pt_attuale;corrisponde a id_pt_attuale-2);
                id_pt_attuale++;
                Pmesh.M0D.col(id_pt_attuale) = baricentro(Pmesh.M0D.col(tetto[k*2+2]), Pmesh.M0D.col(base[2*k+2]), Pmesh.M0D.col(tetto[k*2]));
                tetto[2*k+1] = id_pt_attuale; //baricentro che serve per il parallelogramma successivo
                id_pt_attuale++;
                // id spigoli nel tetto del parallelogramma k 
                Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(tetto[2*k], tetto[2*k+1]); //id_punto_attuale - 1  è il baricentro precedente
                tetto_spigoli[2*k] = id_attuale_spigolo;
                id_attuale_spigolo++;     
                Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(tetto[2*k+1], tetto[2*k+2]);
                tetto_spigoli[2*k+1] = id_attuale_spigolo;
                id_attuale_spigolo++;
                Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(base[2*k], id_pt_attuale-2);
                // unsigned int uno = id_attuale_spigolo;
                id_attuale_spigolo++;
                Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(base[2*k+1], id_pt_attuale-2);
                //unsigned int due = id_attuale_spigolo;
                id_attuale_spigolo++;
                Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(base[2*k+2], id_pt_attuale-2);
                //unsigned int tre = id_attuale_spigolo;
                id_attuale_spigolo++;
                Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(obliquo, id_pt_attuale-2);
                //unsigned int quattro = id_attuale_spigolo;
                id_attuale_spigolo++;
                Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(tetto[2*k], id_pt_attuale-2);
                //unsigned int cinque = id_attuale_spigolo;
                id_attuale_spigolo++;
                Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(tetto[2*k+1], id_pt_attuale-2);
                //unsigned int sei = id_attuale_spigolo;
                id_attuale_spigolo++;
                Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(base[2*k+2], tetto[2*k+1]);
                //unsigned int sette = id_attuale_spigolo;
                id_attuale_spigolo++;       
                crea_triangolo(Pmesh, base[2*k], base[2*k+1], id_pt_attuale-2,base_spigoli[2*k], id_attuale_spigolo-6, id_attuale_spigolo-7); 
                //crea_triangolo(Pmesh, id_attuale_triangoli, base[2*k], base[2*k+1], baricentro1, uno, base_spigoli[2*k], due); 
                //id_attuale_triangoli++;              
                crea_triangolo(Pmesh, base[2*k+1], base[2*k+2], id_pt_attuale-2, base_spigoli[2*k+1], id_attuale_spigolo-5, id_attuale_spigolo-6);
                //crea_triangolo(Pmesh, id_attuale_triangoli, base[2*k+1], base[2*k+2], baricentro1, base_spigoli[2*k+1], tre, due);
                //id_attuale_triangoli++;              
                crea_triangolo(Pmesh, base[2*k], id_pt_attuale-2, obliquo, id_attuale_spigolo-7, id_attuale_spigolo-4, obliquo_spigolo1);
                //crea_triangolo(Pmesh, id_attuale_triangoli, base[2*k], baricentro1, obliquo, uno, quattro, obliquo_spigolo1);
                //id_attuale_triangoli++;              
                crea_triangolo(Pmesh, obliquo, id_pt_attuale-2, tetto[2*k], id_attuale_spigolo-4, id_attuale_spigolo-3, obliquo_spigolo2);
                //crea_triangolo(Pmesh, id_attuale_triangoli, obliquo, baricentro1, tetto[2*k], quattro, cinque, obliquo_spigolo2);
                //id_attuale_triangoli++;               
                crea_triangolo(Pmesh, id_pt_attuale-2, tetto[2*k], tetto[2*k+1], id_attuale_spigolo-3, tetto_spigoli[2*k] , id_attuale_spigolo-2);
                //crea_triangolo(Pmesh, id_attuale_triangoli, baricentro1, tetto[2*k], tetto[2*k+1], cinque, tetto_spigoli[2*k] , sei);
                //id_attuale_triangoli++;
                crea_triangolo(Pmesh, base[2*k+2], id_pt_attuale-2, tetto[2*k+1], id_attuale_spigolo-5, id_attuale_spigolo-2, id_attuale_spigolo-1);
                //crea_triangolo(Pmesh, id_attuale_triangoli, base[2*k+2], baricentro1, tetto[2*k+1], tre, sei, sette);
                //id_attuale_triangoli++;
                
                obliquo = tetto[2*k+1]; // il baricentro del parallelogramma
                obliquo_spigolo1 = id_attuale_spigolo-1; // id del primo spigolo del parallelogramma
                // obliquo_spigolo1 = sette;
                obliquo_spigolo2 = tetto_spigoli[2*k+1]; // id del secondo spigolo del parallelogramma
            }                    
            // creo il triangolo finale del livello
            Pmesh.M0D.col(id_pt_attuale)=baricentro(Pmesh.M0D.col(tetto[tetto.size()-1]), Pmesh.M0D.col(base[base.size()-1]), Pmesh.M0D.col(base[base.size()-3]));
            //unsigned int ultimo_baricentro = id_pt_attuale;      
            Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(obliquo, id_pt_attuale); 
            //unsigned int alfa = id_attuale_spigolo;
            id_attuale_spigolo++;
            Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(tetto[tetto.size()-1], id_pt_attuale); 
            //unsigned int beta = id_attuale_spigolo;
            id_attuale_spigolo++;
            Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(base[base.size()-3], id_pt_attuale); 
            //unsigned int gamma = id_attuale_spigolo;
            id_attuale_spigolo++;
            Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(base[base.size()-2], id_pt_attuale); 
            //unsigned int delta = id_attuale_spigolo;
            id_attuale_spigolo++;
            Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(base[base.size()-1], id_pt_attuale); 
            //unsigned int epsilon = id_attuale_spigolo;
            id_attuale_spigolo++;    
            Pmesh.M1D.col(id_attuale_spigolo) = Vector2i(BC[2*i+1], id_pt_attuale); 
            //unsigned int eta = id_attuale_spigolo;
            id_attuale_spigolo++; 
            crea_triangolo(Pmesh, base[base.size()-3], base[base.size()-2], id_pt_attuale, base_spigoli[base_spigoli.size()-2], id_attuale_spigolo-3, id_attuale_spigolo-4);
            //crea_triangolo(Pmesh, id_attuale_triangoli, base[base.size()-3], base[base.size()-2], ultimo_baricentro, base_spigoli[base_spigoli.size()-2], delta, gamma);
            //id_attuale_triangoli++;
            crea_triangolo(Pmesh, base[base.size()-2], base[base.size()-1], id_pt_attuale, base_spigoli[base_spigoli.size()-1], id_attuale_spigolo-2, id_attuale_spigolo-3);
            //crea_triangolo(Pmesh, id_attuale_triangoli, base[base.size()-2], base[base.size()-1], ultimo_baricentro, base_spigoli[base_spigoli.size()-1], epsilon, delta);
            //id_attuale_triangoli++;
            crea_triangolo(Pmesh, base[base.size()-3], id_pt_attuale, obliquo, id_attuale_spigolo-4 , id_attuale_spigolo-6, obliquo_spigolo1);
            //crea_triangolo(Pmesh, id_attuale_triangoli, base[base.size()-3], ultimo_baricentro, obliquo, gamma, alfa, obliquo_spigolo1);
            //id_attuale_triangoli++;
            crea_triangolo(Pmesh, obliquo, id_pt_attuale, tetto[tetto.size()-1], id_attuale_spigolo-6, id_attuale_spigolo-5 , obliquo_spigolo2);
            //crea_triangolo(Pmesh, id_attuale_triangoli, obliquo, ultimo_baricentro, tetto[tetto.size()-1], alfa, beta, obliquo_spigolo2);
            //id_attuale_triangoli++;
            crea_triangolo(Pmesh, id_pt_attuale, tetto[tetto.size()-1], BC[2*i+1], id_attuale_spigolo-5 ,BC_spigoli[i*2+1], id_attuale_spigolo-1);
            //crea_triangolo(Pmesh, id_attuale_triangoli, ultimo_baricentro, tetto[tetto.size()-1], BC[2*i+1], beta, eta, BC_spigoli[i*2+1]);
            //id_attuale_triangoli++;
            crea_triangolo(Pmesh, id_pt_attuale, BC[2*i+1], base[base.size()-1], id_attuale_spigolo-1 , BC_spigoli[2*i], id_attuale_spigolo-2);
            //crea_triangolo(Pmesh, id_attuale_triangoli, ultimo_baricentro, BC[2*i+1], base[base.size()-1], eta, BC_spigoli[2*i], epsilon);
            //id_attuale_triangoli++;
            id_pt_attuale++;

            base.resize(tetto.size());
            base = tetto;
            base_spigoli.resize(tetto_spigoli.size());
            base_spigoli = tetto_spigoli;
        }
}
return true;
}

bool CamminoMinimo(PolygonalMesh& Pmesh){

    if (Pmesh.nodo_i >= Pmesh.V || Pmesh.nodo_f >= Pmesh.V || Pmesh.nodo_f == Pmesh.nodo_i){
        cerr << "gli id non sono validi, controlla che non siano uguali oppure inserisci numeri piu piccoli di: " << Pmesh.V << endl;
        return false;
    }
    //Associo ad ogni nodo l'id degli spigoli associati 
    vector<vector<unsigned int>> Vettore_spigoli_del_nodo(Pmesh.V);
    for (vector<unsigned int>& numero_spigoli_per_nodo : Vettore_spigoli_del_nodo)
        numero_spigoli_per_nodo.reserve(10);

    for (unsigned int id_spigolo = 0; id_spigolo < Pmesh.E; id_spigolo++) {
        unsigned int id_nodo_A = Pmesh.M1D(0, id_spigolo);
        unsigned int id_nodo_B = Pmesh.M1D(1, id_spigolo);
        Vettore_spigoli_del_nodo[id_nodo_A].push_back(id_spigolo); // Mette al vettore del nodo A l'id dello spigolo
        Vettore_spigoli_del_nodo[id_nodo_B].push_back(id_spigolo);
    }

    //ALGORITMO DI DIJKSTRA per visitare il grafo

    vector<unsigned int> pred(Pmesh.V,-1); //inizilizzo il vettore predecessore a -1
    vector<double> dist(Pmesh.V, numeric_limits<double>::infinity()); //inizilizzo il vettore delle distanze a un numero molto grande
    vector<bool> visitato(Pmesh.V,false);
    using ND = pair<double, int>; //Ogni elemento della coda è una coppia formata da: int (ID di un nodo) e double (distanza)
    priority_queue < ND, vector<ND>, greater<ND>> PQ; // priority_queue<  Tipo,   Contenitore,     Criterio >: classe STL per le code con priorità (distanza)
    // greater(criterio) per avere che il minimo sia estratto per primo ( sia messo in cima ed estratto con top )
    // vector(contenitore) per supportare accesso casuale
    
    pred[Pmesh.nodo_i] = Pmesh.nodo_i;
    dist[Pmesh.nodo_i] = 0.0;
    PQ.push({0.0, Pmesh.nodo_i});    // for(unsigned int i = 0; i < LA.size(); i++){    //     PQ.push({dist[i], i}); Non è necessario
    
    unsigned int u = Pmesh.nodo_i; // inizializzo u a un valore che sicuro non potrà assumere il nodo finale(impostato all'inizio del programma)
    unsigned int contatore = 0;
    // Se u diventa nodo_f significa che è stato prelevato da PQ, di conseguenza è la piu piccola distanza in PQ, quindi d[u] sarà proprio uguale alla distanza minima. 
    while(u != Pmesh.nodo_f){ 
        // contatore per calcolare il reserve sul vettore percorso creato in seguito
        contatore += 1;
        // double p = PQ.top().first; Distanza minimima di u. Serviva per   // if ( p > dist[u] ) continue; si è optato per una versione più ottimale per grafi con moltinodi 
        u = PQ.top().second; // punto
        PQ.pop(); //rimuove ( non restituisce nulla ) l’elemento con priorità più alta 
        // scarta entry obsolete: se esiste già una distanza minore trovata per u, ignoro questa versione ( problema generato dal push in PQ)
        if ( visitato[u] ) continue;        // if ( p > dist[u] ) continue;
        visitato[u] = true;

        // ottimizzazione: se dist[u] è già uguale a dist[nodo_f], allora nodo_f è già stato raggiunto al minimo. Quindi si interrompe il ciclo
        //if ( u != Pmesh.nodo_f && dist[u] == dist[Pmesh.nodo_f] && dist[Pmesh.nodo_f]!= numeric_limits<double>::infinity()) break; // se la distanza associata al nodo u è uguale alla distanza del nodo finale nella lista delle distanze allora per forza è gia la distanza minima per arrivare al nodo finale
        
        // sto iterando su tutti gli spigoli condivisi da u, che mi permettono di arrivare ai nodi vicini.
        for (unsigned int id_spigolo_associato_nodo_u : Vettore_spigoli_del_nodo[u]) {
            unsigned int id_nodoA = Pmesh.M1D(0, id_spigolo_associato_nodo_u);
            unsigned int id_nodoB = Pmesh.M1D(1, id_spigolo_associato_nodo_u);
            //Iterando sugli spigoli costruiti con u, quindi uno dei due nodi che lo crea è proprio u. Inizializzo quindi v come l'altro nodo
            unsigned int v = (u == id_nodoA ? id_nodoB : id_nodoA); 
            //visitato[v]=True se v è stata prelevata da PQ, dunque si conosce d[v](distanza min tra v e nodo_i). Non ha senso quindi confrontarla rispetto a u.
            if ( visitato[v] ) continue;
            //distanza tra Nodo u e nodo v. La norma è O(1), ma sqrt diventa costosa se ripetuta più volte. Per evitare che venga calcolata in casi superflui si fa il controllo precedente. 
            double w = (Pmesh.M0D.col(u) - Pmesh.M0D.col(v)).norm();
            //dist[u] è la distanza più piccola che per ora l'algoritmo è riuscita a calcolare tra il nodo_i e u. 
            double nuova_distanza = dist[u] + w;
            //Controlliamo se abbia più senso raggiungere v passando per u. 
            if (nuova_distanza < dist[v]) {
                // se ha senso allora aggiorno la distanza minima per raggiungerer v con la distanza per arrivare a u + la distanza tra u e v.
                dist[v] = nuova_distanza;
                // Si stabilisce qindi che il nodo precedente per giungere v sarà u.
                pred[v] = u;
                // "decresekey".Non propriamente. Non esistendo un comando in c++ che lo faccia. Si inserisce un elemento in PQ associato al nodo a w. Il risultato è che si avranno piu elementi associati allo stesso nodo w.
                PQ.push({nuova_distanza, v});
            }

        }
    }

    // path contiene id dei vertici del cammino minimo 
    vector<unsigned int> percorso;
    percorso.reserve(contatore);
	Pmesh.percorso.reserve(contatore);
	
    unsigned int nodo_ritroso = Pmesh.nodo_f;
    // riempe al contrario il percodo
    while(nodo_ritroso != Pmesh.nodo_i){
        percorso.push_back(nodo_ritroso);
        nodo_ritroso = pred[nodo_ritroso];
    } 
    percorso.push_back(Pmesh.nodo_i);
	Pmesh.percorso = percorso;

    cout << "la lunghezza del porcorso minimo è di: " << dist[Pmesh.nodo_f] << endl << "Il numero di nodi da percorrere è: " << percorso.size() << endl << "Il percorso minimo è dato da: " << endl;

    for ( unsigned int i = 1; i < percorso.size() + 1 ; i++) // reverse iterator
        cout << percorso[percorso.size()-i] << " ";
    cout << endl ;
    // coloriamo i punti relativi al percorso
    vector<double> ProprietaPunti_path(Pmesh.V, 0.0);

    for(unsigned int punto : percorso){
        ProprietaPunti_path[punto] = 1.0;
    }

    //info utili per utilizzare le funzioni di Vicini
    Gedim::UCDProperty<double> ShortPathProperty;

        ShortPathProperty.Label = "cammino minimo ";
        ShortPathProperty.UnitLabel = "";
        ShortPathProperty.Size = ProprietaPunti_path.size();
        ShortPathProperty.NumComponents = 1;
        ShortPathProperty.Data = ProprietaPunti_path.data();  

    vector<Gedim::UCDProperty<double>> ProprietaPunti;

    ProprietaPunti.push_back( ShortPathProperty );

    Gedim::UCDUtilities utilities;
    utilities.ExportPoints( "./Percorso_minimo_nodi.inp",
                            Pmesh.M0D,
                            ProprietaPunti);

    // coloriamo i lati relativi al percorso e calcoliamo la lunghezza del percorso
    vector<unsigned int> segmenti_Percorso;
    //segmenti_Percorso.reserve(path.size()-1);

    vector<double> ProprietaLatiPercorso(Pmesh.E, 0.0 );


    //itero su ogni coppia di punti adiacenti del percorso 
    for (unsigned int i = 0; i < percorso.size() - 1; i++) {

    unsigned int pt_a = percorso[i];
    unsigned int pt_b = percorso[i+1];

    // cerco lo spigolo che collega pt_a e pt_b
    //itero su ogni id_spigolo

    for (unsigned int j = 0; j<  Pmesh.E; j ++ )
    {
        //id dei punti estremi del segmento
        unsigned int id_A = Pmesh.M1D(0, j);
        unsigned int id_B = Pmesh.M1D(1, j);

        if ((id_A==pt_a && id_B==pt_b)||(id_A==pt_b && id_B==pt_a)) {

            segmenti_Percorso.push_back(j); //inserisco id dello spigolo
            ProprietaLatiPercorso[j] = 1.0;
            break;
        }
    }
    } 

    ShortPathProperty.Label = "Percorso minimo";
    ShortPathProperty.Size = ProprietaLatiPercorso.size();
    ShortPathProperty.NumComponents = 1;
    ShortPathProperty.Data = ProprietaLatiPercorso.data();  

    vector<Gedim::UCDProperty<double>> Segmenti_Proprietà;

    Segmenti_Proprietà.push_back(ShortPathProperty);

    utilities.ExportSegments("./Percorso_minimo_archi.inp",
                            Pmesh.M0D,
                            Pmesh.M1D,
                            {},
                            Segmenti_Proprietà);

    return true;

}


}

// bool CamminoMinimo2(PolygonalMesh& Pmesh){
    
//     unsigned int  id_vertice_iniziale = Pmesh.nodo_i;
//     unsigned int  id_vertice_finale = Pmesh.nodo_f;
     
//      if (id_vertice_iniziale >= Pmesh.V || id_vertice_iniziale < 0 || id_vertice_finale >= Pmesh.V || id_vertice_finale < 0 ){
//         cerr << "gli id non sono validi" << endl;
//         return false;
//      }
 
//      //creiamo la LISTA DI ADIACENZA
//      vector<vector<unsigned int>> LA;
 
//      LA.reserve(Pmesh.V);
 
//      MatrixXd W = MatrixXd::Zero(Pmesh.V, Pmesh.V);
     
//      //riempo la matrice di adiacenza
 
//     for (size_t j = 0 ; j < Pmesh.V ; j ++ ){
 
//          vector<unsigned int> punti_vicini;
//          punti_vicini.reserve(10);
 
//          //considero ogni spigolo e ne prendo i punti 
//          for (size_t i = 0 ; i < Pmesh.E ; i ++){
 
//              unsigned int origine = Pmesh.M1D_triangolini(0,i);
//              unsigned int fine = Pmesh.M1D_triangolini(1,i);
 
//          // se origine del segmento che considero ha id = i , inserisco la fine del segmento nei punti vicini
//          // se fine del segmento che considero ha id = i , inserisco l'inizio del segmento nei punti vicini
 
//              if (origine == j)
//                  punti_vicini.push_back(fine);
 
//              else if (fine == j)
//                  punti_vicini.push_back(origine);   
 
//          }
         
         
//          //per ogni punto inserisco nel vector di vector i suoi vicini 
         
//          LA.push_back(punti_vicini);
//      }
//      //itero sui punti per ottenere matrice di adiacenza 

//     for ( size_t i=0; i< LA.size(); i++){
       
//         //itero sul vettore di punti vicini
//         for (size_t k = 0; k < LA[i].size(); k++) {
//         unsigned int w = LA[i][k];  // vicino del punto i 
//         W(i, w) = (Pmesh.M0D.col(w) - Pmesh.M0D.col(i)).norm();
//         }
//     }

   

    
//     //ALGORITMO DI DIJKSTRA per visitare il grafo

//     //inizilizzo il vettore predecessore a -1
//     vector<unsigned int> pred(Pmesh.V,-1); // inizilizzo il predecessore a -1
 
//     //inizilizzo il vettore delle distanze a un numero molto grande
// 	vector<double> dist(Pmesh.V, numeric_limits<double>::infinity());

//     //Ogni elemento della coda è una coppia formata da: int (ID di un nodo) e double (distanza)

//     using ND = pair<double, int>; //coppia nodo-distanza !!!!!!!CONTROLLARE ORDINE

// 	// priority_queue<  Tipo,   Contenitore,     Criterio >: classe STL per le code con priorità (distanza)
//     // greater(criterio) per avere che il minimo sia estratto per primo ( sia messo in cima ed estratto con top )
//     // vector(contenitore) per supportare accesso casuale

//     priority_queue< ND, vector<ND>, greater<ND>> PQ;
		
//     pred[id_vertice_iniziale] = id_vertice_iniziale;
//     dist[id_vertice_iniziale] = 0.0;

//     for(int i = 0; i < LA.size(); i++){
//         PQ.push({dist[i], i}); 
//     }
//     while(!PQ.empty()){

//         int p = PQ.top().first; // distanza
//         int u = PQ.top().second; // punto

//         PQ.pop(); //rimuove ( non restituisce nulla ) l’elemento con priorità più alta 

//         // Se la distanza più aggiornata per u (dist[u]) è già più piccola di quella che sto estraendo (p) ,  non la uso.
//        // if (dist[u]<p)//
//       //      continue;//
        
//         for(const auto& w : LA[u]){ //foreach
//             if( dist[w] > dist[u] + W(u,w) ) {
//                 //Aggiorna distanza se hai trovato un cammino più corto
//                 dist[w] = dist[u] + W(u,w); 
//                 pred[w] = u;

//                 PQ.push({dist[w],w});
//             }
//         }
// 	}

//     // path contiene id dei vertici del cammino minimo 
//     vector<unsigned int> path;

//     path.reserve(Pmesh.V);

//     int v = id_vertice_finale;
    
//     //esegiuiamo il cammino al contrario cercando sempre il predecessore

//     while(v != id_vertice_iniziale) {
//         path.push_back(v);
//         v = pred[v];
// 	} 
    
//     path.push_back(id_vertice_iniziale);

//     cout << "la lunghezza del porcorso minimo è di: " << dist[Pmesh.nodo_f] << endl << "Il numero di nodi da percorrere è: " << path.size() << endl << "Il percorso minimo è dato da: " << endl;

//     for ( unsigned int i = 1; i < path.size() + 1 ; i++) // reverse iterator
//         cout << path[path.size()-i] << " ";
//     cout << endl ;
    
		
//     // coloriamo i punti relativi al percorso
//     vector<double> ProprietaPunti_path(Pmesh.V, 0.0);

// 	for(unsigned int punto : path){
// 		ProprietaPunti_path[punto] = 1.0;
//     }

//     //info utili per utilizzare le funzioni di Vicini
//     Gedim::UCDProperty<double> ShortPathProperty;

// 			ShortPathProperty.Label = "cammino minimo ";
// 			ShortPathProperty.UnitLabel = "";
// 			ShortPathProperty.Size = ProprietaPunti_path.size();
// 			ShortPathProperty.NumComponents = 1;
// 			ShortPathProperty.Data = ProprietaPunti_path.data();  
//             vector<Gedim::UCDProperty<double>> ProprietaPunti;

//             ProprietaPunti.push_back( ShortPathProperty );
        
//             Gedim::UCDUtilities utilities;
//             utilities.ExportPoints("./prova2_nodi",
//                                         Pmesh.M0D,
//                                     ProprietaPunti);
        
        
//             // coloriamo i lati relativi al percorso e calcoliamo la lunghezza del percorso
//             vector<unsigned int> segmenti_Percorso;
//             //segmenti_Percorso.reserve(path.size()-1);
        
//             vector<double> ProprietaLatiPercorso(Pmesh.E, 0.0 );
        
//             double lunghezzaPercorso = 0.0;
        
//             //itero su ogni coppia di punti adiacenti del percorso 
//             for (unsigned int i = 0; i < path.size() - 1; i++) {
        
//                 unsigned int pt_a = path[i];
//                 unsigned int pt_b = path[i + 1];
        
//                 lunghezzaPercorso += W(pt_a , pt_b);
        
//                 // cerco lo spigolo che collega pt_a e pt_b
//                 //itero su ogni id_spigolo
        
//                 for (unsigned int j = 0; j<  Pmesh.E; j ++ )
//                 {
//                     //id dei punti estremi del segmento
//                     unsigned int id_A = Pmesh.M1D_triangolini(0, j);
//                     unsigned int id_B = Pmesh.M1D_triangolini(1, j);
        
//                     if ((id_A==pt_a && id_B==pt_b)||(id_A==pt_b && id_B==pt_a)) {
        
//                         segmenti_Percorso.push_back(j); //inserisco id dello spigolo
//                         ProprietaLatiPercorso[j] = 1.0;
//                     }
//                 }
//             } 
//             ShortPathProperty.Label = "Percorso minimo";
//             ShortPathProperty.Size = ProprietaLatiPercorso.size();
//             ShortPathProperty.NumComponents = 1;
//             ShortPathProperty.Data = ProprietaLatiPercorso.data();  
        
//             vector<Gedim::UCDProperty<double>> Segmenti_Proprietà;
        
//             Segmenti_Proprietà.push_back(ShortPathProperty);
        
//             utilities.ExportSegments("./prova2",
//                                           Pmesh.M0D,
//                                           Pmesh.M1D_triangolini,
//                                           {},
//                                           Segmenti_Proprietà);
        
//          return true;
        
//         }
        
        


// }
