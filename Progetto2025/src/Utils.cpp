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
bool input_solido_platonico(PolygonalMesh& mesh, int argc ,char* argv[]){

if (argc !=5)
    {
        cout << "Error: the number of input is not correct" << endl;
        return false;
    }

    int p, q, b, c, V, E, F, T;  
    int flag = 0;
    istringstream convert0(argv[1]);
    convert0 >> p; //numero di lati per faccia
    istringstream convert1(argv[2]);
    convert1 >> q; //numero di facce che si incontrano in ogni vertice
    istringstream convert2(argv[3]);
    convert2 >> b; //parametro triangolazione geodetica
    istringstream convert3(argv[4]);
    convert3 >> c; //parametro triangolazione geodetica

    mesh.p = p; // riempiamo il contenitore mesh con i dati di input
    mesh.q = q;

if (p==3){ //facce triangolari
    switch (q)
    {
    case 3: //tetraedro
        mesh.nomefile0 = "./PolygonalMesh/Cell0Ds_tetraedro.csv";
        mesh.nomefile1 = "./PolygonalMesh/Cell1Ds_tetraedro.csv"; 
        mesh.nomefile2 = "./PolygonalMesh/Cell2Ds_tetraedro.csv";
        break;
    case 4: //ottaedro 
        mesh.nomefile0 = "./PolygonalMesh/Cell0Ds_ottaedro.csv";
        mesh.nomefile1 = "./PolygonalMesh/Cell1Ds_ottaedro.csv";
        mesh.nomefile2 = "./PolygonalMesh/Cell2Ds_ottaedro.csv";
        break;
    case 5: //icosaedro
        mesh.nomefile0 = "./PolygonalMesh/Cell0Ds_icosaedro.csv";
        mesh.nomefile1 = "./PolygonalMesh/Cell1Ds_icosaedro.csv";
        mesh.nomefile2 = "./PolygonalMesh/Cell2Ds_icosaedro.csv";
        break;
    default:
        cerr << " p ammissibile ma q no" << endl;
        return false;
        break;
    }

flag = 1;
} 

if (q==3 && p!=3){ //duale di un solido a facce triangolari

    switch (p)
    {
    case 4: //ottaedro -> tetraedro duale
        mesh.nomefile0 = "../PolygonalMesh/Cell0Ds_tetraedro.csv";
        mesh.nomefile1 = "../PolygonalMesh/Cell1Ds_tetraedro.csv";
        mesh.nomefile2 = "../PolygonalMesh/Cell2Ds_tetraedro.csv";
        break;
    case 5: //dodecaedro -> icosaedro duale
        mesh.nomefile0 = "../PolygonalMesh/Cell0Ds_icosaedro.csv";
        mesh.nomefile1 = "../PolygonalMesh/Cell1Ds_icosaedro.csv";
        mesh.nomefile2 = "../PolygonalMesh/Cell2Ds_icosaedro.csv";
        break;
    default:
        cerr << "q ammissibile ma p no" << endl;
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

if ((b==0 && c >= 1) ||(c==0 && b >= 1)){
    mesh.classe = 1;
	mesh.d=max(b,c);

    switch (q) {
    case 3:
        V = 2*T + 2;
        E = 6*T;
        F = 4*T;
        break;
    case 4:
        V = 4*T + 2;
        E = 12*T;
        F = 8*T;
        break;
    case 5:
        V = 10*T + 2;
        E = 30*T;
        F = 20*T;
        break;
    }
    mesh.V = V;
    mesh.E = E;
    mesh.F = F;
    // mesh.n_spigoli = (mesh.d-1)*(6*mesh.q/(6-3*mesh.q+2*mesh.p)); //quanti spigoli verranno creati sulla superficie triangolata
    
    // unsigned int sum = mesh.n_spigoli*(mesh.d-1); // punti sui nuovi spigoli

    // for (unsigned int i = 2; i < mesh.d; i++) //Aggiunge anche punti interni 
    // {
    //     sum += i;
    // }
    
    // mesh.pt_totali = sum;
    flag = 1;
}

if (b==c && b >= 1){
    mesh.classe = 2;
    mesh.d=b;
    flag = 1;
}

if (flag == 0){
    cerr <<"b e c Non sono triangolazioni possibili" << endl;
    return false;
}


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


    // dobbimao fare la matrice di dimensione 3xnumnero di lunti totale
    mesh.M0D = Eigen::MatrixXd::Zero(3, mesh.V); 
// voglio stampare la lista

    for (string& li : lista_dim) // itero una lista
    {
        replace(li.begin(), li.end(), ';', ' ');
        istringstream convertitore(li);
        unsigned int id,marker;
        convertitore >> id >> marker >> mesh.M0D(0,id) >> mesh.M0D(1,id) >> mesh.M0D(2 ,id)   ;
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
          vector<unsigned int> linea1;
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
          vector<unsigned int> linea2;
          linea2.resize(n_vertici_spigoli);
          for (unsigned int j = 0; j < n_vertici_spigoli ; j++)
          {
            convertitore >> linea2[j];
          }

          mesh.M2D_spigoli.push_back(linea2);
          
        }     
        
    return true;
    }

bool Inizializzazione_vertici( PolygonalMesh& Pmesh )
{   
    unsigned int n = Pmesh.Dim0D; // numero di punti (serve per assegnare ID nuovi ai punti intermedi)
    unsigned int b = Pmesh.d;// numero di suddivisioni per lato, b+1 = numero di punti per spigolo
    //double l = Pmesh.lunghezza_lato_triangolino = (Pmesh.M0D.col(Pmesh.M1D(0,0))-Pmesh.M0D.col(Pmesh.M1D(1,0))).norm()/(b); // lunghezza per distanziare ogni punto nello spigolo, prendendo coord dei primi due punti dello spigolo zero, facendone la norma per trovare lunghezza A-B e la dividiamo per il numero di segmenti 
    //cambiamo la norma perche costosa 
    
    const auto s = VectorXd::LinSpaced(b+1,0.0,1.0);
    //matrice che conterrà l'ID del j-esimo punto sul i-esimo spigolo.   
    Pmesh.M_pt_spigoli= MatrixXi::Zero(Pmesh.Dim1D,b+1); 
   
    for (unsigned int i = 0; i < Pmesh.Dim1D; i++) // lo faccio per ogni spigolo
    {  
        Pmesh.M_pt_spigoli(i,0) = Pmesh.M1D(0,i); // ID punto iniziale
        Pmesh.M_pt_spigoli(i,b) = Pmesh.M1D(1,i); // ID punto finale

        Vector3d Acoord = Pmesh.M0D.col(Pmesh.M1D(0,i)); //estrae coordinate primo estremo da M0D
        Vector3d Bcoord = Pmesh.M0D.col(Pmesh.M1D(1,i)); //estrae coordinate secondo estremo da M0D

        //Creiamo i punti intermedi lungo lo spigolo
        for (unsigned int j = 1; j < b; j++) // ci muoviamo lungo interno saltando il primo e lultimo 
        {               
            unsigned int id = n+i*(b-1)+j-1; //id punti intermedi  :n è l'offset per partire dopo i vertici iniziali, i * (e - 2) tiene conto dei punti dei precedenti spigoli, j è il punto specifico su questo spigolo
            Pmesh.M_pt_spigoli(i,j) = id; //inserisce id punti intermedi nella matrice degli spigoli TOTALI
            //Vector3d versore = (Bcoord - Acoord).normalized(); //direzione spigolo
            //Vector3d coord = Acoord + versore * double(j) * l; 
            Vector3d coord = Acoord + s(j)* (Bcoord-Acoord); //coordinate punti intermedi
            Pmesh.M0D.col(id) = coord; //inserisce coordinate punti intermedi nella matrice globale
        }
    }	

    unsigned int id_spigoli= 0;
	Pmesh.M1D_triangolini = MatrixXi::Zero(2,Pmesh.E);
	Pmesh.M1D_spigoli_intermedi= MatrixXi::Zero(Pmesh.Dim1D, b);
	
    for (unsigned int k=0; k < Pmesh.Dim1D; k++)
    {
        for (unsigned int h=0; h < b; h++)
        {
            //inizialmente riempiamo M1D_triangolini inserendo per ogni id_spigoli gli id dei suoi vertici
            Pmesh.M1D_triangolini(0,id_spigoli) = Pmesh.M_pt_spigoli(k,h);
            Pmesh.M1D_triangolini(1,id_spigoli) = Pmesh.M_pt_spigoli(k,h+1);

            //ID dello spigolino h-esimo dello spigolo originale k-esimo.
            //Per ogni spigolo originale k, teniamo traccia dei suoi b spigolini

			Pmesh.M1D_spigoli_intermedi(k,h)=id_spigoli; 
            id_spigoli ++;
        }	 

    }	
    //Dopo aver eseguito tutta la funzione:
    //M0D conterrà i vertici originali + tutti i punti interpolati sugli spigoli
    //M_pt_spigoli ti dirà quali punti sono su ciascuno spigolo (per triangolare le facce)
}

Vector3d baricentro(const Vector3d& p0, const Vector3d& p1, const Vector3d& p2){
	Vector3d bar = (p0+p1+p2)/3.0;
	return bar;
}

bool Inizializzazione_punti_interni(PolygonalMesh& Pmesh){ //Iteriamo su tutte le facce Pmesh.M2D_vertici
    unsigned int d = Pmesh.d;
    unsigned int num_facce = Pmesh.Dim2D;
    unsigned int n_nuovi_punti = num_facce * (d - 2) * (d - 1) / 2;
    unsigned int id_pt_attuale = Pmesh.V - n_nuovi_punti; 
// E poi:
    unsigned int id_attuale_triangoli = 0;
    unsigned int id_attuale_spigoli = Pmesh.Dim1D*d; // numero di spigoli già creati (dalla triangolazione geodetica) + quelli che creeremo ora

    cout << "id_attuale_spigoli: " << id_attuale_spigoli << endl;

    Pmesh.M2D = MatrixXi::Zero(6,Pmesh.F); // 3 ID vertici, 3 ID spigoli per colonna che descrivono un triangolino ( F numero totale di triangolini )

    // ciclo sulle facce

    for (unsigned int faccia_id =0;  faccia_id < num_facce;  faccia_id++)
    {   
        // ID dei 3 vertici della faccia triangolare.
        unsigned int A_id = Pmesh.M2D_vertici[faccia_id][0];
        unsigned int B_id = Pmesh.M2D_vertici[faccia_id][1];
        unsigned int C_id = Pmesh.M2D_vertici[faccia_id][2];

        // coordinate 3D di quei vertici.
        //Vector3d A = Pmesh.M0D.col(A_id); // considero la colonna della matrice delle coordinate 3*ID_punti
        //Vector3d B = Pmesh.M0D.col(B_id); 
        //Vector3d C = Pmesh.M0D.col(C_id);

        //creo vettori per spostarsi lungo la triangolazione
		/*
        //ottengo vettore che si sposta in orizzontale
        Vector3d versore_orizzontale = (B-A) / (B-A).norm() ;
        Vector3d vettore_orizzontale = versore_orizzontale * l ;
        //ottengo vettore che si sposta in obliquo
        Vector3d versore_obliquo = (C-A) / (C-A).norm() ;
        Vector3d vettore_obliquo = versore_obliquo * l ;
		*/
		
       // //M2D_spigoli[faccia_id] è un vettore di 3 interi: gli ID degli spigoli. in questo modo trovo l'ID di ogni spigolo per faccia
        unsigned int AB_id = Pmesh.M2D_spigoli[faccia_id][0];
        unsigned int BC_id = Pmesh.M2D_spigoli[faccia_id][1];
        unsigned int CA_id = Pmesh.M2D_spigoli[faccia_id][2];


        //inserisco nel vettore tutti i punti sullo spigolo 
        VectorXi AB = Pmesh.M_pt_spigoli.row(AB_id);
        VectorXi BC = Pmesh.M_pt_spigoli.row(BC_id);
        VectorXi CA = Pmesh.M_pt_spigoli.row(CA_id);
		
        //inserisco nel vettore tutti gli id degli spigolini ( tra i punti nel vettore sopra ) dello spigolo 
		VectorXi AB_spigoli = Pmesh.M1D_spigoli_intermedi.row(AB_id);
		VectorXi BC_spigoli = Pmesh.M1D_spigoli_intermedi.row(BC_id);
		VectorXi CA_spigoli = Pmesh.M1D_spigoli_intermedi.row(CA_id);
		
        for (unsigned int t=0; t < d+1; t++)
        {
            cout << CA[t] << " ";
        }
        cout << endl;

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
            cout << "flagder: " << CA[0] << endl;
            VectorXi CA_reversed = CA.reverse().eval();
            CA = CA_reversed;
			VectorXi CA_spigoli_reversed = CA_spigoli.reverse().eval();
            CA_spigoli= CA_spigoli_reversed;
        }

        cout << "AB: " << AB[0] << " → " << AB[AB.size()-1] << " (expected A=" << A_id << ", B=" << B_id << ")" << endl;
        cout << "BC: " << BC[0] << " → " << BC[BC.size()-1] << " (expected B=" << B_id << ", C=" << C_id << ")" << endl;
        cout << "CA: " << CA[0] << " → " << CA[CA.size()-1] << " (expected C=" << C_id << ", A=" << A_id << ")" << endl;


        //Controlla che Pmesh.M2D_spigoli[faccia_id] abbia ID coerenti per ogni faccia 
        cout << "Faccia " << faccia_id << ": spigoli " << AB_id << ", " << BC_id << ", " << CA_id << endl;
       
        VectorXi base = AB; // vettore di tutti i punti tra A e B
		VectorXi base_spigoli = AB_spigoli; // vettore di id di spigolini tra A e B
        VectorXi tetto;
		VectorXi tetto_spigoli;

		unsigned int flag_spigoli = 0;
        for (unsigned int i = 0; i < d-1; i ++){

            Pmesh.M2D(0,id_attuale_triangoli)=base[0];//AB dovrà essere sovrascritto dalla nuova base
            Pmesh.M2D(1,id_attuale_triangoli)=base[1];
            Pmesh.M2D(2,id_attuale_triangoli)=CA[i+1];
            
            Pmesh.M2D(3,id_attuale_triangoli)=base_spigoli[0];
            
            Pmesh.M1D_triangolini(0,id_attuale_spigoli)=base[1];
            Pmesh.M1D_triangolini(1,id_attuale_spigoli)=CA[i+1];

            Pmesh.M2D(4,id_attuale_triangoli)=id_attuale_spigoli; //inserisce l'id spigolo (è) il secondo
			id_attuale_spigoli++;

			
            Pmesh.M2D(5,id_attuale_triangoli)=CA_spigoli[i];

			id_attuale_triangoli++;
            
            unsigned int dim = base.size();
            tetto.resize(dim-1);
            tetto[0]=CA[i+1];
			tetto[dim-2]=BC[i+1];
            tetto_spigoli.resize(dim-2);
			const auto s = VectorXd::LinSpaced(d-i,0.0,1.0);// linspaced da 0 a 1 con d-i punti

            // itero sui parallelogrammi dipo aver creato il primo triangolo 
            for (unsigned int j = 1; j < d - i; j ++) // con i + j == d , avremmo k=0 : Il punto cade su un lato (perché è una combinazione convessa di soli due vertici)
            {   
				
                if (j == d-i-1){ // uultimo parallelogramma 
                    //tetto[j]=BC[i+1];
					flag_spigoli = 1;
                }
                else{
                    //Vector3d pt = A+(i+1)*vettore_obliquo +j*vettore_orizzontale;
					Vector3d pt = Pmesh.M0D.col(tetto[0])+s(j)*(Pmesh.M0D.col(tetto[dim-2])-Pmesh.M0D.col(tetto[0]));
                    Pmesh.M0D.col(id_pt_attuale) = pt;
                    tetto[j]=id_pt_attuale;
                    id_pt_attuale++;
                }

                Pmesh.M2D(0,id_attuale_triangoli)=tetto[j-1];//AB dovrà essere sovrascritto dalla nuova base
                Pmesh.M2D(1,id_attuale_triangoli)=base[j];
                Pmesh.M2D(2,id_attuale_triangoli)=tetto[j];
				
				Pmesh.M2D(3,id_attuale_triangoli)=id_attuale_spigoli-1;


				Pmesh.M1D_triangolini(0,id_attuale_spigoli)=tetto[j-1];
				Pmesh.M1D_triangolini(1,id_attuale_spigoli)=tetto[j];
				tetto_spigoli[j-1] = id_attuale_spigoli;
                
                /*if (id_attuale_spigoli < 12){
                    cout << "scemo di guerra" << endl;
                    cout << id_attuale_spigoli << endl;
                }*/

				

				Pmesh.M2D(5,id_attuale_triangoli)=id_attuale_spigoli;
                id_attuale_spigoli++;
				
                Pmesh.M1D_triangolini(0,id_attuale_spigoli)=base[j]; //crea spoigolo obliquo interno
				Pmesh.M1D_triangolini(1,id_attuale_spigoli)=tetto[j];
				
				Pmesh.M2D(4,id_attuale_triangoli)=id_attuale_spigoli;
				id_attuale_spigoli++;

                id_attuale_triangoli++;

                Pmesh.M2D(0,id_attuale_triangoli)=base[j];//AB dovrà essere sovrascritto dalla nuova base
                Pmesh.M2D(1,id_attuale_triangoli)=base[j+1];
                Pmesh.M2D(2,id_attuale_triangoli)=tetto[j];
				
				Pmesh.M2D(5,id_attuale_triangoli)=id_attuale_spigoli-1;
				Pmesh.M2D(3,id_attuale_triangoli)=base_spigoli[j];
				
				if (flag_spigoli){
					Pmesh.M2D(4,id_attuale_triangoli)=BC_spigoli[i];
				}
				else{
				    Pmesh.M1D_triangolini(0,id_attuale_spigoli)=base[j+1]; //crea spoigolo obliquo interno
				    Pmesh.M1D_triangolini(1,id_attuale_spigoli)=tetto[j];
				    Pmesh.M2D(4,id_attuale_triangoli)=id_attuale_spigoli;
                    id_attuale_spigoli++;
				}
				id_attuale_triangoli++;
            }
            // lo riscrivo per poterlo riutilizzare nel ciclo successivo
			flag_spigoli= 0;
            base.resize(dim-1);
            base = tetto;
            base_spigoli.resize(dim-2);
            base_spigoli = tetto_spigoli;
        }
        Pmesh.M2D(0,id_attuale_triangoli)=CA[d-1];
        Pmesh.M2D(1,id_attuale_triangoli)=BC[d-1];
        Pmesh.M2D(2,id_attuale_triangoli)=CA[d];
			
        Pmesh.M2D(4,id_attuale_triangoli)=BC_spigoli[d-1]; //inserisce l'id spigolo 
		Pmesh.M2D(3,id_attuale_triangoli)=base_spigoli[0];
        Pmesh.M2D(5,id_attuale_triangoli)=CA_spigoli[d-1];

        id_attuale_triangoli++;

}
// stampa cella M0D

/*for (unsigned int i = 0; i < Pmesh.M_pt_spigoli.rows(); i++)
{
    for (unsigned int g=0; g < Pmesh.M_pt_spigoli.cols(); g++)
    {
        cout << Pmesh.M_pt_spigoli(i,g) << " ";
    }
    cout << endl;
}*/
cout << "M2D: " << endl;
for (unsigned int i = 0; i < Pmesh.M2D.cols(); i++)
{
    cout << Pmesh.M2D(0,i) << " " << Pmesh.M2D(1,i) << " " << Pmesh.M2D(2,i) << " "<< Pmesh.M2D(3,i) << " " << Pmesh.M2D(4,i) << " " << Pmesh.M2D(5,i) << endl;
}
/*
cout << "M0D: " << endl;
for (unsigned int i = 0; i < Pmesh.M0D.cols(); i++)
{
    cout << Pmesh.M0D(0,i) << " " << Pmesh.M0D(1,i) << " " << Pmesh.M0D(2,i) << endl;
}
cout << "M1D: " << endl;
for (unsigned int i = 0; i < Pmesh.M1D_triangolini.cols(); i++)
{
    cout << Pmesh.M1D_triangolini(0,i) << " " << Pmesh.M1D_triangolini(1,i) << endl;
}
*/
return true;
}

bool Proiezione_sfera(PolygonalMesh& Pmesh)
{
    //Proiezione sferica dei punti
     for (unsigned int i = 0; i < Pmesh.M0D.cols(); i++)
     {
        Vector3d coord = Pmesh.M0D.col(i);
        double r = coord.norm();
        Pmesh.M0D.col(i) = coord / r; // Normalizza le coordinate
    }
    return true;
}

bool stampa_geodetico(PolygonalMesh& Pmesh)
{
    ofstream file_C0D("Cell0Ds.txt"); 
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

    ofstream file_C1D("Cell1Ds.txt");
    if (!file_C1D.is_open())
    {
        cerr << "Error opening file Cell1Ds.txt" << endl;
        return false;
    }
    file_C1D << "ID;Vertice1;Vertice2" << endl;
    for (unsigned int i = 0; i < Pmesh.M1D_triangolini.cols(); i++)
    {
        file_C1D << i << ";" << Pmesh.M1D_triangolini(0,i) << ";" << Pmesh.M1D_triangolini(1,i) << endl;
    }
    file_C1D.close();

    ofstream file_C2D("Cell2Ds.txt");
    if (!file_C2D.is_open())
    {
        cerr << "Error opening file Cell2Ds.txt" << endl;
        return false;
    }
    file_C2D << "ID;Punti;Spigoli" << endl;
    for (unsigned int i = 0; i < Pmesh.M2D.cols(); i++)
    {
        file_C2D << i << ";";
        for (unsigned int j = 0; j < 3; j++)
        {
            file_C2D << Pmesh.M2D(j,i) << ";";
        }
        for (unsigned int j = 3; j < 6; j++)
        {
            file_C2D << Pmesh.M2D(j,i) << ";";
        }
        file_C2D.seekp(-1, std::ios::end); // vai all'ultimo carattere
        file_C2D.put(' '); // lo sovrascrivi (es. con uno spazio)
        file_C2D << endl; // vai a capo
    }
    file_C2D.close();

    ofstream file_C3D("Cell3Ds.txt");
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
	Dmesh.Dim0D = Pmesh.F;
	Dmesh.M1D_triangolini = -1* MatrixXi::Ones(2,Pmesh.E); //righe: flag, primo bar, secondo bar colonna id spigolo
	Dmesh.Dim1D = Pmesh.E;
	Dmesh.M2D = -1*MatrixXi::Ones(12,Pmesh.V);
	Dmesh.Dim2D = Pmesh.V;
	
	//MatrixXi spigoli_baricentri = -1* MatrixXd::ones(2,Pmesh.E); //righe: flag, primo bar, secondo bar colonna id spigolo
	unsigned int counter_if = 0;
		unsigned int counter_else = 0;
	for (unsigned int i=0; i < Pmesh.F; i++){
		Dmesh.M0D.col(i) = baricentro(Pmesh.M0D.col(Pmesh.M2D(0,i)), Pmesh.M0D.col(Pmesh.M2D(1,i)), Pmesh.M0D.col(Pmesh.M2D(2,i)));
		//cout << baricentro(Pmesh.M0D.col(Pmesh.M2D(0,i)), Pmesh.M0D.col(Pmesh.M2D(1,i)), Pmesh.M0D.col(Pmesh.M2D(2,i))) << endl;
		for (unsigned int j=3; j<6; j++){
			if (Dmesh.M1D_triangolini(0,Pmesh.M2D(j,i))== -1)
			{
				Dmesh.M1D_triangolini(0,Pmesh.M2D(j,i)) = i;
				counter_if++;
			}
			else{
				if (Dmesh.M1D_triangolini(1,Pmesh.M2D(j,i)) != -1){
					cout << "scemo di guerra sovrascrive 1" << endl;
				}
				Dmesh.M1D_triangolini(1,Pmesh.M2D(j,i)) = i;
				cout << Dmesh.M1D_triangolini(1,Pmesh.M2D(j,i)) << endl;
				counter_else++;
			}
		}
	}

	// cout << "M1D" << endl;
	// for (unsigned int i = 0; i < Pmesh.M1D_triangolini.cols(); i++)
	// {
    // cout << Dmesh.M1D_triangolini(0,i) << " " << Dmesh.M1D_triangolini(1,i) << endl;
	// }
	/*cout << "M0D: " << endl;
	for (unsigned int i = 0; i < Pmesh.M0D.cols(); i++)
	{
    cout << Dmesh.M0D(0,i) << " " << Dmesh.M0D(1,i) << " " << Dmesh.M0D(2,i) << endl;
	}
	*/
	
    MatrixXi id_spigoli = -1*MatrixXi::Ones(6,Pmesh.F);
	
	unsigned int id_vertice;
	unsigned int id_spigolo_precedente;
	unsigned int id_spigolo_successivo;
	unsigned int id_vertice_iniziale;
	
	for (unsigned int k =0; k < Pmesh.E; k++){

		unsigned int A_id = Pmesh.M1D_triangolini(0,k);// k che è l'id dello spigolo (e cpincide con la i di prima)  richiama gli id dei due vertici delo spigolo
		unsigned int B_id = Pmesh.M1D_triangolini(1,k);
		
		for (unsigned int j=0; j<6; j++){
			if (id_spigoli(j,A_id)== -1){
				id_spigoli(j,A_id)=k;
				break;
			}
		}
		for (unsigned int j=0; j<6; j++){
			if (id_spigoli(j,B_id)== -1){
				id_spigoli(j,B_id)=k;
				break;
			}
		}
	}
	
	for (unsigned int h = 0; h < Pmesh.V; h++){
		id_spigolo_precedente = id_spigoli(0,h);
		Dmesh.M2D(6,h) = id_spigolo_precedente;
		
		id_vertice_iniziale = Dmesh.M1D_triangolini(0,id_spigolo_precedente);
		Dmesh.M2D(0,h) = id_vertice_iniziale;
		
		//id_vertice = Dmesh.M1D_triangolini(1,id_spigolo_precedente);
		//Dmesh.M2D(1,h) = id_vertice;
		
		unsigned int contatore = 0;//vertice
		unsigned int contatore_spigoli = 6;
		
		while (id_vertice == id_vertice_iniziale){
			for(unsigned int j=1;j<6;j++){
				id_spigolo_successivo = id_spigoli(j,h);
				
				if (id_spigolo_successivo != id_spigolo_precedente){
					
					unsigned int A_id = Dmesh.M1D_triangolini(0,id_spigolo_successivo);
					unsigned int B_id = Dmesh.M1D_triangolini(1,id_spigolo_successivo);
					
					if (A_id == id_vertice){
						Dmesh.M2D(contatore+1, h) = A_id;
						Dmesh.M2D(contatore_spigoli+1,h) = id_spigolo_successivo;
						
						id_vertice = B_id;
						id_spigolo_precedente = id_spigolo_successivo;
						contatore++;
						contatore_spigoli++;
						break;
					}
					
					if (B_id == id_vertice){
						Dmesh.M2D(contatore+1, h) = B_id;
						Dmesh.M2D(contatore_spigoli+1,h) = id_spigolo_successivo;
						
						id_vertice = A_id;
						id_spigolo_precedente = id_spigolo_successivo;
						contatore++;
						contatore_spigoli++;
						break;
					}
				}
			}
		}
		
	}
	return true;
}
bool Inizializzazione_punti_interni_classe(PolygonalMesh& Pmesh){
	unsigned int d = Pmesh.d;
    unsigned int num_facce = Pmesh.Dim2D;
    unsigned int n_nuovi_punti = num_facce * pow(d,2);
    unsigned int id_pt_attuale = Pmesh.V - n_nuovi_punti;
    unsigned int id_attuale_triangoli = 0;
    unsigned int id_attuale_spigoli = Pmesh.Dim1D*2*d; // numero di spigoli già creati (dalla triangolazione geodetica) + quelli che creeremo ora

    //cout << "id_attuale_spigoli: " << id_attuale_spigoli << endl;
	Pmesh.M2D = MatrixXi::Zero(6,Pmesh.F);
	
	for (unsigned int faccia_id =0;  faccia_id < num_facce;  faccia_id++)
    {   
        // ID dei 3 vertici della faccia triangolare.
        unsigned int A_id = Pmesh.M2D_vertici[faccia_id][0];
        unsigned int B_id = Pmesh.M2D_vertici[faccia_id][1];
        unsigned int C_id = Pmesh.M2D_vertici[faccia_id][2];

        // coordinate 3D di quei vertici.
        Vector3d A = Pmesh.M0D.col(A_id); // considero la colonna della matrice delle coordinate 3*ID_punti
        Vector3d B = Pmesh.M0D.col(B_id); 
        Vector3d C = Pmesh.M0D.col(C_id);

        //creo vettori per spostarsi lungo la triangolazione
		/*
        //ottengo vettore che si sposta in orizzontale
        Vector3d versore_orizzontale = (B-A) / (B-A).norm() ;
        Vector3d vettore_orizzontale = versore_orizzontale * l ;
        //ottengo vettore che si sposta in obliquo
        Vector3d versore_obliquo = (C-A) / (C-A).norm() ;
        Vector3d vettore_obliquo = versore_obliquo * l ;
		*/
		
       // //M2D_spigoli[faccia_id] è un vettore di 3 interi: gli ID degli spigoli. in questo modo trovo l'ID di ogni spigolo per faccia
        unsigned int AB_id = Pmesh.M2D_spigoli[faccia_id][0];
        unsigned int BC_id = Pmesh.M2D_spigoli[faccia_id][1];
        unsigned int CA_id = Pmesh.M2D_spigoli[faccia_id][2];


        //inserisco nel vettore tutti i punti sullo spigolo 
        VectorXi AB = Pmesh.M_pt_spigoli.row(AB_id);
        VectorXi BC = Pmesh.M_pt_spigoli.row(BC_id);
        VectorXi CA = Pmesh.M_pt_spigoli.row(CA_id);
		
        //inserisco nel vettore tutti gli id degli spigolini ( tra i punti nel vettore sopra ) dello spigolo 
		VectorXi AB_spigoli = Pmesh.M1D_spigoli_intermedi.row(AB_id);
		VectorXi BC_spigoli = Pmesh.M1D_spigoli_intermedi.row(BC_id);
		VectorXi CA_spigoli = Pmesh.M1D_spigoli_intermedi.row(CA_id);

        for (unsigned int t=0; t < d+1; t++)
        {
            cout << CA[t] << " ";
        }
        cout << endl;

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
            cout << "flagder: " << CA[0] << endl;
            VectorXi CA_reversed = CA.reverse().eval();
            CA = CA_reversed;
			VectorXi CA_spigoli_reversed = CA_spigoli.reverse().eval();
            CA_spigoli= CA_spigoli_reversed;
        }

        cout << "AB: " << AB[0] << " → " << AB[AB.size()-1] << " (expected A=" << A_id << ", B=" << B_id << ")" << endl;
        cout << "BC: " << BC[0] << " → " << BC[BC.size()-1] << " (expected B=" << B_id << ", C=" << C_id << ")" << endl;
        cout << "CA: " << CA[0] << " → " << CA[CA.size()-1] << " (expected C=" << C_id << ", A=" << A_id << ")" << endl;


        //Controlla che Pmesh.M2D_spigoli[faccia_id] abbia ID coerenti per ogni faccia 
        cout << "Faccia " << faccia_id << ": spigoli " << AB_id << ", " << BC_id << ", " << CA_id << endl;
       
        VectorXi base = AB; // vettore di tutti i punti tra A e B
		VectorXi base_spigoli = AB_spigoli;
        VectorXi tetto;
		VectorXi tetto_spigoli;
		VectorXi Baricentri;
		unsigned int flag_spigoli = 0;
		VectorXd s;
		VectorXd ausiliario_pt;
		
		for (unsigned int i = 0; i < d; i ++){
			s = LinSpaced(d-i, 0.0,1.0);
			for (unsigned int l=0; l< d-i;l++){// iteriamo sulla lunghezza
			
				coordi = Pmesh.M0D.col(CA[(i+1)*2]);
				coordf = Pmesh.M0d.col(BC[i+1)*2);
				
				tetto = coordi +s(l)*(coordf-coordi);
			Pmesh.M0D.col(id_pt_attuale) = baricentro()
				
				
				baricentri
				}
			baricentri[]
	}
}

