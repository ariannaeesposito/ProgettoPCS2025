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
bool input_solido_platonico(PolygonalMesh& mesh, int argc ,char* argv[]){

if (argc != 7) 
    {
        cout << "Error: the number of input is not correct" << endl;
        return false;
    }

    unsigned int p, q, b, c, V, E, F, T , id_vertice_iniziale , id_vertice_finale ; 
    int flag = 0;

    istringstream convert0(argv[1]);
    convert0 >> p; //numero di lati per faccia
    istringstream convert1(argv[2]);
    convert1 >> q; //numero di facce che si incontrano in ogni vertice
    istringstream convert2(argv[3]);
    convert2 >> b; //parametro triangolazione geodetica
    istringstream convert3(argv[4]);
    convert3 >> c; //parametro triangolazione geodetica

    istringstream convert4(argv[5]);
    convert4 >> id_vertice_iniziale; // vertice iniziale per il cammino minimo

    istringstream convert5(argv[6]);
    convert5 >> id_vertice_finale; // vertice iniziale per il cammino minimo


    mesh.p = p; // riempiamo il contenitore mesh con i dati di input
    mesh.q = q;
    mesh.id_vertice_iniziale =  id_vertice_iniziale;
    mesh.id_vertice_finale =  id_vertice_finale;

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

    switch (q) {
    case 3:
        V = 4+6*(2*b-1)+4*((3*pow(b,2)-3*b)/2+1);
        E = 12*b+4*(9*pow(b,2)+3*b)/2;
        F = 12*b*(b+1);
        break;
    case 4:
        V = 6+12*(2*b-1)+8*((3*pow(b,2)-3*b)/2+1);
        E = 24*b+8*(9*pow(b,2)+3*b)/2;
        F = 24*b*(b+1);
        break;
    case 5:
        V = 12+30*(2*b-1)+20*((3*pow(b,2)-3*b)/2+1);
        E = 60*b+20*(9*pow(b,2)+3*b)/2;
        F = 60*b*(b+1);
        break;
    }
    mesh.V = V;
    mesh.E = E;
    mesh.F = F;

    //cout << mesh.V << " " << mesh.E << " " << mesh.F << endl;

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
    mesh.M0D = MatrixXd::Zero(3, mesh.V); 
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
    unsigned int b;

    if (Pmesh.classe == 1)
        b = Pmesh.d;// numero di suddivisioni per lato, b+1 = numero di punti per spigolo
    else 
        b = Pmesh.d*2;// numero di suddivisioni per lato, b+1 = numero di punti per spigolo

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
    return true;
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
        int A_id = Pmesh.M2D_vertici[faccia_id][0];
        int B_id = Pmesh.M2D_vertici[faccia_id][1];
        int C_id = Pmesh.M2D_vertici[faccia_id][2];

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
        int AB_id = Pmesh.M2D_spigoli[faccia_id][0];
        int BC_id = Pmesh.M2D_spigoli[faccia_id][1];
        int CA_id = Pmesh.M2D_spigoli[faccia_id][2];


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

bool crea_triangolo(PolygonalMesh& Pmesh,const unsigned int& id_triangolo ,const unsigned int& id_pt_1, const unsigned int& id_pt_2, const unsigned int& id_pt_3, const unsigned int& id_sp_1, const unsigned int& id_sp_2, const unsigned int& id_sp_3){
    Pmesh.M2D(0, id_triangolo) = id_pt_1;
    Pmesh.M2D(1, id_triangolo) = id_pt_2;
    Pmesh.M2D(2, id_triangolo) = id_pt_3;
    Pmesh.M2D(3, id_triangolo) = id_sp_1;
    Pmesh.M2D(4, id_triangolo) = id_sp_2;
    Pmesh.M2D(5, id_triangolo) = id_sp_3;
    return true;
}


bool Inizializzazione_punti_interni_classe2(PolygonalMesh& Pmesh){
    
	unsigned int d = Pmesh.d;
    unsigned int num_facce = Pmesh.Dim2D;
    unsigned int n_nuovi_punti = num_facce *  (pow(d,2)+((d-2)*(d-1))/2);
    unsigned int id_pt_attuale = Pmesh.V - n_nuovi_punti;
    unsigned int id_attuale_triangoli = 0;
    unsigned int id_attuale_spigolo = Pmesh.Dim1D*2*d; // numero di spigoli già creati (dalla triangolazione geodetica) + quelli che creeremo ora



    //cout << "id_attuale_spigoli: " << id_attuale_spigoli << endl;
	Pmesh.M2D = MatrixXi::Zero(6,Pmesh.F);
	
	for (unsigned int faccia_id =0;  faccia_id < num_facce;  faccia_id++)
    {   
        // ID dei 3 vertici della faccia triangolare.
        int A_id = Pmesh.M2D_vertici[faccia_id][0];
        int B_id = Pmesh.M2D_vertici[faccia_id][1];
        int C_id = Pmesh.M2D_vertici[faccia_id][2];

        // coordinate 3D di quei vertici.
        // Vector3d A = Pmesh.M0D.col(A_id); // considero la colonna della matrice delle coordinate 3*ID_punti
        // Vector3d B = Pmesh.M0D.col(B_id); 
        // Vector3d C = Pmesh.M0D.col(C_id);


		
       // //M2D_spigoli[faccia_id] è un vettore di 3 interi: gli ID degli spigoli. in questo modo trovo l'ID di ogni spigolo per faccia
        int AB_id = Pmesh.M2D_spigoli[faccia_id][0];
        int BC_id = Pmesh.M2D_spigoli[faccia_id][1];
        int CA_id = Pmesh.M2D_spigoli[faccia_id][2];

        //inserisco nel vettore tutti i punti sullo spigolo 
        VectorXi AB = Pmesh.M_pt_spigoli.row(AB_id);
        VectorXi BC = Pmesh.M_pt_spigoli.row(BC_id);
        VectorXi CA = Pmesh.M_pt_spigoli.row(CA_id).eval();
		
        //inserisco nel vettore tutti gli id degli spigolini ( tra i punti nel vettore sopra ) dello spigolo 
		VectorXi AB_spigoli = Pmesh.M1D_spigoli_intermedi.row(AB_id);
		VectorXi BC_spigoli = Pmesh.M1D_spigoli_intermedi.row(BC_id);
		VectorXi CA_spigoli = Pmesh.M1D_spigoli_intermedi.row(CA_id).eval();

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
		VectorXd s;

	
  

        //itero sui livelli del triangolo 
		for (unsigned int i = 0; i < d; i ++){

            tetto = VectorXi::Zero(d*2-2*i-1);
            tetto_spigoli = VectorXi::Zero(tetto.size()-1);
            //stampa tetto 

            tetto[0]= CA[2*i+2];
            tetto[tetto.size()-1]=BC[2*i+2];
            
            //linspace sul tetto di i 
			s = VectorXd::LinSpaced(d-i, 0.0, 1.0);

            Vector3d coord_iniziale = Pmesh.M0D.col(CA[(i+1)*2]);

            Vector3d coord_finale = Pmesh.M0D.col(BC[(i+1)*2]);

			for (unsigned int l=1; l < d-i-1; l++)
            {
                // iteriamo sulla lunghezza se = 2 non lo fa
                Pmesh.M0D.col(id_pt_attuale)= coord_iniziale +s(l)*(coord_finale-coord_iniziale);
				tetto[2*l] = id_pt_attuale; // tetto continene id dei punti sul tetto di i 
                id_pt_attuale++;
			}

            //PRIMO TRIANGOLINO 
            //inserisco baricentro nella matrice dei punti e aumento l'id da dare ai punti  
            Pmesh.M0D.col(id_pt_attuale)=baricentro(Pmesh.M0D.col(base[0]), Pmesh.M0D.col(base[2]), Pmesh.M0D.col(tetto[0]));
            id_pt_attuale++;

            //inserisco nuovi vertici che partono dal baricentro 
            Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(tetto[0], id_pt_attuale-1);
            id_attuale_spigolo++;
            
            Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(base[0], id_pt_attuale-1);
            id_attuale_spigolo++;
            
            Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(base[2], id_pt_attuale-1);
            id_attuale_spigolo++;

            Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(CA[2*i+1], id_pt_attuale-1);
            id_attuale_spigolo++;
            
            Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(base[1], id_pt_attuale-1);
            id_attuale_spigolo++;

            for (unsigned int k=0; k < d-i-1; k++) //itero sui parallelogrammi verso destra
             {     
               
                cout << " k= " << k << endl;

                //primo baricentro del parallelograma
                Pmesh.M0D.col(id_pt_attuale)=baricentro(Pmesh.M0D.col(base[2*k+2]), Pmesh.M0D.col(tetto[2*k+2]), Pmesh.M0D.col(tetto[k*2]));
                tetto[2*k+1]=id_pt_attuale;
                id_pt_attuale++;
                cout<<id_pt_attuale<<endl;
                //secondo baricentro del parallelograma
                Pmesh.M0D.col(id_pt_attuale)=baricentro(Pmesh.M0D.col(tetto[(k+1)*2]), Pmesh.M0D.col(base[2*k+2]), Pmesh.M0D.col(base[k*2+4]));
                id_pt_attuale++;
                
                //INSERIAMO GLI SPIGOLI NELLA MATRICE M1D_TRIANGOLINI

                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(id_pt_attuale-3, id_pt_attuale-2);
                id_attuale_spigolo++;

                //spigoli del tetto
                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(tetto[2*k], id_pt_attuale-2);
                tetto_spigoli[2*k] = id_attuale_spigolo;
                id_attuale_spigolo++;
                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(tetto[2*k+2], id_pt_attuale-2);
                tetto_spigoli[2*k+1] = id_attuale_spigolo;
                id_attuale_spigolo++;


            
                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(base[2*k+2], id_pt_attuale-2);
                id_attuale_spigolo++;

                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(id_pt_attuale-1, id_pt_attuale-2);
                id_attuale_spigolo++;

                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(base[2*k+2], id_pt_attuale-1);
                id_attuale_spigolo++;

                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(base[2*k+4], id_pt_attuale-1);
                id_attuale_spigolo++;

                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(tetto[2*k+2], id_pt_attuale-1);
                id_attuale_spigolo++;

                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(base[2*k+3], id_pt_attuale-1);
                id_attuale_spigolo++;

                if(k==d-i-2){

                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(BC[2*i+1], id_pt_attuale-1);
                id_attuale_spigolo++;

                }

                /*crea_triangolo(Pmesh, id_attuale_triangoli, base[2*k], base[2*k+1], id_pt_attuale-2, base_spigoli[2*k], id_attuale_spigolo-6, id_attuale_spigolo-7);
                id_attuale_triangoli++;
                crea_triangolo(Pmesh, id_attuale_triangoli, base[2*k+1], base[2*k+2], id_pt_attuale-2, base_spigoli[2*k+1], id_attuale_spigolo-5, id_attuale_spigolo-6);
                id_attuale_triangoli++;
                crea_triangolo(Pmesh, id_attuale_triangoli, base[2*k], id_pt_attuale-2,CA[i*2+1], id_attuale_spigolo-7, id_attuale_spigolo-4, CA_spigoli[i*2]);
                id_attuale_triangoli++;
                crea_triangolo(Pmesh, id_attuale_triangoli,CA[i*2+1], id_pt_attuale-2, tetto[2*k], id_attuale_spigolo-4, id_attuale_spigolo-3, CA_spigoli[i*2+1]);
                id_attuale_triangoli++;
                crea_triangolo(Pmesh, id_attuale_triangoli, id_pt_attuale-2, id_pt_attuale-1, tetto[2*k], id_attuale_spigolo-2, id_attuale_spigolo-9, id_attuale_spigolo-3);
                id_attuale_triangoli++;
                crea_triangolo(Pmesh, id_attuale_triangoli, base[2*k+2], id_pt_attuale-1, id_pt_attuale-2, id_attuale_spigolo-1, id_attuale_spigolo-2, id_attuale_spigolo-5);
                id_attuale_triangoli++;*/

                //CA[i*2+1] = id_pt_attuale-1;
                

                //CA_spigoli[2*i] = id_attuale_spigolo-1;
               // CA_spigoli[2*i+1] = id_attuale_spigolo-8;

            }
            
            if (i==d-1){

                Pmesh.M1D_triangolini.col(id_attuale_spigolo) = Vector2i(BC[2*i+1], id_pt_attuale-1);
                id_attuale_spigolo++;

                }

            /*crea_triangolo(Pmesh, id_attuale_triangoli, id_base_1, id_base_mezzo, id_pt_attuale-1, base_spigoli[dim-i*2-3], id_attuale_spigolo-5, id_attuale_spigolo-6);
                id_attuale_triangoli++;
            crea_triangolo(Pmesh, id_attuale_triangoli, id_base_mezzo, id_base_2, id_pt_attuale-1, base_spigoli[dim-i*2-2], id_attuale_spigolo-4, id_attuale_spigolo-5);
                id_attuale_triangoli++;
            crea_triangolo(Pmesh, id_attuale_triangoli, id_pt_attuale-1, id_base_2, BC[i*2+1], id_attuale_spigolo-4, BC_spigoli[i*2], id_attuale_spigolo-3);
                id_attuale_triangoli++;
            crea_triangolo(Pmesh, id_attuale_triangoli,id_pt_attuale-1, BC[i*2+1], id_tetto, id_attuale_spigolo-3, BC_spigoli[i*2+1], id_attuale_spigolo-2);
                id_attuale_triangoli++;
            crea_triangolo(Pmesh, id_attuale_triangoli, CA[i*2+1], id_pt_attuale-1, id_tetto, id_attuale_spigolo-1, id_attuale_spigolo-2, CA_spigoli[i*2+1]);
                id_attuale_triangoli++;
            crea_triangolo(Pmesh, id_attuale_triangoli, CA[i*2+1], id_pt_attuale-1, id_base_1, id_attuale_spigolo-1, id_attuale_spigolo-6, CA_spigoli[i*2]);
                id_attuale_triangoli++;*/

            //base.resize(dim-(i+1)*2-2);
            //base_spigoli.resize(dim-(i+1)*2-3);

            base.resize(tetto.size());
    
            base = tetto;

            base_spigoli.resize(tetto_spigoli.size());
            base_spigoli = tetto_spigoli;
            // dim = base.size();

            cout << "fine riga: " << i << endl;

        }
} 
return true;
}

bool CamminoMinimo(PolygonalMesh& Pmesh){
    
   unsigned int  id_vertice_iniziale = Pmesh.id_vertice_iniziale;
   unsigned int  id_vertice_finale = Pmesh.id_vertice_finale;
    
    if (id_vertice_iniziale >= Pmesh.V || id_vertice_iniziale < 0 || id_vertice_finale >= Pmesh.V || id_vertice_finale < 0 ){
       cerr << "gli id non sono validi" << endl;
       return false;
    }

    //creiamo la LISTA DI ADIACENZA
    vector<vector<unsigned int>> LA;
    LA.reserve(Pmesh.V);
    MatrixXd W = MatrixXd::Zero(Pmesh.V, Pmesh.V);
    
    //riempo la matrice di adiacenza

   /*for (size_t j = 0 ; j < Pmesh.V ; j ++ ){

        vector<unsigned int> punti_vicini;

       //considero ogni spigolo e ne prendo i punti 
        for (size_t i = 0 ; i < Pmesh.E ; i ++){
        cout<<"uff"<< j <<endl;

            unsigned int origine = Pmesh.M1D_triangolini(0,i);
            unsigned int fine = Pmesh.M1D_triangolini(1,i);

        // se origine del segmento che considero ha id = i , inserisco la fine del segmento nei punti vicini
        // se fine del segmento che considero ha id = i , inserisco l'inizio del segmento nei punti vicini

            if (origine == j)
                punti_vicini.push_back(fine);

            else if (fine == j)
                punti_vicini.push_back(origine);
        }
    
        //per ogni punto inserisco nel vector di vector i suoi vicini 
        LA.push_back(punti_vicini);
    }*/ 
    

for (size_t i = 0; i < Pmesh.E; ++i) {
    unsigned int origine = Pmesh.M1D_triangolini(0, i);
    unsigned int fine    = Pmesh.M1D_triangolini(1, i);

    if (origine < Pmesh.V && fine < Pmesh.V) {
        LA[origine].push_back(fine);
        LA[fine].push_back(origine);  // grafo non orientato
    } else {
        cerr << "⚠️ Errore: origine o fine fuori dai limiti per lo spigolo " << i << endl;
    }
}

    //itero sui punti

    for ( size_t i=0; i< LA.size(); i++){
       
        //itero sul vettore di punti vicini
        for (size_t k = 0; k < LA[i].size(); ++k) {

           

        unsigned int w = LA[i][k];  // vicino del punto i 
        if (i < Pmesh.V && w < Pmesh.V) {
        W(i, w) = (Pmesh.M0D.col(w) - Pmesh.M0D.col(i)).norm();
    } else {
        cerr << " Indice fuori dai limiti: i = " << i << ", w = " << w << std::endl;
    }
    }
    }

    for (size_t i = 0; i < LA.size(); ++i) {
    std::cout << "Nodo " << i << " ha vicini: ";
    for (unsigned int w : LA[i]) {
        std::cout << w << " ";
    }
    std::cout << std::endl;
}

    
    //ALGORITMO DI DIJKSTRA per visitare il grafo

    //inizilizzo il vettore predecessore a -1
    vector<unsigned int> pred(Pmesh.V,-1); // inizilizzo il predecessore a -1
 
    //inizilizzo il vettore delle distanze a un numero molto grande
	vector<double> dist(Pmesh.V, numeric_limits<double>::infinity());

    //Ogni elemento della coda è una coppia formata da: int (ID di un nodo) e double (distanza)

    using ND = pair<double, int>; //coppia nodo-distanza !!!!!!!CONTROLLARE ORDINE

	// priority_queue<  Tipo,   Contenitore,     Criterio >: classe STL per le code con priorità (distanza)
    // greater(criterio) per avere che il minimo sia estratto per primo ( sia messo in cima ed estratto con top )
    // vector(contenitore) per supportare accesso casuale

    priority_queue< ND, vector<ND>, greater<ND>> PQ;
		
    pred[id_vertice_iniziale] = id_vertice_iniziale;
    dist[id_vertice_iniziale] = 0;

    for(int i = 0; i < LA.size(); i++){
        PQ.push({dist[i], i}); 
    }
			
    while(!PQ.empty()){

        int p = PQ.top().first; // punto
        int u= PQ.top().second; // distanza

        PQ.pop(); //rimuove ( non restituisce nulla ) l’elemento con priorità più alta 

        // Se la distanza più aggiornata per u (dist[u]) è già più piccola di quella che sto estraendo (p) ,  non la uso.
        if (dist[u]<p)
            continue;
        
        for(const auto& w : LA[u]){ //foreach
            if( dist[w] > dist[u] + W(u,w) ) {
                //Aggiorna distanza se hai trovato un cammino più corto
                dist[w] = dist[u] + W(u,w); 
                pred[w] = u;

                PQ.push({w, dist[w]});
            }
        }
	}

    // path contiene id dei vertici del cammino minimo 
    vector<unsigned int> path;
    path.reserve(Pmesh.V);

    int v = id_vertice_finale;

    while(v != id_vertice_iniziale) {
        path.push_back(v);
        v = pred[v];
	} 
	
    path.push_back(id_vertice_iniziale);
		
    // coloriamo i punti relativi al percorso
    vector<double> ProprietaPunti_path(Pmesh.V, 0.0);

	for(unsigned int punto : path){
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
	utilities.ExportPoints("./Cell0Ds.inp",
								Pmesh.M0D,
							ProprietaPunti);


    // coloriamo i lati relativi al percorso e calcoliamo la lunghezza del percorso
    vector<unsigned int> segmenti_Percorso;
    //segmenti_Percorso.reserve(path.size()-1);

    vector<double> ProprietaLatiPercorso(Pmesh.V, 0.0 );

    cout << path.size() << endl;

    double lunghezzaPercorso = 0.0;

    //itero su ogni coppia di punti adiacenti del percorso 
    for (unsigned int i = 0; i < path.size() - 1; ++i) {

        unsigned int pt_a = path[i];
        unsigned int pt_b = path[i + 1];

        lunghezzaPercorso += (Pmesh.M0D.col(pt_a) - Pmesh.M0D.col(pt_b)).norm();

        // cerco lo spigolo che collega pt_a e pt_b
        //itero su ogni id_spigolo

        for (unsigned int j = 0; j<  Pmesh.M1D_triangolini.cols(); j ++ )
        {
            //id dei punti estremi del segmento
            unsigned int id_A = Pmesh.M1D_triangolini(0, j);
            unsigned int id_B = Pmesh.M1D_triangolini(1, j);

            if ((id_A == pt_a && id_B== pt_b) || (id_B == pt_a && id_A == pt_b)) {

                segmenti_Percorso.push_back(j); //inserisco id dello spigolo
                ProprietaLatiPercorso[j] = 1.0;
                break; // trovato lo spigolo, esco dal ciclo interno
            }
        }
    } 

    ShortPathProperty.Label = "Percorso minimo";
    ShortPathProperty.Size = ProprietaLatiPercorso.size();
    ShortPathProperty.NumComponents = 1;
    ShortPathProperty.Data = ProprietaLatiPercorso.data();  

    vector<Gedim::UCDProperty<double>> Segmenti_Proprietà;

    Segmenti_Proprietà.push_back(ShortPathProperty);

    utilities.ExportSegments("./Cell1Ds.inp",
								  Pmesh.M0D,
								  Pmesh.M1D_triangolini,
                    {},
                    Segmenti_Proprietà);

 return true;

}




}

