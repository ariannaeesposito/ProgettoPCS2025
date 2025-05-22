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
        mesh.nomefile0 = "../PolygonalMesh/Cell0Ds_tetraedro.csv";
        mesh.nomefile1 = "../PolygonalMesh/Cell1Ds_tetraedro.csv"; 
        mesh.nomefile2 = "../PolygonalMesh/Cell2Ds_tetraedro.csv";
        break;
    case 4: //ottaedro 
        mesh.nomefile0 = "../PolygonalMesh/Cell0Ds_ottaedro.csv";
        mesh.nomefile1 = "../PolygonalMesh/Cell1Ds_ottaedro.csv";
        mesh.nomefile2 = "../PolygonalMesh/Cell2Ds_ottaedro.csv";
        break;
    case 5: //icosaedro
        mesh.nomefile0 = "../PolygonalMesh/Cell0Ds_icosaedro.csv";
        mesh.nomefile1 = "../PolygonalMesh/Cell1Ds_icosaedro.csv";
        mesh.nomefile2 = "../PolygonalMesh/Cell2Ds_icosaedro.csv";
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

//struct per comparare ( nella mappa srve qualcosa che abbia implementato in se loperatore di disuguaglianza)

/*struct Vector3dComparator {
    bool operator()(const Vector3d& a, const Vector3d& b) const {
        for (int i = 0; i < 3; ++i) {
            if (a[i] < b[i]) return true;
            if (a[i] > b[i]) return false;
        }
        return false;
    }
};*/

// unsigned int accedimappa(PoligonalLibrary::TriangularMesh& mesh, PolygonalMesh& Pmesh, Vector3d& coord){

// 	double tol = pow(10,-15); //definiamo tolleranza

//     //controlliamo se gia presente il puntio 

// 	for (const auto& coppia : mesh.coordinate_punti){ //iteriamo sulla mappa
// 		double diff = (coppia.first - coord).norm();
//     //calcoliamo la differnaza tra i due set di coordinate
// 		if (diff < tol){ //se sono vicine allora le consideriamo le stesse
// 			return coppia.second;  //restituiamo l'id associato a quelle coordinate ( gia presente )
// 		}
// 	}
// 	//se non c'è nella mappa lo creiamo nuovo 

// 	unsigned int dim = Pmesh.Dim0D; 
//     mesh.coordinate_punti[coord] = dim;       // assegna nuovo ID al punto. Inserisce coord come nuova chiave nella mappa coordinate_punti, con valore dim (cioè l’ID).
//     Pmesh.Dim0D += 1;

//     Pmesh.M0D(0, dim) = coord(0);    //Salva le coordinate del nuovo punto nella matrice M0D
//     Pmesh.M0D(1, dim) = coord(1);
//     Pmesh.M0D(2, dim) = coord(2);

//     return dim;
    
    
  /*unsigned int dim = Pmesh.Dim0D; //dim corrente

	mesh.coordinate_punti[coord] = dim+1; //cell0ds già piena stiamo, stiamo aggiungendo alla mappa una chiave e l'id relativo
	Pmesh.Dim0D += 1;  //aumentiamo il conteggio dei punti già creati
	Pmesh.M0D(0, dim+1)=coord(0);  //mettiamo il nuovo punto dentro la polygonal mesh
	Pmesh.M0D(1, dim+1)=coord(1);
	Pmesh.M0D(2, dim+1)=coord(2);
	return dim+1;*/



bool Inizializzazione_vertici( PolygonalMesh& Pmesh , TriangularMesh& Tmesh)
{   
    unsigned int n = Pmesh.Dim0D; // numero di punti (serve per assegnare ID nuovi ai punti intermedi)
    
    //d è il numero di suddivisioni per lato
    unsigned int e = Pmesh.d+1;// numero di punti per ogni spigolo

    Tmesh.M_pt_spigoli= MatrixXd::Zero(Pmesh.Dim1D,e); //matrice che conterrà l'ID del j-esimo punto sul i-esimo spigolo.
    // stampa M_pt_spigoli
   
    Pmesh.lunghezza_lato_triangolino = (Pmesh.M0D.col(Pmesh.M1D(0,0))-Pmesh.M0D.col(Pmesh.M1D(1,0))).norm()/(e-1); // lunghezza per distanziare ogni punto nello spigolo, prendendo coord dei primi due punti dello spigolo zero, facendone la norma per trovare lunghezza A-B e la dividiamo per il numero di segmenti 
    
    //assumiamo che tutti gli spigoli abbiano lunghezza uguale (giusto per solidi platonici)

    for (unsigned int i = 0; i < Pmesh.Dim1D; i++) // lo faccio per ogni spigolo

    {  
        Tmesh.M_pt_spigoli(i,0) = Pmesh.M1D(0,i); // ID punto iniziale
        Tmesh.M_pt_spigoli(i,e-1) = Pmesh.M1D(1,i); // ID punto finale

        //Prendiamo le coordinate degli estremi
        Vector3d Acoord=Pmesh.M0D.col(Pmesh.M1D(0,i)); //coordinate primo estremo 
        Vector3d Bcoord=Pmesh.M0D.col(Pmesh.M1D(1,i)); //coordinate secondo estremo 

        //Creiamo i punti intermedi lungo lo spigolo
        for (unsigned int j = 1; j < e-1; j++) // ci muoviamo lungo interno saltando il primo e lultimo 
        {               

            unsigned int id = n+i*(e-2)+j-1; //id punti intermedi  :n è l'offset per partire dopo i vertici iniziali, i * (e - 2) tiene conto dei punti dei precedenti spigoli, j è il punto specifico su questo spigolo
            Tmesh.M_pt_spigoli(i,j) = id; //inserisce id punti intermedi nella matrice degli spigoli TOTALI
            Vector3d versore = (Bcoord - Acoord).normalized(); //direzione spigolo

            Vector3d coord = Acoord + versore * double(j) * Pmesh.lunghezza_lato_triangolino; //coordinate punti intermedi

            Pmesh.M0D.col(id) = coord; //inserisce coordinate punti intermedi nella matrice globale
        }
        //creiamo gli spigolimedi
        
    }
    Pmesh.M1D_triangolini = MatrixXd::Zero(2,Pmesh.E);
        for (unsigned int k=0; k < Pmesh.Dim1; k++)
        {
            for (unsigned int h=0; h < e-1; h++)
            {
                Pmesh.M1D_triangolini(0,id) = Tmesh.M_pt_spigoli(k,h);
                Pmesh.M1D_triangolini(1,id) = Tmesh.M_pt_spigoli(k,h+1);
            }
        }
    //Dopo aver eseguito tutta la funzione:
    //M0D conterrà i vertici originali + tutti i punti interpolati sugli spigoli
    //M_pt_spigoli ti dirà quali punti sono su ciascuno spigolo (per triangolare le facce)
}

/*bool Inizializzazione_punti_interni(PolygonalMesh& Pmesh, TriangularMesh& Tmesh){ //Iteriamo su tutte le facce Pmesh.M2D_vertici

    unsigned int d = Pmesh.d;
    unsigned int n_attuale = Pmesh.M0D.cols(); // n. punti attuali ( colonne della matrice delle coordinate)
    unsigned int id_punto = n_attuale;

    vector<Vector3d> nuovi_punti;

    //ciclo sulle facce

    for (const auto& faccia : Pmesh.M2D_vertici)
    {   
        // ID dei 3 vertici della faccia triangolare.
        unsigned int A_id = faccia[0];
        unsigned int B_id = faccia[1];
        unsigned int C_id = faccia[2];

        //coordinate 3D di quei vertici.
        Vector3d A = Pmesh.M0D.col(A_id);
        Vector3d B = Pmesh.M0D.col(B_id);
        Vector3d C = Pmesh.M0D.col(C_id);

        for (unsigned int i = 1; i < d; ++i)
        {
            for (unsigned int j = 1; j < d - i; ++j) // con i + j == d , avremmo k=0 : Il punto cade su un lato (perché è una combinazione convessa di soli due vertici)
            {
                unsigned int k = d - i - j;

                // parametri della combinazione convessa , sommano a 1 ; coordinate baricentriche
                double a = double(i) / d; 
                double b = double(j) / d;
                double c = double(k) / d;

                Vector3d P = a * A + b * B + c * C;
                nuovi_punti.push_back(P);
            }
        }
    }

    // Aggiungiamo questi punti a M0D
    unsigned int n_nuovi = nuovi_punti.size();
    Pmesh.M0D.conservativeResize(3, id_punto + n_nuovi);
    for (unsigned int i = 0; i < n_nuovi; ++i)
    {
        Pmesh.M0D.col(id_punto + i) = nuovi_punti[i];
    }

    cout << "Aggiunti " << n_nuovi << " punti interni" << endl;
    return true;
}*/

bool Inizializzazione_punti_interni(PolygonalMesh& Pmesh, TriangularMesh& Tmesh){ //Iteriamo su tutte le facce Pmesh.M2D_vertici
    double l = Pmesh.lunghezza_lato_triangolino;
    unsigned int d = Pmesh.d;
    unsigned int num_facce = Pmesh.Dim2D;
    unsigned int n_nuovi_punti = num_facce * (d - 2) * (d - 1) / 2;
    unsigned int id_pt_attuale = Pmesh.V - n_nuovi_punti; 
// E poi:
    unsigned int id_attuale_triangoli = 0;
    unsigned int id_attuale_spigoli = Pmesh.Dim1D;

    Pmesh.M2D = MatrixXd::Zero(6,Pmesh.F);

    //ciclo sulle facce
    cout << "prova" << endl;
    cout << Pmesh.Dim2D << endl;

    for (unsigned int faccia_id =0;  faccia_id < num_facce;  faccia_id++)
    {   
        // ID dei 3 vertici della faccia triangolare.
        unsigned int A_id = Pmesh.M2D_vertici[faccia_id][0];
        unsigned int B_id = Pmesh.M2D_vertici[faccia_id][1];
        unsigned int C_id = Pmesh.M2D_vertici[faccia_id][2];


        //coordinate 3D di quei vertici.
        Vector3d A = Pmesh.M0D.col(A_id); 
        Vector3d B = Pmesh.M0D.col(B_id);
        Vector3d C = Pmesh.M0D.col(C_id);
        Vector3d versore_orizzontale = (B-A) / (B-A).norm() ;
        Vector3d vettore_orizzontale = versore_orizzontale * l ;
        Vector3d versore_obliquo = (C-A) / (C-A).norm() ;
        Vector3d vettore_obliquo = versore_obliquo * l ;
       
        unsigned int AB_id = Pmesh.M2D_spigoli[faccia_id][0]; //M2D_spigoli[faccia_id] è un vettore di 3 interi: gli ID degli spigoli.
        unsigned int BC_id = Pmesh.M2D_spigoli[faccia_id][1];
        unsigned int CA_id = Pmesh.M2D_spigoli[faccia_id][2];
        VectorXd AB = Tmesh.M_pt_spigoli.row(AB_id);
        VectorXd BC = Tmesh.M_pt_spigoli.row(BC_id);
        VectorXd CA = Tmesh.M_pt_spigoli.row(CA_id);

        cout << "check 1" << endl;
        // orientiamo i vertici 
        if (AB[0]==B_id){
            AB = AB.reverse();
        }
        if (BC[0]==C_id){
            BC = BC.reverse();
        }
        if (CA[0]==C_id){
            CA = CA.reverse();
        }

        VectorXd base = AB; // vettore di tutti i punti tra A e B
        VectorXd tetto;


        for (unsigned int i = 0; i < d; i ++){

            Pmesh.M2D(0,id_attuale_triangoli)=base[0];//AB dovrà essere sovrascritto dalla nuova base
            Pmesh.M2D(1,id_attuale_triangoli)=base[1];
            Pmesh.M2D(2,id_attuale_triangoli)=CA[i+1];
            

            Pmesh.M1D_triangolini(0,id_attuale_spigoli)=base[1];
            Pmesh.M1D_triangolini(1,id_attuale_spigoli)=CA[i+1];

            Pmesh.M2D(2,id_attuale_triangoli)=id_attuale_spigoli; //inserisce l'id spigolo (è) il secondo

            id_attuale_triangoli++;
            id_attuale_spigoli++;

            unsigned int dim = base.size();
            tetto.resize(dim-1);
            tetto[0]=CA[i+1];
        
            for (unsigned int j = 1; j < d - i; j ++) // con i + j == d , avremmo k=0 : Il punto cade su un lato (perché è una combinazione convessa di soli due vertici)
            {   
                if (j == d-i-1){
                    tetto[j]=BC[i+1];
                }
                else{
                    Vector3d pt = A+(i+1)*vettore_obliquo +j*vettore_orizzontale;
                    Pmesh.M0D.col(id_pt_attuale) = pt;
                    tetto[j]=id_pt_attuale;
                    id_pt_attuale++;
                }
                
                Pmesh.M2D(0,id_attuale_triangoli)=tetto[j-1];//AB dovrà essere sovrascritto dalla nuova base
                Pmesh.M2D(1,id_attuale_triangoli)=base[j];
                Pmesh.M2D(2,id_attuale_triangoli)=tetto[j];
                id_attuale_triangoli++;
                Pmesh.M2D(0,id_attuale_triangoli)=base[j];//AB dovrà essere sovrascritto dalla nuova base
                Pmesh.M2D(1,id_attuale_triangoli)=base[j+1];
                Pmesh.M2D(2,id_attuale_triangoli)=tetto[j];
                id_attuale_triangoli++;>
            }
            base.conservativeResize(dim-1);
            base = tetto;
        }


}
// stampa cella M0D
cout << "M0D: " << endl;
for (unsigned int i = 0; i < Pmesh.M0D.cols(); i++)
{
    cout << Pmesh.M0D(0,i) << " " << Pmesh.M0D(1,i) << " " << Pmesh.M0D(2,i) << endl;
}

return true;

}



}
