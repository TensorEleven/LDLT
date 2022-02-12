#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

class Matrice
{
public:
    void rechercheProfil(ifstream &file);
    void recuperationMatrice();
    // void recuperationMatrice();

private:
    int dim, profil; 
    float** A;
    float* AP;
    int* p, *l;
};

int main()
{
    Matrice X; //déclaration d'une matrice X
    ifstream file;
    file.open("Matrice.txt"); //ouvrir le fichier
    X.rechercheProfil(file); //initialisation de la matrice 
    // X.recuperationMatrice();
}

void Matrice::rechercheProfil(ifstream &file)
{
    //récupération de la matrice dans le fichier
    file >> dim;

    //Allocation dynamique pour la matrice
    A = new float *[dim];

    for (int i = 0; i < dim; i++)
    {
        A[i] = new float[dim];
    }

    //Initialisation de la matrice
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            A[i][j] = 0;
        }
    }

    //Allocation dynamique du vecteur
    p = new int [dim];
    l = new int [dim];
    AP = new float [profil];

    for (int i = 0; i < dim; i++)
    {
        p[i] = 0;
        l[i] = 0;
        AP[i] = 0;
    }

    ///Récupération de la matrice dans le fichier
    int limite = 1;
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < limite; j++)
        {
            file >> A[i][j];
            cout<< A[i][j]<<" ";
        }
        cout<<endl;
        limite++;
    }

    //Recherche profil
    profil = 0;
    cout<<endl;
    cout<<"Le profil de la matrice A est:"<<endl<<"profil(A) = {"<<endl;
    for(int i=0; i<dim; i++)    //boucle pour la recherche du premier terme non nul de la ligne i 
    {
        for(int j=0; j<=i; j++)  //avec condition j<i pour avoir le demi profil
        {
            if(A[i][j] != 0)    //condition d'arrêt de la recherche
            {
                int ligne = 0;
                p[i] = j;       //récupération de l'indice de la première colonne non nulle
                // cout<<p[i]<<"   ";
                for (int k=j; k<=i; k++)
                {
                    cout<< " (" << i << "," << k <<")";     //affichage des éléments du profil
                    ligne = ligne+1;    //
                    if(k != dim-1)
                        cout<<",";
                    profil++;   //compteur élément du profil
                }
                l[i] = ligne;   //compteur élément du profil pour chaque ligne
                break; 
            }
        }
        cout<<endl;
    }
    cout<<" }"<<endl;
    cout<<"nombre element profil: "<< profil<< endl;
}

void Matrice::recuperationMatrice() 
{
    int k = 0;
    for(int i = 0; i<dim; i++)
    {
        int j = p[i];
        cout<<"pi: "<<p[i]<<"   li: "<<l[i]<<endl;
        for(int m = 0; m<l[i]; m++)
        {
            cout<<k<<": ";
            cout<<A[i][j]<<endl;
            AP[k] = A[i][j];
            k++;
            j++;
        }
        cout<<endl;
    }
}