/*
 *  auteur : RATEFIARISON Harivony Lalatiana
 */

#include <iostream> //inclure bibliothèque standard
#include <fstream>  //inclure le bibliothèque des flus de fichier
#include <string>   //inclusion du bibliothèque sur les chaine de charactère
#include <cmath>    //inclusion du bibliothèque pour les fonctions mathématiques

using namespace std;//esppace de nom standard

int taille (0);        // initialisation de la taille de la matrice et de vecteur

float **newMat(int rows, int cols){ //allocation dynamique de tableau rows*cols
    float **mat(NULL);
    mat = new float*[rows];

    for (int i=0;i<rows;i++)        
        mat[i] = new float[cols];
    return mat;
}

float *newVec(int n){               //alllouer un vecteur de taille n
    float* b;
    b  = new float [n];
    return b;
}

template <class T>                 //fonction pour echanger deux valeur a et b
void swap (T a, T b){              
    T tmp;
    tmp = a;
    a = b;
    a = tmp;
}

void displayMatrix(float** mat){    //afficher la matrice
    for (int i=0;i<taille;i++){
        for(int j=0; j<taille; j++)
            cout << mat[i][j] << " ";
        cout << endl;
    }
}

void displayVec (float* v){         //afficher un vecteur
    for(int i=0; i<taille; i++){
        cout << "[" << v[i] << "]" << endl;
    }
}

class Lsolver{
//METHODS
    public:
    //constructor
    Lsolver(){}
    Lsolver(string file);
    //destructor
    ~Lsolver();
    //display
    void displayMatTrian();
    void displayResult();
    //resolution
    void gaussElim();
    void solveMatrice();

    //accessors    
    float** getMatrix(){return matrice;}
    float* getB(){return b;}
    float* getSolution(){return solution;}
    int getDim(){return dim;}

//ATTRIBUTES
    private:
    float** matrice;            //matrice premier membre
    float*  b;                  //second membre du systeme
    float*  solution;           //solution
    int dim;
};

/*CONSTRUCTOR*/
Lsolver::Lsolver(string file){     

    ifstream fichier(file);
    if (fichier){
        fichier >> dim;
        taille = dim;
        matrice = newMat(dim,dim);  //allouer un nouveau tableau pour la matrice
        b = newVec(dim);            //alouer un nouveau tableau pour un vector
        
        for(int i=0; i<dim;i++){
            for(int j=0; j<dim;j++){
                fichier >> matrice [i][j];
            }
        }  
        
        for(int i=0; i<dim;i++){
                fichier >> b[i];
        }
    }
    fichier.close();
}

/*DESTRUCTOR*/

Lsolver::~Lsolver(){
    for(int y = 0; y < dim; y++)    //delete colonne
        delete[] matrice[y];
    delete[] matrice;               //supprimé ligne
    delete[] b;                     //supprimer le tableau contenant le second membre
    delete[] solution;              //supprimer le tableau contenant la solution
}

void Lsolver::displayMatTrian(){    //afficher le systeme triangularisé
    for (int i=0;i<dim;i++){        
        for(int j=0; j<dim; j++){
            if (matrice[i][j]==0)
                cout << "      ";
            else
                if(matrice[i][j]>0&&matrice[i][j-1]!=0)
                    cout <<  "+ " << matrice[i][j] << ".x" <<j+1 << " ";
                else
                    cout << matrice[i][j] << ".x" <<j+1 << " ";
        }
        cout << " = " << b[i] <<endl;
    }
}

void Lsolver::displayResult(){
    for(int i=0; i<dim; i++){
        cout << "x"<< i+1 << " = " << solution[i] << endl;
    }
}

void Lsolver::gaussElim(){
    solution = newVec(dim);
    
    //repeter l'operation jusqu'à dim-1
    for (int k=0; k<dim;k++){
        //trouver max
        float piv = abs(matrice[k][k]);
        int lpiv = k;
        
        for(int i=k;i<dim;i++){
            if(piv<abs(matrice[i][k])){
                piv = abs(matrice[i][k]);   //trouver le max |pivot|
                lpiv = i;                   //sauvegarder la igne
            }
        }
        //echanger la ligne k+1 au ligne lpiv
        for (int i=0; i<dim; i++){
            std::swap<float>(matrice[k][i],matrice[lpiv][i]);
            std::swap<float>(b[k],b[lpiv]);
        }
        /*
        cout << "debug :" << endl;
        displayMatrix(matrice);
        */
        //reduire par method de pivot
        for(int i=k+1;i<dim;i++){
            for (int j=k+1; j<dim; j++){  //or for (int j=k;j<dim;j++)
                if(matrice[i][j]==0)
                    continue;
                matrice[i][j] = matrice [i][j] - (matrice[i][k]/matrice[k][k])*matrice [k][j];
            }
            b[i] = b[i] - (matrice[i][k]/matrice[k][k])*b[k];
            matrice[i][k]=0;
        }
    }
}

void Lsolver::solveMatrice(){
    //résolution du systeme
    for(int i=dim-1;i>=0;i--){
        float sum=0;
        for(int j = i+1;j<dim;j++)
            sum += (matrice[i][j]*solution[j]);   
        solution[i] = (b[i]-sum)/(matrice[i][i]); 
    }
}

int main(){
    Lsolver s("matrice");   //initialiser le systeme depuis le fichier "matrice"

    /*le contenu du fichier "matrice" :
    3
    1 2 3
    2 4 5
    7 8 9
    4
    9
    29
    */

    cout << "Resolution d'un systeme linéaire A.x=b" << endl;
    cout << "La matrice A :" << endl;
    displayMatrix(s.getMatrix());
    cout << endl;
    
    cout << "Le second membre b" << endl;
    displayVec(s.getB());
    cout << endl;

    //Resolution
    cout << "Le probleme triangularisé :" << endl;
    s.gaussElim();
    s.solveMatrice();
    s.displayMatTrian();
    cout << endl;

    cout << "La solution est :" << endl;
    s.displayResult();
    
    return 0;   //fonction executé avec succes
}