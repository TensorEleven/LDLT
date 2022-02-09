/*
 *  RATEFIARISON 
 *  Harivony Lalatiana
 * 
 *  Factorisation LDLt
 * 
 */

#include <iostream>     // bibliothèque standard d'entrE et de sortie
#include <fstream>      // bibliothèque d'entrE et de sortie pour les fichier
#include <string>       // bibliotheque de manipulation de chaine de charactere

using namespace std;    // utiliser l'espace de nom standard std


int dim(10);            // definir la taille du systeme 

/*
    Allocation/Affichage de matrice et de vecteur
*/

float **newMat(int rows, int cols);     // creer une matrice rows x cols
float *newVec(int n);                   // creer un vecteur de taille n
void displayMat(float** mat);           // afficher une matrice de mat
void displayVec (float* v);             // afficher un vecteur 

/*
 *
 *  CLASS: LDLt
 * 
 * 
*/

class LDLT{
// Methods
public:
    LDLT();
    ~LDLT();
    void loadAb();
    void computeP_i();
    void factrizeA();

    // getter
    float** getA(){return A;}
    float* getb(){return b;}
    float** getD(){return D;}
    float* getX(){return x;}

    void solveTriangInf(float **mat, float *vec); 
    void transpose(float **mat);
	void solveTriangSup(float **mat, float *vec);
    void solveTriangSupNoTranspose(float **mat, float *vec);
    void profil(float** mat, string name);
    void getProfil(float* mat, string name);
    void solve();

// Attributs
private:
    float** A;
    float** L;
    float** D;
    float* x;    // la solution

    float* b;    // second membre
    float* p_i;  // indice des premieres element non null
};

/*
 *
 *  MAIN
 * 
 * 
*/

int main (){
    // initialiser une instance de LDLt 
    LDLT solver;

    // Interface Utilisateur
    cout << "____________________________________" << endl;
    cout << "\n  Factorisation LDLt et Resolution " << endl;
    cout << "____________________________________" << endl << endl;

    cout << "Factoriser une matrice symetric positive A\n tq: A = L.D.Lt"<< endl<< endl;
    // factorisation et resolution
    solver.solve();

    // succes de l'execution
    return 0;
}

/***** MATRICE - VECTEUR ******/

float **newMat(int rows, int cols){ //allocation dynamique de tableau rows*cols
    float **mat(NULL);
    mat = new float*[rows];

    for (int i=0;i<rows;i++)        
        mat[i] = new float[cols];
    return mat;
}

float *newVec(int n){               //alllouer un vecteur de dim n
    float* b;
    b  = new float [n];
    return b;
}

void displayMat(float** mat){    //afficher la matrice
    for (int i=0;i<dim;i++){
        for(int j=0; j<dim; j++)
            cout << mat[i][j] << " ";
        cout << endl;
    }
}

void displayVec (float* v){         //afficher un vecteur
    for(int i=0; i<dim; i++){
        cout << "[" << v[i] << "]" << endl;
    }
}


/*
 *
 *  DEFINITION : METHODE LDLT
 * 
 * 
*/

LDLT::LDLT(){
    // Allcation des matrice et des vecteurs
    A = newMat(dim,dim);
    L = newMat(dim,dim);
    D = newMat(dim,dim);
    b = newVec(dim);
    x = newVec(dim);
}

// Chargger la matrice depuis un fichier
void LDLT::loadAb(){
    string fileName;
    ifstream inFile;
    
    //cout << "Entrez le fichier a charger : "; cin >> fileName;

    // ouverture du fichier, 
    inFile.open("matrix.txt");

    if(inFile){
        //recuperer la taille de la matrice
        inFile >> dim;

        // recuperer la partie inférieure avec le duagonal de la matrice
        for(int i=0; i<dim;i++){
            for(int j=0; j<=i;j++){
                inFile >> A[i][j];
            }
        }  
        
        // recuperer le second membre
        for(int i=0; i<dim;i++){
                inFile >> b[i];
        }
    }
    
    // Alerter l'erreur d'ouverture du fichier
    else
        cout << "erreur de lecture" << endl;
}

// charger le profil p_i de A
void LDLT::computeP_i(){
    p_i = newVec(dim);

    for(int i=0;i<dim;i++){
        bool logicNull = true;
        for(int j=0;j<=i; j++){
            if(A[i][j]==0){
                if(logicNull)
                    cout <<"  ";
                else
                    cout <<". ";
            }
            else {
                cout << "x ";
                if(logicNull)
                    p_i[i] = j;        // enregistrer la valeur de j comme etant le premmier element non null de la ligne
                logicNull = false;     // les zeros qui suive ne sont plus null logiquement
            }
        }
        cout << endl;
    }
    cout << endl;
}

// determiner la valeur de L et D avec A = L.D.Lt
void LDLT::factrizeA(){
    float sum = 0;
    for(int i=0;i<dim;i++){
        L[i][i] = 1;

        // commencer les calcule à partir du la première élément non null
        for(int j=p_i[i];j<i;j++){
            sum=0;
            for(int k=p_i[i];k<=j-1;k++){
                // eviter de depenser le calcule pour des zero
                if(L[i][k]!=0)
                    sum += L[i][k]*D[k][k]*L[j][k];
            }

            // L conserve le profil de A
            // if(A[i][j]==0)              // 
            //     L[i][j] = 0;
            // else
                L[i][j] = (A[i][j] - sum)/D[j][j];
        }

            sum = 0;

        for(int k = p_i[i];k<=(i-1);k++){
            // eviter de depenser le calcule pour des zero
            if(L[i][k]!=0)
                sum += L[i][k]*D[k][k]*L[i][k];
        }
        D[i][i] = A[i][i] - sum;
    }
}

// resolution d'un system a diagonal inferieur
void LDLT::solveTriangInf(float **mat, float *vec){
    float s(0);
    int i(0), j(0);
    for(i=0; i<dim; i++){    /// de haut en bas
        for(j=0, s=0; j<i; j++)
            s += (mat[i][j]*x[j]);
        x[i] = (vec[i]-s)/mat[i][i];
    }

}

// methode de transposition
void LDLT::transpose(float **mat){
    for(int i = 0; i<dim; i++)
        for (int j = 0; j<dim; j++){
            mat[i][j] = mat[j][i];
        }  
}

// resolution d'un system a diagonal superieur
void LDLT::solveTriangSup(float **mat, float *vec){
    float s(0);
    int i(0), j(0);
    for(i=dim-1; i>=0; i--){    /// de bas en haut
        for(j=i+1, s=0; j<int(dim); j++)
            s += (mat[i][j]*x[j]);
        x[i] = (vec[i]-s)/mat[i][i];
    }
}

// Résoudre sans faire de transpose
void LDLT::solveTriangSupNoTranspose(float **mat, float *vec){
    float s(0);
    int i(0), j(0);
    for(i=dim-1; i>=0; i--){        // de bas en haut
        for(j=i+1, s=0; j<int(dim); j++)
            s += (mat[j][i]*x[j]);  // récupérer la symétrie au lieu de mat[i][j]
        x[i] = (vec[i]-s)/mat[i][i];
    }
}

void LDLT::solve(){
    // charger les donner dempuis les fichier
    loadAb();
    
    // Interface
    cout << "La matrice A :" << endl;
    displayMat(A);
    
    cout << "\nLe second membre :" << endl;
    displayVec(b);

    cout << "\nProfil de A :" << endl;
    computeP_i();         // tracer le profil

    
    factrizeA();          // calculer les valeurs de L et D 
    
    // A et L n bien le mem prfil
    profil(A,"A");
    profil(L,"L");


    solveTriangInf(L, b); //on resout L.x = b
    solveTriangInf(D, x); // on resout D.x = x

    //transpose(L); // on transpose L pour avoir Lt
    //solveTriangSup(L, x); // on resout Lt.x = x

    solveTriangSupNoTranspose(L, x); // Résoudre Lt.x = x sans passer par le transpose
    cout << "La solution du systeme est :"<< endl;
    displayVec(x);                   // afficher le resultat
}

void LDLT::profil(float** mat, string name){
    cout << "Profil de la matrice " << name << " :"<< endl;
    for(int i=0;i<dim;i++){
        // eviter de faire un boucle de plus
        bool logicNull = true;
        for(int j=0;j<=i; j++){
            if(mat[i][j]==0){
                if(logicNull)
                    continue;
                else
                    cout <<"("<< i+1 << ";" << j+1 << "), ";
            }
            else {
                cout <<"("<< i+1 << ";" << j+1 << "), ";
                logicNull = false; // les zeros qui suive ne sont plus null logiquement
            }
        }
        cout << endl;
    }
    cout << endl;
}

Lsolver::~Lsolver(){
	delete[] x;
	delete[] b;

	for(int i=0; i< dim; i++){
        delete[] D[i];
        delete[] A[i];
        delete[] L[i];
    }
    delete[] D;
    delete[] A;
    delete[] L;
}