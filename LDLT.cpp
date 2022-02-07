/*
 *  RATEFIARISON 
 *  Harivony Lalatiana
 * 
 *  Factorisation LDLt
 * 
 */

#include <iostream>
#include <fstream>
#include <string>

using namespace std;


int dim(10);

float **newMat(int rows, int cols);
float *newVec(int n);
void displayMat(float** mat);
void displayVec (float* v);

class LDLT{
// Methods
public:
    LDLT();

    void loadAb();
    void factrizeA();

    // getter
    float** getA(){return A;}
    float* getb(){return b;}
    float** getD(){return D;}
    float* getX(){return x;}

    void solveTriangInf(float **mat, float *vec); 
    void transpose(float **mat);
	void solveTriangSup(float **mat, float *vec);

    void solve();

// Attributs
    float** A;
    float** L;
    float** D;
    float* x;   // la solution

    float* b;   // second membre
};

int main (){
    LDLT solver;
    solver.solve();

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


/***** LDLT ******/

LDLT::LDLT(){
    A = newMat(dim,dim);
    L = newMat(dim,dim);
    D = newMat(dim,dim);
    b = newVec(dim);
    x = newVec(dim);
    for(int i=0;i<dim;i++){
        
        b[i] = 0;
        x[i] = 0;

        for(int j=0;j<dim;j++){
            A[i][j] = 0;
            L[i][j] = 0;
            D[i][j] = 0;
        }
    }
}

void LDLT::loadAb(){
    string fileName;
    ifstream inFile;
    
    //cout << "Entrez le fichier a charger : "; cin >> fileName;

    inFile.open("matrix.txt");
    if(inFile){
        inFile >> dim;

        for(int i=0; i<dim;i++){
            for(int j=0; j<dim;j++){
                inFile >> A[i][j];
            }
        }  
        
        for(int i=0; i<dim;i++){
                inFile >> b[i];
        }
    }
    else
        cout << "erreur de lecture" << endl;
}

void LDLT::factrizeA(){
    float sum = 0;
    for(int i=0;i<dim;i++){
        L[i][i] = 1;

        for(int j=0;j<i;j++){
            sum=0;
            for(int k=0;k<=j-1;k++){
                sum += L[i][k]*D[k][k]*L[j][k];
            }
            L[i][j] = (A[i][j] - sum)/D[j][j];
            // get D[i]
        }

            sum = 0;

            for(int k=0;k<=(i-1);k++){
                sum += L[i][k]*D[k][k]*L[i][k];
            }
            D[i][i] = A[i][i] - sum;
    }
}

void LDLT::solveTriangInf(float **mat, float *vec){
    float s(0);
    int i(0), j(0);
    for(i=0; i<dim; i++){    /// de haut en bas
        for(j=0, s=0; j<i; j++)
            s += (mat[i][j]*x[j]);
        x[i] = (vec[i]-s)/mat[i][i];
    }

}

void LDLT::transpose(float **mat){
    for(int i = 0; i<dim; i++)
        for (int j = 0; j<dim; j++){
            mat[i][j] = mat[j][i];
        }  
}

void LDLT::solveTriangSup(float **mat, float *vec){
    float s(0);
    int i(0), j(0);
    for(i=dim-1; i>=0; i--){    /// de bas en haut
        for(j=i+1, s=0; j<int(dim); j++)
            s += (mat[i][j]*x[j]);
        x[i] = (vec[i]-s)/mat[i][i];
    }
}

void LDLT::solve(){
    loadAb();
    factrizeA();

    // computeZ();
    // computeY();
    // computeX();
    solveTriangInf(L, b); //on resout L.x = b
    solveTriangInf(D, x); // on resout D.x = x
    transpose(L); // on transpose L pour avoir Lt
    solveTriangSup(L, x); // on resout Lt.x = x
    displayVec(x);
}
