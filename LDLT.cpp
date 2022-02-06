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
    void computeZ();
    void computeY();
    void computeX();
    void solve();

// Attributs
    float** A;
    float** L;
    float** Lt;
    float* D;
    float* x;   // la solution
    float* y;   
    float* z;

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
    Lt = newMat(dim,dim);
    D = newVec(dim);
    b = newVec(dim);
    x = newVec(dim);
    y = newVec(dim);
    z = newVec(dim);
    for(int i=0;i<dim;i++){
        D[i] = 0;
        b[i] = 0;
        x[i] = 0;
        y[i] = 0;
        z[i] = 0;

        for(int j=0;j<dim;j++){
            A[i][j] = 0;
            L[i][j] = 0;
            Lt[i][j] = 0;
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
        Lt[i][i] = 1;

        for(int j=0;j<=i;j++){
            sum=0;
            for(int k=0;k<j-1;k++){
                sum += D[k]*L[j][k]*L[j][k];
            }
            
            // get D[i]
            D[j] = A[j][j] - sum;
            //cout << D[j] << " ";

            // get L[i][j]
            sum = 0;
            for(int k=0;k<(j);k++){
                sum += L[i][k]*D[k]*L[j][k];
            }

            L[i][j] = (A[i][j] - sum)/D[j];
            //cout << L[i][j] << " ";
        }
        //cout << endl;
    }
}

void LDLT::computeZ(){
    float sum =0;
    for(int i=0;i<dim;i++){
        sum = 0;
        for(int k=0;k<i-1;k++){
            sum+= L[i][k]*z[k];
        }
        z[i] = b[i] - sum;
    }
} 

void LDLT::computeY(){
    for(int i=0;i<dim;i++){
        y[i] = z[i]/D[i];
    }
}

void LDLT::computeX(){
    float sum = 0;
    for (int i=dim-1;i>=0;i--){
        for (int k=i;k>=0;k--){
            sum += L[k][i]*x[k];
        }
        x[i] = y[i] - sum;
    }
}

void LDLT::solve(){
    loadAb();
    factrizeA();
    computeZ();
    computeY();
    computeX();
    displayVec(x);
}

    //
    /*
    {
     solver.loadAb();
    // cout << "La matrice A :" << endl;
    // displayMat(solver.A);
    // cout << "Le vecteur b :" << endl;
    // displayVec(solver.b);

    // // solver.factrizeA();
    // // solver.computeZ();

    // // cout << " D :" << endl;
    // // // displayVec(solver.D);

    // // cout << " z :" << endl;
    // // displayVec(solver.z);

    // solver.computeY();
    // // cout << " y :" << endl;
    // // displayVec(solver.y);

    // solver.computeX();
    // cout << " x :" << endl;
    // displayVec(solver.x);
    // // cout << endl;
    // // displayMat(solver.L);
    // // cout << endl;
    // // displayMat(solver.Lt);
    // 
    }*/