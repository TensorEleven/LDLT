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
#include <vector>

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
    void APi_nDiag();
    void factrizeA();

    // getter
    float** getA(){return A;}
    float* getb(){return b;}
    float** getD(){return D;}
    float* getX(){return x;}
    float getAp(int i, int j);
    void setAp(int i, int j, float val);

    void solveTriangInf();
    void solveDiag();
    void transpose(float **mat);
	void solveTriangSup();
    void solveTriangSupNoTranspose(float **mat, float *vec);
    void profil(float** mat, string name);
    void getProfil(float* mat, string name);
    void solve();
    void gaussSolve();
    vector<float> AP = {};
	vector<int> nDiag = {};
	int *l, *p;

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
void LDLT::factrizeA()
{     
    float s(0);
    for (int i(0); i < dim; i++)
    {
        for (int j(0); j < i; j++)
        {
            s = 0;
            for (int k(0); k < j; k++){
                s+= getAp(i, k) * AP[nDiag[k]] * getAp(j,k);
            }
            setAp(i,j, (getAp(i,j) - s) / AP[nDiag[j]]);
        }
        s = 0;
        for (int k(0); k < i; k++){
            s+= AP[nDiag[k]] * getAp(i,k) * getAp(i,k);
        }
        AP[nDiag[i]] = AP[nDiag[i]] - s;
    }
}

// solve L.x = b
void LDLT::solveTriangInf(){
    float s(0);
    for (int i(0); i < dim; i++)
    {
        s = 0;
        for (int j(0); j < i; j++)
        {
            s += getAp(i,j) * x[j];
        }
        x[i] = (b[i] - s);
    }
}

// solve D.x = x
void LDLT::solveDiag(){
    for (int i(0); i < dim; i++){
        x[i] = x[i] / AP[nDiag[i]];
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
void LDLT::solveTriangSup(){
    float s(0);
    for (int i(dim - 1); i >= 0; i--){
        s = 0;
        for (int j(dim - 1); j > i; j--)
        {
            s += getAp(j, i) * x[j];
        }
        x[i] = (x[i] - s);
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

    APi_nDiag();    

    solveTriangInf();
    solveDiag();
    solveTriangSup();

    cout << "La solution du systeme est :"<< endl;
    displayVec(x);                   // afficher le resultat
}

// afficher le profi de la matrice mat
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

void LDLT::APi_nDiag(){
	int d(0);
	
	for(int i=0; i<dim; i++){
		for(int j=0; j<=i; j++){
			if(A[i][j] != 0){
				for(int k=j; k<=i; k++){
					 AP.push_back(A[i][k]);
					 d++;
					 if(k == i){
						nDiag.push_back(d - 1);
					 }					 
				}
				break;			
			}
		}
	}
	
	l = new int[dim];
	p = new int[dim];
	l[0] = 0;
	p[0] = 0;
	for(int i = 1; i<dim; i++){
		l[i] = nDiag[i] - nDiag[i-1] -1;
		p[i] = i - l[i] + 1 - 1;
	}
}

float LDLT::getAp(int i, int j){
	
	if(j >= p[i] && j<=i){
		return AP[nDiag[i] - i + j];
	}	
	else 
	return 0;
}

void LDLT::setAp(int i, int j, float val){
	if(j >= p[i]) AP[nDiag[i] - i + j] = val;
}

LDLT::~LDLT(){
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