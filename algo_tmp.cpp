#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>

using namespace std;

class Profil {
	public:
		Profil(string _filename);
		~Profil();
		double **newMat(int dim);
		void getData();
		void display(double ** Matrix, int lig, int col, string name);
		void display1D(double *Matrix, int taille, string title);
		void APi_nDiag(); 
		void VectorDouble(vector<double> & double_);
		void VectorInt(vector<int> & int_);
		void factorisation();
		double getAp(int i, int j);
		void setAp(int i, int j, double val);
		void solve();
		
		
		double **A, *b, *x;
		int dim;
		vector<double> AP = { };
			
	private:
		vector<int> nDiag = { };
		string filename;
		int *l, *p;
};

void Profil::solve(){
/// Solve L.z = b <=> L.x = b
    x = new double[dim];
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
    
/// Solve D.y = z <=> D.x = x
    for (int i(0); i < dim; i++)
    {
        x[i] = x[i] / AP[nDiag[i]];
    }
   
/// Solve Lt.x = y <=> Lt.x = y
    for (int i(dim - 1); i >= 0; i--)
    {
        s = 0;
        for (int j(dim - 1); j > i; j--)
        {
            s += getAp(j, i) * x[j];
        }
        x[i] = (x[i] - s);
    }
}

double Profil::getAp(int i, int j){
	
	if(j >= p[i] && j<=i){
		return AP[nDiag[i] - i + j];
	}	
	else 
	return 0;
}
void Profil::setAp(int i, int j, double val){
	if(j >= p[i]) AP[nDiag[i] - i + j] = val;
}
void Profil::factorisation()
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

void Profil::VectorDouble(vector<double> & double_){
	vector<double>::iterator it;
	
	for(it = double_.begin(); it != double_.end(); it++)    {
          cout<< *it <<" ";     
	}
}
void Profil::VectorInt(vector<int> & int_){
	vector<int>::iterator it;
	
	for(it = int_.begin(); it != int_.end(); it++)    {
          cout<< *it <<" ";     
	}
}

void Profil::APi_nDiag(){
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

void Profil::display1D(double *Matrix, int taille, string title){
	cout<<endl<<"__________"<<title <<"______________"<<endl;
	for(int i=0;i<taille;i++){
		cout<< Matrix[i] <<" ";
	}	
	cout<<endl;	
}
void Profil::display(double ** Matrix, int lig, int col, string name){
	cout<<"__________"<<name <<"______________"<<endl;
	for(int i=0;i<lig;i++){
		for(int j=0; j<col; j++){
			cout<< Matrix[i][j] <<"\t";			
		}
		cout<<endl;
	}		
}
void Profil::getData(){
	string temp;
	ifstream file(filename);
	if(file.is_open()){
		///take the dimension of the matrix
		file >> dim;
		getline(file, temp);
		///construct the matrix
		A = newMat(dim);
		
		///read the matrix
	
		for(int i = 0; i < dim; i++){
			for(int j = 0; j <dim ; j++)
			{				
				file >> A[i][j];				
			}
			getline(file, temp);
			
		}
		// to fill the second member
        getline(file, temp);
        b = new double[dim];
        for (int i = 0; i < dim; i++)
            file >> b[i];
	}
}

double** Profil::newMat(int dim){	
	double **M(NULL);
	M = new double*[dim];
	for(int i=0; i<dim; i++){
		M[i] = new double[dim];
	}
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			M[i][j] = 0;
		}
	}
	return M;
}

Profil::~Profil(){
	
}

Profil::Profil(string _filename){
	filename = _filename;
}
int main(){
	
	Profil profil("Matrix.txt");
	cout<<"This is the resolution of the system Ax= b"<<endl;
	profil.getData();
	
	profil.display(profil.A, profil.dim, profil.dim, "A");
	profil.display1D(profil.b, profil.dim, "b");
	cout<<endl <<"___________________________________________________"<<endl;
	
	cout<<endl<<"Building the AP from the matrix A"<<endl;
	profil.APi_nDiag();
	cout<<"the AP before the factorisation of the matrix A " <<endl;
	profil.VectorDouble(profil.AP);
	
	cout<<endl <<"___________________________________________________"<<endl;
	cout<<endl<<"The factorisation of the matrix "<<endl;
	profil.factorisation();
	
	cout<<"the AP after the factorisation of the matrix A " <<endl;
	profil.VectorDouble(profil.AP);
	
	cout<<endl <<"___________________________________________________"<<endl;
	cout<<endl<<"Solving the matrix from AP factorised"<<endl;
	profil.solve();
	cout<<"The solution of the system"<<endl;
	profil.display1D(profil.x, profil.dim, "x");

	return 0;
}
