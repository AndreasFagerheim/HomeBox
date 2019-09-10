#include <iostream>

using namespace std;
#include <cstdlib>
#include <cmath>

#include <vector>
#include <time.h>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <armadillo>
using namespace std;
ofstream ofile;



/* Generell metode for Gaussian ellimenering gitt tre vektorer for diagoanlene i en tridiagonal matrise */
void alg_tri_general(double *a_vec, double *b_vec, double *c_vec, double *ftilde, double *btilde, double *beta, double *f_vec, int n){

    beta[1] = b_vec[1];
    f_vec[1] = btilde[1];
    double temp;

 /* regner ut forward substitution*/
    for (int i = 2; i <n+1; i++){
        temp = a_vec[i-1]/beta[i-1];
        beta[i] = b_vec[i] - c_vec[i-1] * temp;
        f_vec[i] = btilde[i] - f_vec[i-1] * temp;
    }

    ftilde[n] = f_vec[n]/beta[n];
    /* regner ut backward substitution*/
       for (int i = 2; i <n+1; i++){
           ftilde[i] = (f_vec[i] - ftilde[i+1]*c_vec[i])/beta[i];
       }
}
/* Spesiell metode for Gaussian ellimenering hvor hver diagonale vektor tilsvarer en oppgitt verdi a = c=-1, b= 2*/
void alg_tri_special( double *ftilde, double *btilde, double *beta, double *f_vec, int n){

    f_vec[1] = btilde[1];
    double temp;

 /* regner ut forward substitution*/
    for (int i = 2; i <n+1; i++){
        temp = 1/beta[i-1];
        beta[i] = 2 -  temp;
        f_vec[i] = btilde[i] + f_vec[i-1] * temp;
    }

    ftilde[n] = f_vec[n]/beta[n];
    /* regner ut backward substitution*/
       for (int i = 2; i <n+1; i++){
           ftilde[i] = (f_vec[i] + ftilde[i+1])/beta[i];
       }
}


double fx(double x){
    return 100*exp(-10*x);
}
double ux(double x){
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}





void problem_a(int n){
    // lager vectorer med n+2 elementer med plass til boundery conditions
    double *a = new double[n+2];
    double *b = new double[n+2];
    double *c = new double[n+2];

    //høyresiden
    double *ftilde = new double[n+2];
    ftilde[0] = 0;
    ftilde[n+1] = 0;                    // boundery conditions
    double *btilde = new double[n+2];   //

    double *analytisk = new double[n+2];
   // double *numerisk = new double[n+2];

    double *x = new double[n+2];        // steg
    double h = 1.0/(n+1.0);             // steglengde
    double h2 = h*h;

    double *beta = new double[n+1];
    double *f_vec = new double[n+1];
    // fyller vector vi vet verdien til

    for(int i = 0;i<n;i++){
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
    }
    for(int i = 0;i<n+2;i++){
        x[i] = i *h;
        //cout<< x[i];
        btilde[i] = fx(x[i])*h2;        //definisjon fra oppgaven
        analytisk[i] = ux(x[i]);
        cout << analytisk[i];

    }
    double start = clock();
    alg_tri_general(a,b,c,ftilde,btilde,beta,f_vec, n);
    double stopp = clock();
    double time = (stopp-start);
    //cout <<"start"<<start <<"stopp=" <<stopp;
    //cout << "Tid: tridiagonal generell gaussian:"<<time;
    //cout << a;
    ofile.open("data1a.txt");
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << n << endl;
    ofile << "       x:             analytisk:          ftilde(x):  " << endl;
    for (int i=0;i<n+2;i++) {
       ofile << setw(15) << setprecision(8) << x[i];
       ofile << setw(15) << setprecision(8) << analytisk[i];
       ofile << setw(15) << setprecision(8) << ftilde[i] << endl;
    }

    ofile.close();

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] ftilde;
    delete [] beta;
    delete [] x;
    delete [] btilde;
    delete [] f_vec;

}
//---------------------------------------------------------------------------------------------------------------------------------------------
void problem_b(int n){
    // lager vectorer med n+2 elementer med plass til boundery conditions


    //høyresiden
    double *ftilde = new double[n+2];
    ftilde[0] = 0;
    ftilde[n+1] = 0;                    // boundery conditions
    double *btilde = new double[n+2];   //

    double *analytisk = new double[n+2];
   // double *numerisk = new double[n+2];

    double *x = new double[n+2];        // steg
    double h = 1.0/(n+1.0);             // steglengde
    double h2 = h*h;

    double *beta = new double[n+1];
    double *f_vec = new double[n+1];
    // fyller vector vi vet verdien til


    for(int i = 0;i<n+2;i++){
        x[i] = i *h;
        //cout<< x[i];
        btilde[i] = fx(x[i])*h2;        //definisjon fra oppgaven
        analytisk[i] = ux(x[i]);
        cout << analytisk[i];

    }
    double start = clock();
    alg_tri_special(ftilde,btilde,beta,f_vec, n);
    double stopp = clock();
    double time = (stopp-start);
    //cout <<"start"<<start <<"stopp=" <<stopp;
    //cout << "Tid: tridiagonal generell gaussian:"<<time;
    //cout << a;
    ofile.open("data1a.txt");
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << n << endl;
    ofile << "       x:             analytisk:          ftilde(x):  " << endl;
    for (int i=0;i<n+2;i++) {
       ofile << setw(15) << setprecision(8) << x[i];
       ofile << setw(15) << setprecision(8) << analytisk[i];
       ofile << setw(15) << setprecision(8) << ftilde[i] << endl;
    }

    ofile.close();

    delete [] ftilde;
    delete [] beta;
    delete [] x;
    delete [] btilde;
    delete [] f_vec;

}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void problem_e(int n){
    arma::mat A = arma::zeros<arma::mat>(n+2,n+2);
    arma::vec a(n+1);
    arma::vec b(n+2);
    arma::vec c(n+1);

    a.fill(-1);
    b.fill(2);
    c.fill(-1);
    arma::vec x(n+2);
    arma::vec btilde(n+2);
    double h = 1.0/(n+1);
    double h2 = h*h;
    for (int i = 0;i<n+2; i ++){
        x[i] = i*h;
        btilde[i] = h2*fx(x[i]);
    }
    // lager tridiagonale matrisen A

    for(int i= 0; i<n+2;i++){
        A[i,i] = b(i);
       if(i>0){
           A[i,i-1] = a(i);
       }
       if(i<n+2){
               A[i,i+1] = c(i);
           }
    }



    arma::vec v = arma::solve(A,btilde);
    arma::mat L;
    arma::mat U;
    arma::lu(L,U,A);

}

int main(int argc, char* args[]) {

    int n = 10;

    if (argc > 1) {
        n = atof(args[1]);
        cout<<n;
    }

    problem_a(n);

    return 0;
}
