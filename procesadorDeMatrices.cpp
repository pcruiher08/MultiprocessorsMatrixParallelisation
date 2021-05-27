#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;

void multiplicaDosMatrices(double** resultado, double** A, double** B, int verticalResultado, int horizontalResultado, int verticalA, int horizontalA, int verticalB, int horizontalB){


    for(int i = 0; i < verticalResultado; i++){
        for(int j = 0; j < horizontalResultado; j++){
            for(int k = 0; k < horizontalA; k++){
                //cout<< *(*(A+i)+k) << " * " << *(*(B+j)+k) << " = " << *(*(A+i)+k) * *(*(B+j)+k)<<endl;
                *(*(resultado+i)+j) += *(*(A+i)+k) * *(*(B+j)+k);
            }
        }
    }

}

void imprimeMatriz(double** matriz, int vertical, int horizontal){
    for(int i = 0; i < vertical; i++){
        for(int j = 0; j < horizontal; j++){
            //cout<<matriz[i][j]<<" ";
            cout<<*(*(matriz+i)+j)<<" ";
            //cout<<*( matriz + i * horizontal + j) << " ";
        }
        cout<<endl;
    }
}

int main(){
clock_t start, end;
int tiempoTranscurrido = 0;
ifstream lecturaA;
ifstream lecturaB;
lecturaA.open("matrixA1048576.txt");
lecturaB.open("matrixB1048576.txt");

int horizontalA, verticalA; 
cout<<"Cuanto mide la matriz A horizontalmente? "; 
cin>>horizontalA;
cout<<"Cuanto mide la matriz A verticalmente? ";
cin>>verticalA;

int horizontalB, verticalB; 
cout<<"Cuanto mide la matriz B horizontalmente? "; 
cin>>horizontalB;
cout<<"Cuanto mide la matriz B verticalmente? ";
cin>>verticalB;

if(horizontalA != verticalB){
    cout<<"No se pueden multiplicar las matrices"<<endl;
    return 0;
}

//se separa memoria para las matrices y se leen

//double matrizA[verticalA][horizontalA];
//double matrizB[verticalB][horizontalB];
double **matrizA = (double**)malloc(horizontalA * sizeof(double));
for(int i = 0; i < horizontalA; i++){
    matrizA[i] = (double*)malloc(verticalA * sizeof(double));
}

double **matrizB = (double**)malloc(horizontalB * sizeof(double));
for(int i = 0; i < horizontalB; i++){
    matrizB[i] = (double*)malloc(verticalB * sizeof(double));
}

//se prepara la lectura para la matriz resultante

int verticalResultado, horizontalResultado;
verticalResultado = verticalA;
horizontalResultado = horizontalB;

double **resultado = (double**)malloc(horizontalResultado * sizeof(double));
for(int i = 0; i < horizontalResultado; i++){
    resultado[i] = (double*)malloc(verticalResultado * sizeof(double));
}

double read; 

//se lee normal sin transponer
for(int i = 0; i < verticalA; i++){
    for(int j = 0; j < horizontalA; j++){
        lecturaA>>read;
        matrizA[i][j] = read;
        //*(matrizA + i * horizontalA + j) = read;
    }
}

//se lee transpuesta para hacer operaciones con memoria contigua
for(int i = 0; i < verticalB; i++){
    for(int j = 0; j < horizontalB; j++){
        lecturaB>>read;
        matrizB[j][i] = read;
        //*(matrizB + j * verticalA + i) = read;
    }
}

swap(verticalB, horizontalB);

//imprimeMatriz((double **)matrizA, verticalA, horizontalA);
//se imprime transpuesta tambien
//imprimeMatriz((double **)matrizB, verticalB, horizontalB);
//cout<<"vertical "<<verticalResultado << " horizontal "<<horizontalResultado<<endl;

start = clock();
multiplicaDosMatrices(resultado, matrizA, matrizB, verticalResultado, horizontalResultado, verticalA, horizontalA, verticalB, horizontalB);
end = clock();
tiempoTranscurrido = end - start;
//imprimeMatriz((double **)resultado, verticalResultado, horizontalResultado);
cout<<"Se obtuvo en "<<tiempoTranscurrido<<" ms"<<endl;

/*
for(int i = 0; i < vertical; i++){
    for(int j = 0; j < horizontal; j++){
        cout<<matriz[i][j]<<" ";
    }
    cout<<endl;
}
*/

lecturaA.close();
lecturaB.close();



return 0;
}