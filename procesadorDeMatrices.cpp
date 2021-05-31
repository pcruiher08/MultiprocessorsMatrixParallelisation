#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <chrono>
#include <iomanip>
#include <stdint.h>

using namespace std;

void multiplicaDosMatrices(double** resultado, double** A, double** B, int verticalResultado, int horizontalResultado, int verticalA, int horizontalA, int verticalB, int horizontalB) {
    for (int i = 0; i < verticalResultado; i++) {
        for (int j = 0; j < horizontalResultado; j++) {
            *(*(resultado + i) + j) = 0;
            for (int k = 0; k < horizontalA; k++) {
                //cout<< *(*(A+i)+k) << " * " << *(*(B+j)+k) << " = " << *(*(A+i)+k) * *(*(B+j)+k)<<endl;
                *(*(resultado + i) + j) += (*(*(A + i) + k)) * (*(*(B + j) + k));
            }
        }
    }
}


void multiplicaDosMatricesOMP(double** resultado, double** A, double** B, int verticalResultado, int horizontalResultado, int verticalA, int horizontalA, int verticalB, int horizontalB) {
#pragma omp parallel for
    for (int i = 0; i < verticalResultado; i++) {
        for (int j = 0; j < horizontalResultado; j++) {
            *(*(resultado + i) + j) = 0;
//#pragma omp parallel for firstprivate(i, j)
            for (int k = 0; k < horizontalA; k++) {
                //cout<< *(*(A+i)+k) << " * " << *(*(B+j)+k) << " = " << *(*(A+i)+k) * *(*(B+j)+k)<<endl;
                *(*(resultado + i) + j) += *(*(A + i) + k) * *(*(B + j) + k);
            }
        }
    }
}

void imprimeMatriz(double** matriz, int vertical, int horizontal, ofstream& archivoResultante) {
    for (int i = 0; i < vertical; i++) {
        for (int j = 0; j < horizontal; j++) {
            //cout<<matriz[i][j]<<" ";
            //cout<<*(*(matriz+i)+j)<<" ";
            archivoResultante << fixed << setprecision(10) << *(*(matriz + i) + j) << endl;

            //cout<<*( matriz + i * horizontal + j) << " ";
        }
        //cout<<endl;
    }
}

uint64_t nanos()
{
    uint64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::
        now().time_since_epoch()).count();
    return ns;
}

int main() {
    clock_t start, end;
    uint64_t inicio, fin;
    int tiempoTranscurrido = 0;
    uint64_t tiempoEnNanosegundos = 0;
    ifstream lecturaA;
    ifstream lecturaB;
    ofstream archivoResultante;
    //lecturaA.open("matrixA9.txt");
    //lecturaB.open("matrixB9.txt");
    lecturaA.open("matrixA1048576.txt");
    lecturaB.open("matrixB1048576.txt");
    //archivoResultante.open("matrizResultante3.txt");
    archivoResultante.open("matrizResultante1024.txt");

    int horizontalA, verticalA;
    cout << "Cuanto mide la matriz A horizontalmente? ";
    cin >> horizontalA;
    cout << "Cuanto mide la matriz A verticalmente? ";
    cin >> verticalA;

    int horizontalB, verticalB;
    cout << "Cuanto mide la matriz B horizontalmente? ";
    cin >> horizontalB;
    cout << "Cuanto mide la matriz B verticalmente? ";
    cin >> verticalB;

    if (horizontalA != verticalB) {
        cout << "No se pueden multiplicar las matrices" << endl;
        return 0;
    }

    //se separa memoria para las matrices y se leen

    //double matrizA[verticalA][horizontalA];
    //double matrizB[verticalB][horizontalB];

    cout << "----------Obtencion de las matrices----------" << endl;

    start = clock();
    inicio = nanos();

    double** matrizA = (double**)malloc(horizontalA * sizeof(double));
    for (int i = 0; i < horizontalA; i++) {
        matrizA[i] = (double*)malloc(verticalA * sizeof(double));
    }

    double** matrizB = (double**)malloc(horizontalB * sizeof(double));
    for (int i = 0; i < horizontalB; i++) {
        matrizB[i] = (double*)malloc(verticalB * sizeof(double));
    }

    //se prepara la lectura para la matriz resultante

    int verticalResultado, horizontalResultado;
    verticalResultado = verticalA;
    horizontalResultado = horizontalB;

    double** resultado = (double**)malloc(horizontalResultado * sizeof(double));
    for (int i = 0; i < horizontalResultado; i++) {
        resultado[i] = (double*)malloc(verticalResultado * sizeof(double));
    }

    double read;

    //se lee normal sin transponer
    for (int i = 0; i < verticalA; i++) {
        for (int j = 0; j < horizontalA; j++) {
            lecturaA >> read;
            matrizA[i][j] = read;
            //cout << read << endl;
            //*(matrizA + i * horizontalA + j) = read;
        }
    }

    //se lee transpuesta para hacer operaciones con memoria contigua
    for (int i = 0; i < verticalB; i++) {
        for (int j = 0; j < horizontalB; j++) {
            lecturaB >> read;
            matrizB[j][i] = read;
            //*(matrizB + j * verticalA + i) = read;
        }
    }

    swap(verticalB, horizontalB);

    fin = nanos();
    end = clock();

    tiempoTranscurrido = end - start;
    tiempoEnNanosegundos = fin - inicio;

    cout << "Se configuiraron las matrices en " << tiempoTranscurrido << " ms" << endl;
    cout << "Se configuiraron las matrices en " << tiempoEnNanosegundos << " ns" << endl;


    //imprimeMatriz((double **)matrizA, verticalA, horizontalA);
    //se imprime transpuesta tambien
    //imprimeMatriz((double **)matrizB, verticalB, horizontalB);
    //cout<<"vertical "<<verticalResultado << " horizontal "<<horizontalResultado<<endl;

    cout << "----------Calculo de la matriz resultante----------" << endl;

    cout << "-----Serial-----" << endl;

    start = clock();
    inicio = nanos();

    multiplicaDosMatrices(resultado, matrizA, matrizB, verticalResultado, horizontalResultado, verticalA, horizontalA, verticalB, horizontalB);

    fin = nanos();
    end = clock();

    tiempoTranscurrido = end - start;
    tiempoEnNanosegundos = fin - inicio;

    cout << "Se obtuvo el resultado en " << tiempoTranscurrido << " ms" << endl;
    cout << "Se obtuvo el resultado en " << tiempoEnNanosegundos << " ns" << endl;

    //imprimeMatriz((double**)resultado, verticalResultado, horizontalResultado, archivoResultante);
    
    cout << "-----Open MP-----" << endl;

    start = clock();
    inicio = nanos();

    multiplicaDosMatricesOMP(resultado, matrizA, matrizB, verticalResultado, horizontalResultado, verticalA, horizontalA, verticalB, horizontalB);

    fin = nanos();
    end = clock();

    tiempoTranscurrido = end - start;
    tiempoEnNanosegundos = fin - inicio;

    cout << "Se obtuvo el resultado en " << tiempoTranscurrido << " ms" << endl;
    cout << "Se obtuvo el resultado en " << tiempoEnNanosegundos << " ns" << endl;
    
    cout << "----------Guardando el resultado----------" << endl;

    start = clock();
    inicio = nanos();

    imprimeMatriz((double**)resultado, verticalResultado, horizontalResultado, archivoResultante);

    fin = nanos();
    end = clock();

    tiempoTranscurrido = end - start;
    tiempoEnNanosegundos = fin - inicio;

    cout << "Se guardo el archivo en " << tiempoTranscurrido << " ms" << endl;
    cout << "Se guardo el archivo en " << tiempoEnNanosegundos << " ns" << endl;


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
    archivoResultante.close();


    return 0;
}