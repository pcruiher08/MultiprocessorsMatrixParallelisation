#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <chrono>
#include <iomanip>
#include <stdint.h>
#include <intrin.h>

using namespace std;

void multiplicaDosMatrices(double* resultado, double* A, double* B, int verticalResultado, int horizontalResultado, int verticalA, int horizontalA, int verticalB, int horizontalB) {
	for (int i = 0; i < verticalResultado; i++) {
		for (int j = 0; j < horizontalResultado; j++) {
			//*(*(resultado + i) + j) = 0;
			resultado[i * verticalResultado + j] = 0;
			for (int k = 0; k < horizontalA; k++) {
				//*(*(resultado + i) + j) += (*(*(A + i) + k)) * (*(*(B + j) + k));
				resultado[i * verticalResultado + j] += A[i * verticalResultado + k] * B[j * verticalResultado + k];

			}
		}
	}
}


void multiplicaDosMatricesOMP(double* resultado, double* A, double* B, int verticalResultado, int horizontalResultado, int verticalA, int horizontalA, int verticalB, int horizontalB) {
#pragma omp parallel for
	for (int i = 0; i < verticalResultado; i++) {
		for (int j = 0; j < horizontalResultado; j++) {
			resultado[i * verticalResultado + j] = 0;
			//*(*(resultado + i) + j) = 0;
			for (int k = 0; k < horizontalA; k++) {
				//*(*(resultado + i) + j) += *(*(A + i) + k) * *(*(B + j) + k);
				resultado[i * verticalResultado + j] += A[i * verticalResultado + k] * B[j * verticalResultado + k];

			}
		}
	}
}

void multiplicaDosMatricesIntrin(double* resultado, double* A, double* B, int verticalResultado, int horizontalResultado, int verticalA, int horizontalA, int verticalB, int horizontalB) {
	for (int i = 0; i < verticalResultado; i++) {
		for (int j = 0; j < horizontalResultado; j++) {
			resultado[i * verticalResultado + j] = 0;

			//*(*(resultado + i) + j) = 0;
			__m256d a_reg, b_reg, c_reg;
			double* aux;
			aux = (double*)malloc(sizeof(double) * 4);
			for (int k = 0; k < horizontalA; k += 4) {
				/*
				a_reg = _mm256_load_pd(*(A + i) + k);
				b_reg = _mm256_load_pd(*(B + j) + k);
				c_reg = _mm256_mul_pd(a_reg, b_reg);
				_mm256_store_pd(aux, c_reg);
				*(*(resultado + i) + j) += aux[0] + aux[1] + aux[2] + aux[3];
				*/
				a_reg = _mm256_load_pd(&A[i * verticalResultado + k]);
				b_reg = _mm256_load_pd(&B[j * verticalResultado + k]);
				c_reg = _mm256_mul_pd(a_reg, b_reg);
				_mm256_store_pd(aux, c_reg);

				resultado[i * verticalResultado + j] += aux[0] + aux[1] + aux[2] + aux[3];


			}
			free(aux);
		}
	}
}

void imprimeMatriz(double* matriz, int vertical, int horizontal, ofstream& archivoResultante) {
	for (int i = 0; i < vertical; i++) {
		for (int j = 0; j < horizontal; j++) {
			archivoResultante << fixed << setprecision(10) << matriz[i * vertical + j] << endl;
			//cout << matriz[i * vertical + j] << " ";
			//cout<<matriz[i][j]<<" ";
			//cout<<*(*(matriz+i)+j)<<" ";
			//archivoResultante << fixed << setprecision(10) << *(*(matriz + i) + j) << endl;

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
	int tiempoTranscurridoS1 = 0;
	uint64_t tiempoEnNanosegundosS1 = 0;
	int tiempoTranscurridoS2 = 0;
	uint64_t tiempoEnNanosegundosS2 = 0;
	int tiempoTranscurridoS3 = 0;
	uint64_t tiempoEnNanosegundosS3 = 0;
	int tiempoTranscurridoS4 = 0;
	uint64_t tiempoEnNanosegundosS4 = 0;
	int tiempoTranscurridoS5 = 0;
	uint64_t tiempoEnNanosegundosS5 = 0;
	int tiempoTranscurridoO1 = 0;
	uint64_t tiempoEnNanosegundosO1 = 0;
	int tiempoTranscurridoO2 = 0;
	uint64_t tiempoEnNanosegundosO2 = 0;
	int tiempoTranscurridoO3 = 0;
	uint64_t tiempoEnNanosegundosO3 = 0;
	int tiempoTranscurridoO4 = 0;
	uint64_t tiempoEnNanosegundosO4 = 0;
	int tiempoTranscurridoO5 = 0;
	uint64_t tiempoEnNanosegundosO5 = 0;
	int tiempoTranscurridoI1 = 0;
	uint64_t tiempoEnNanosegundosI1 = 0;
	int tiempoTranscurridoI2 = 0;
	uint64_t tiempoEnNanosegundosI2 = 0;
	int tiempoTranscurridoI3 = 0;
	uint64_t tiempoEnNanosegundosI3 = 0;
	int tiempoTranscurridoI4 = 0;
	uint64_t tiempoEnNanosegundosI4 = 0;
	int tiempoTranscurridoI5 = 0;
	uint64_t tiempoEnNanosegundosI5 = 0;
	ifstream lecturaA;
	ifstream lecturaB;
	ofstream archivoResultante;
	//lecturaA.open("matrixA16.txt");
	//lecturaB.open("matrixB16.txt");
	lecturaA.open("matrixA1048576.txt");
	lecturaB.open("matrixB1048576.txt");

	//lecturaA.open("matrixA9.txt");
	//lecturaB.open("matrixB9.txt");

	//archivoResultante.open("matrizResultante4.txt");
	archivoResultante.open("matrizResultante1024.txt");
	//archivoResultante.open("matrizChiquita.txt");
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
	/*
	double** matrizA = (double**)malloc(horizontalA * sizeof(double));
	for (int i = 0; i < horizontalA; i++) {
		matrizA[i] = (double*)malloc(verticalA * sizeof(double));
	}
	double** matrizB = (double**)malloc(horizontalB * sizeof(double));
	for (int i = 0; i < horizontalB; i++) {
		matrizB[i] = (double*)malloc(verticalB * sizeof(double));
	}
	*/

	double* matrizA = (double*)malloc(horizontalA * verticalA * sizeof(double));
	double* matrizB = (double*)malloc(horizontalA * verticalB * sizeof(double));


	//se prepara la lectura para la matriz resultante

	int verticalResultado, horizontalResultado;
	verticalResultado = verticalA;
	horizontalResultado = horizontalB;
	/*
	double** resultado = (double**)malloc(horizontalResultado * sizeof(double));
	for (int i = 0; i < horizontalResultado; i++) {
		resultado[i] = (double*)malloc(verticalResultado * sizeof(double));
	}
	*/
	double* resultadoS = (double*)malloc(horizontalResultado * verticalResultado * sizeof(double));
	double* resultadoO = (double*)malloc(horizontalResultado * verticalResultado * sizeof(double));
	double* resultadoI = (double*)malloc(horizontalResultado * verticalResultado * sizeof(double));


	double read;

	//se lee normal sin transponer
	for (int i = 0; i < verticalA; i++) {
		for (int j = 0; j < horizontalA; j++) {
			lecturaA >> read;
			matrizA[i * verticalA + j] = read;
			//matrizA[i][j] = read;
			//cout << read << endl;
			//*(matrizA + i * horizontalA + j) = read;
		}
	}

	//se lee transpuesta para hacer operaciones con memoria contigua
	for (int i = 0; i < verticalB; i++) {
		for (int j = 0; j < horizontalB; j++) {
			lecturaB >> read;
			matrizB[j * verticalB + i] = read; // hay que revisar que se este transponiendo con una matriz pequenia
			//matrizB[j][i] = read;
			//*(matrizB + j * verticalA + i) = read;
		}
	}

	swap(verticalB, horizontalB);
	//imprimeMatriz(matrizA, verticalA, horizontalA, archivoResultante);
	//imprimeMatriz(matrizB, verticalB, horizontalB, archivoResultante);


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

	for (int i = 0; i < 5; i++) {

		cout << "-------Iteracion " << i + 1 << "-------" << endl;

		cout << "-----Serial-----" << endl;

		start = clock();
		inicio = nanos();

		multiplicaDosMatrices(resultadoS, matrizA, matrizB, verticalResultado, horizontalResultado, verticalA, horizontalA, verticalB, horizontalB);

		fin = nanos();
		end = clock();

		tiempoTranscurrido = end - start;
		tiempoEnNanosegundos = fin - inicio;

		switch (i) {
			case 0:
				tiempoTranscurridoS1 = end - start;
				tiempoEnNanosegundosS1 = fin - inicio;
			case 1:
				tiempoTranscurridoS2 = end - start;
				tiempoEnNanosegundosS2 = fin - inicio;
			case 2:
				tiempoTranscurridoS3 = end - start;
				tiempoEnNanosegundosS3 = fin - inicio;
			case 3:
				tiempoTranscurridoS4 = end - start;
				tiempoEnNanosegundosS4 = fin - inicio;
			case 4:
				tiempoTranscurridoS5 = end - start;
				tiempoEnNanosegundosS5 = fin - inicio;
		}

		cout << "Se obtuvo el resultado en " << tiempoTranscurrido << " ms" << endl;
		cout << "Se obtuvo el resultado en " << tiempoEnNanosegundos << " ns" << endl;

		//imprimeMatriz(resultado, verticalResultado, horizontalResultado, archivoResultante);

		cout << "-----Open MP-----" << endl;

		start = clock();
		inicio = nanos();

		multiplicaDosMatricesOMP(resultadoO, matrizA, matrizB, verticalResultado, horizontalResultado, verticalA, horizontalA, verticalB, horizontalB);

		fin = nanos();
		end = clock();

		tiempoTranscurrido = end - start;
		tiempoEnNanosegundos = fin - inicio;

		switch (i) {
		case 0:
			tiempoTranscurridoO1 = end - start;
			tiempoEnNanosegundosO1 = fin - inicio;
		case 1:
			tiempoTranscurridoO2 = end - start;
			tiempoEnNanosegundosO2 = fin - inicio;
		case 2:
			tiempoTranscurridoO3 = end - start;
			tiempoEnNanosegundosO3 = fin - inicio;
		case 3:
			tiempoTranscurridoO4 = end - start;
			tiempoEnNanosegundosO4 = fin - inicio;
		case 4:
			tiempoTranscurridoO5 = end - start;
			tiempoEnNanosegundosO5 = fin - inicio;
		}

		cout << "Se obtuvo el resultado en " << tiempoTranscurrido << " ms" << endl;
		cout << "Se obtuvo el resultado en " << tiempoEnNanosegundos << " ns" << endl;

		cout << "-----Intrinsecas-----" << endl;

		start = clock();
		inicio = nanos();

		multiplicaDosMatricesIntrin(resultadoI, matrizA, matrizB, verticalResultado, horizontalResultado, verticalA, horizontalA, verticalB, horizontalB);

		fin = nanos();
		end = clock();

		tiempoTranscurrido = end - start;
		tiempoEnNanosegundos = fin - inicio;

		switch (i) {
		case 0:
			tiempoTranscurridoI1 = end - start;
			tiempoEnNanosegundosI1 = fin - inicio;
		case 1:
			tiempoTranscurridoI2 = end - start;
			tiempoEnNanosegundosI2 = fin - inicio;
		case 2:
			tiempoTranscurridoI3 = end - start;
			tiempoEnNanosegundosI3 = fin - inicio;
		case 3:
			tiempoTranscurridoI4 = end - start;
			tiempoEnNanosegundosI4 = fin - inicio;
		case 4:
			tiempoTranscurridoI5 = end - start;
			tiempoEnNanosegundosI5 = fin - inicio;
		}

		cout << "Se obtuvo el resultado en " << tiempoTranscurrido << " ms" << endl;
		cout << "Se obtuvo el resultado en " << tiempoEnNanosegundos << " ns" << endl;

		cout << "-------Comprobacion-------" << endl;

		//TODO: check if it is correct
		for (int j = 0; j < verticalResultado; j++) {
			for (int k = 0; k < horizontalResultado; k++) {
				if (resultadoS[j * verticalResultado + k] != resultadoO[j * verticalResultado + k] || resultadoS[j * verticalResultado + k] != resultadoI[j * verticalResultado + k]) {
					cout << "Error en alguna de las matrices" << endl;
					return 0;
				}
			}
		}

		cout << "Las matrices fueron calculadas correctamente" << endl;

	}

	cout << "----------Guardando el resultado----------" << endl;

	start = clock();
	inicio = nanos();

	imprimeMatriz(resultadoS, verticalResultado, horizontalResultado, archivoResultante);

	fin = nanos();
	end = clock();

	tiempoTranscurrido = end - start;
	tiempoEnNanosegundos = fin - inicio;

	cout << "Se guardo el archivo en " << tiempoTranscurrido << " ms" << endl;
	cout << "Se guardo el archivo en " << tiempoEnNanosegundos << " ns" << endl;

	lecturaA.close();
	lecturaB.close();
	archivoResultante.close();
	/*
	for (int i = 0; i < horizontalA; i++) {
		free(matrizA[i]);
	}
	free(matrizA);
	for (int i = 0; i < horizontalA; i++) {
		free(matrizB[i]);
	}
	free(matrizB);
	for (int i = 0; i < horizontalA; i++) {
		free(resultado[i]);
	}
	free(resultado);
	*/
	free(matrizA);
	free(matrizB);
	free(resultadoS);
	free(resultadoO);
	free(resultadoI);
	return 0;
}
