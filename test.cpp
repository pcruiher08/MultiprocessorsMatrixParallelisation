#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <chrono>
#include <iomanip>
#include <stdint.h>
#include <intrin.h>
#include <cstdlib>
#include <vector>

using namespace std;

void liberaMemoria(vector<double*> espaciosDeMemoria) {
	for (int i = 0; i < espaciosDeMemoria.size(); i++) {
		free(espaciosDeMemoria[i]);
	}
	cout << "Se ha liberado la memoria" << endl;
}

void multiplicaDosMatrices(double* resultado, double* A, double* B, int verticalResultado, int horizontalResultado, int verticalA, int horizontalA, int verticalB, int horizontalB) {
	for (int i = 0; i < verticalResultado; i++) {
		for (int j = 0; j < horizontalResultado; j++) {
			resultado[i * verticalResultado + j] = 0;
			for (int k = 0; k < horizontalA; k++) {
				resultado[i * horizontalResultado + j] += A[i * horizontalResultado + k] * B[j * horizontalResultado + k];

			}
		}
	}
}


void multiplicaDosMatricesOMP(double* resultado, double* A, double* B, int verticalResultado, int horizontalResultado, int verticalA, int horizontalA, int verticalB, int horizontalB) {
#pragma omp parallel for
	for (int i = 0; i < verticalResultado; i++) {
		for (int j = 0; j < horizontalResultado; j++) {
			resultado[i * verticalResultado + j] = 0;
			for (int k = 0; k < horizontalA; k++) {
				resultado[i * verticalResultado + j] += A[i * verticalResultado + k] * B[j * verticalResultado + k];

			}
		}
	}
}

bool multiplicaDosMatricesIntrin(double* resultado, double* A, double* B, int verticalResultado, int horizontalResultado, int verticalA, int horizontalA, int verticalB, int horizontalB) {
	for (int i = 0; i < verticalResultado; i++) {
		for (int j = 0; j < horizontalResultado; j++) {
			resultado[i * verticalResultado + j] = 0;

			__m256d a_reg, b_reg, c_reg;
			double* aux;
			aux = (double*)malloc(sizeof(double) * 4);
			

			if (aux == NULL) {
				cout << "No se pudo separar memoria interna para el calculo de intrinsecas, se va a liberar la memoria previamente separada y terminara la ejecución del programa" << endl;
				return false;
			}

			aux[0] = 0;
			aux[1] = 0;
			aux[2] = 0;
			aux[3] = 0;

			for (int k = 0; k < horizontalA; k += 4) {
				a_reg = _mm256_load_pd(&A[i * verticalResultado + k]);
				b_reg = _mm256_load_pd(&B[j * verticalResultado + k]);
				c_reg = _mm256_mul_pd(a_reg, b_reg);
				_mm256_store_pd(aux, c_reg);

				resultado[i * verticalResultado + j] += aux[0] + aux[1] + aux[2] + aux[3];


			}
			free(aux);
		}
	}
	return true;
}

void imprimeMatriz(double* matriz, int vertical, int horizontal, ofstream& archivoResultante) {
	for (int i = 0; i < vertical; i++) {
		for (int j = 0; j < horizontal; j++) {
			archivoResultante << fixed << setprecision(10) << matriz[i * vertical + j] << endl;
		}
	}
}

uint64_t nanos()
{
	uint64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::
		now().time_since_epoch()).count();
	return ns;
}

int main() {
	//inicializacion de las variables que se van a usar para relojes y tiempos para calculo de promedio
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

	//se inician los archivos de lectura y escritura
	ifstream lecturaA;
	ifstream lecturaB;
	ofstream archivoResultante;
	lecturaA.open("matrixA9.txt");
	lecturaB.open("matrixB9.txt");//CAMBIAR NOMBRES DE MATRICES
	archivoResultante.open("matrizResultante3.txt");

	//se pide al usuario ingresar los datos sobre las longitudes de las matrices
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

	//se revisa que se puedan multiplicar las matrices, si no se puede, se acaba el programa
	if (horizontalA != verticalB) {
		cout << "No se pueden multiplicar las matrices" << endl;
		return 0;
	}

	//se separa memoria para las matrices y se leen


	cout << "----------Obtencion de las matrices----------" << endl;

	start = clock();
	inicio = nanos();
	

	//se prepara un vector de apuntadores tipo double que se va a usar para liberar memoria despues con la funcion liberaMemoria
	vector<double*> espaciosDeMemoria;

	//se preparan las matrices en donde se va a leer de los archivos
	double* matrizA = (double*)malloc(horizontalA * verticalA * sizeof(double));

	//se separa memoria para las matrices con la siguiente logica
	if (matrizA == NULL) {
		cout << "No se pudo separar memoria para la matriz A, se va a terminar la ejecución del programa" << endl;

		return 0;
	} else {
		espaciosDeMemoria.push_back(matrizA);
	}

	double* matrizB = (double*)malloc(horizontalB * verticalB * sizeof(double));

	if (matrizB == NULL) {
		cout << "No se pudo separar memoria para la matriz B, se va a liberar la memoria previamente separada y terminara la ejecución del programa" << endl;
		liberaMemoria(espaciosDeMemoria);
		return 0;
	}
	else {
		espaciosDeMemoria.push_back(matrizB);
	}

	//se prepara la lectura para las matrices resultantes
	int verticalResultado, horizontalResultado;
	verticalResultado = verticalA;
	horizontalResultado = horizontalB;


	double* resultadoS = (double*)malloc(horizontalResultado * verticalResultado * sizeof(double));

	if (resultadoS == NULL) {
		cout << "No se pudo separar memoria para la matriz resultante de la ejecucion serial, se va a liberar la memoria previamente separada y terminara la ejecución del programa" << endl;
		liberaMemoria(espaciosDeMemoria);

		return 0;
	}
	else {
		espaciosDeMemoria.push_back(resultadoS);
	}

	double* resultadoO = (double*)malloc(horizontalResultado * verticalResultado * sizeof(double));

	if (resultadoO == NULL) {
		cout << "No se pudo separar memoria para la matriz resultante de la ejecucion OMP, se va a liberar la memoria previamente separada y terminara la ejecución del programa" << endl;
		liberaMemoria(espaciosDeMemoria);
		return 0;
	}
	else {
		espaciosDeMemoria.push_back(resultadoO);
	}

	double* resultadoI = (double*)malloc(horizontalResultado * verticalResultado * sizeof(double));

	if (resultadoI == NULL) {
		cout << "No se pudo separar memoria para la matriz resultante de la ejecucion con intrinsecas, se va a liberar la memoria previamente separada y terminara la ejecución del programa" << endl;
		liberaMemoria(espaciosDeMemoria);
		return 0;
	}
	else {
		espaciosDeMemoria.push_back(resultadoI);
	}

	double read;

	//se lee normal sin transponer
	for (int i = 0; i < verticalA; i++) {
		for (int j = 0; j < horizontalA; j++) {
			if (!lecturaA.eof()) {
				lecturaA >> read;
				matrizA[i * horizontalA + j] = read;
			} else {
				//datos incorrectos en la medida de la matriz que se introdujo
				cout << "La medida indicada de la matriz A es mayor a la longitud del archivo, se va a liberar la memoria previamente separada y terminara la ejecucion del programa" << endl;

				liberaMemoria(espaciosDeMemoria);


				return 0;
			}
		}
	}
	if (!lecturaA.eof()) {
		//datos incorrectos en la medida de la matriz que se introdujo
		cout << "La medida indicada de la matriz A es menor a la longitud del archivo, se va a liberar la memoria previamente separada y terminara la ejecucion del programa" << endl;

		liberaMemoria(espaciosDeMemoria);


		return 0;
	}

	//se lee transpuesta para hacer operaciones con memoria contigua
	for (int i = 0; i < verticalB; i++) {
		for (int j = 0; j < horizontalB; j++) {
			if (!lecturaB.eof()) {
				lecturaB >> read;
				matrizB[j * verticalB + i] = read;
			} else {
				//datos incorrectos en la medida de la matriz que se introdujo
				cout << "La medida indicada de la matriz B es mayor a la longitud del archivo, se va a liberar la memoria previamente separada y terminara la ejecucion del programa" << endl;

				liberaMemoria(espaciosDeMemoria);


				return 0;
			}
		}
	}
	if (!lecturaB.eof()) {
		//datos incorrectos en la medida de la matriz que se introdujo
		cout << "La medida indicada de la matriz B es menor a la longitud del archivo, se va a liberar la memoria previamente separada y terminara la ejecucion del programa" << endl;

		liberaMemoria(espaciosDeMemoria);


		return 0;
	}

	//se cambian los ejes de la matriz B porque se guarda transpuesta
	swap(verticalB, horizontalB);


	fin = nanos();
	end = clock();

	tiempoTranscurrido = end - start;
	tiempoEnNanosegundos = fin - inicio;

	cout << "Se configuiraron las matrices en " << tiempoTranscurrido << " ms" << endl;
	cout << "Se configuiraron las matrices en " << tiempoEnNanosegundos << " ns" << endl;


	cout << "----------Calculo de la matriz resultante----------" << endl;

	for (int i = 0; i < 5; i++) {

		cout << "-------Iteracion " << i + 1 << "-------" << endl;

		cout << "-----Serial-----" << endl;

		start = clock();
		inicio = nanos();
		//aqui se mandan a multiplicar las matrices con el metodo serial
		multiplicaDosMatrices(resultadoS, matrizA, matrizB, verticalResultado, horizontalResultado, verticalA, horizontalA, verticalB, horizontalB);

		fin = nanos();
		end = clock();

		tiempoTranscurrido = end - start;
		tiempoEnNanosegundos = fin - inicio;
		//aqui se revisa el tiempo que tarda cada uno de los procesos de multiplicacion
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
		//aqui se mandan a multiplicar las matrices con el metodo OMP
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
		//aqui se mandan a multiplicar las matrices con el metodo de intrinsecas

		bool memoriaCorrecta = multiplicaDosMatricesIntrin(resultadoI, matrizA, matrizB, verticalResultado, horizontalResultado, verticalA, horizontalA, verticalB, horizontalB);

		if (!memoriaCorrecta) {
			liberaMemoria(espaciosDeMemoria);
			return 0;
		}

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

		start = clock();
		inicio = nanos();
		//aqui se estan comparando las matrices para buscar diferencias con presicion de 10^(-10)
		for (int j = 0; j < verticalResultado; j++) {
			for (int k = 0; k < horizontalResultado; k++) {
				if (abs(resultadoS[j * verticalResultado + k] - resultadoO[j * verticalResultado + k]) > 10.0e-10 || abs(resultadoS[j * verticalResultado + k] - resultadoI[j * verticalResultado + k]) > 10.0e-10) {
					std::cout << std::fixed;
					std::cout << std::setprecision(20);
					cout << resultadoS[j * verticalResultado + k] << " " << resultadoO[j * verticalResultado + k] << " " << resultadoI[j * verticalResultado + k] << " " << endl;

					cout << (abs(resultadoS[j * verticalResultado + k] - resultadoO[j * verticalResultado + k])) << endl;
					cout << (abs(resultadoS[j * verticalResultado + k] - resultadoI[j * verticalResultado + k])) << endl;

					cout << (abs(resultadoS[j * verticalResultado + k] - resultadoO[j * verticalResultado + k]) > 10.0e-10) << endl;
					cout << (abs(resultadoS[j * verticalResultado + k] - resultadoI[j * verticalResultado + k]) > 10.0e-10) << endl;
					
			
					cout << "Error en alguna de las matrices" << endl;
					return 0;
				}
			}
		}


		cout << "Las matrices fueron calculadas correctamente" << endl;

		fin = nanos();
		end = clock();

		tiempoTranscurrido = end - start;
		tiempoEnNanosegundos = fin - inicio;


		cout << "Se comprobaron las matrices en " << tiempoTranscurrido << " ms" << endl;
		cout << "Se comprobaron las matrices en " << tiempoEnNanosegundos << " ns" << endl;

	}

	cout << "----------Guardando el resultado----------" << endl;

	start = clock();
	inicio = nanos();
	//se escribe el archivo de la matriz resultante
	imprimeMatriz(resultadoS, verticalResultado, horizontalResultado, archivoResultante);

	fin = nanos();
	end = clock();

	tiempoTranscurrido = end - start;
	tiempoEnNanosegundos = fin - inicio;

	cout << "Se guardo el archivo en " << tiempoTranscurrido << " ms" << endl;
	cout << "Se guardo el archivo en " << tiempoEnNanosegundos << " ns" << endl;

	//se cierran los archivos porque ya no se van a utilizar
	lecturaA.close();
	lecturaB.close();
	archivoResultante.close();

	//aqui se calculan los tiempos promedio de las ejecuciones
	double promedioS = (tiempoTranscurridoS1 + tiempoTranscurridoS2 + tiempoTranscurridoS3 + tiempoTranscurridoS4 + tiempoTranscurridoS5) / 5;
	double promedioO = (tiempoTranscurridoO1 + tiempoTranscurridoO2 + tiempoTranscurridoO3 + tiempoTranscurridoO4 + tiempoTranscurridoO5) / 5;
	double promedioI = (tiempoTranscurridoI1 + tiempoTranscurridoI2 + tiempoTranscurridoI3 + tiempoTranscurridoI4 + tiempoTranscurridoI5) / 5;

	double pervsSO = promedioO / promedioS;
	double pervsSI = promedioI / promedioS;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << "-----------------------------------Tabla de Comparación de Tiempo-----------------------------------" << endl;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << endl;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << "|        Corrida        |         Serial         |         Open MP        |       Intrinsecas      |" << endl;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << "|           1           |" << std::setw(23) << tiempoTranscurridoS1 << " |" << std::setw(23) << tiempoTranscurridoO1 << " |" << std::setw(23) << tiempoTranscurridoI1 << " |" << endl;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << "|           2           |" << std::setw(23) << tiempoTranscurridoS2 << " |" << std::setw(23) << tiempoTranscurridoO2 << " |" << std::setw(23) << tiempoTranscurridoI2 << " |" << endl;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << "|           3           |" << std::setw(23) << tiempoTranscurridoS3 << " |" << std::setw(23) << tiempoTranscurridoO3 << " |" << std::setw(23) << tiempoTranscurridoI3 << " |" << endl;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << "|           4           |" << std::setw(23) << tiempoTranscurridoS4 << " |" << std::setw(23) << tiempoTranscurridoO4 << " |" << std::setw(23) << tiempoTranscurridoI4 << " |" << endl;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << "|           5           |" << std::setw(23) << tiempoTranscurridoS5 << " |" << std::setw(23) << tiempoTranscurridoO5 << " |" << std::setw(23) << tiempoTranscurridoI5 << " |" << endl;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << "|        Promedio       |" << std::setw(23) << promedioS << " |" << std::setw(23) << promedioO << " |" << std::setw(23) << promedioI << " |" << endl;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << "|      % vs Serial      |           -            |" << std::setw(23) << pervsSO << " |" << std::setw(23) << pervsSI << " |" << endl;

	cout << "----------------------------------------------------------------------------------------------------" << endl;

	cout << endl;

	cout << "El metodo mas rapido fue " << ((promedioS < promedioO) ? "Serial" : ((promedioO < promedioI) ? "Open MP" : "Intrinsecas")) << endl;

	cout << endl;
	liberaMemoria(espaciosDeMemoria);
	


	return 0;
}