// sp.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "readMatrix.h"

#include <iostream> 
#include <fstream>
#include <vector>
#include <map>
#include <string>
using namespace std;

#define MALLOC( ptr, type, size )                        \
	ptr = (type*)malloc((size)* sizeof(type));                 \
if (ptr == 0) {												\
	fprintf(stderr, "!!!! malloc failed for: %s\n", #ptr); \
	exit(-1);                                                \
}
#define FREE( ptr ) free( ptr )

#define INP_MATRIX_RANK 5
#define POWER 7
#define FILENAME "D:\mat3b.mtx"
#define MATRIX_POWER_JOB 0
#define MATRIX_EXPONENT_JOB 1
#define EXPONENT_PRECISION 10


class CLParser
{
public:

	CLParser(int argc_, char * argv_[], bool switches_on_ = false);
	~CLParser(){}

	string get_arg(int i);
	string get_arg(string s);

private:

	int argc;
	vector<string> argv;

	bool switches_on;
	map<string, string> switch_map;
};

CLParser::CLParser(int argc_, char * argv_[], bool switches_on_)
{
	argc = argc_;
	argv.resize(argc);
	copy(argv_, argv_ + argc, argv.begin());
	switches_on = switches_on_;

	//map the switches to the actual
	//arguments if necessary
	if (switches_on)
	{
		vector<string>::iterator it1, it2;
		it1 = argv.begin();
		it2 = it1 + 1;

		while (true)
		{
			if (it1 == argv.end()) break;
			if (it2 == argv.end()) break;

			if ((*it1)[0] == '-')
				switch_map[*it1] = *(it2);

			it1++;
			it2++;
		}
	}
}

string CLParser::get_arg(int i)
{
	if (i >= 0 && i<argc)
		return argv[i];

	return "";
}

string CLParser::get_arg(string s)
{
	if (!switches_on) return "";

	if (switch_map.find(s) != switch_map.end())
		return switch_map[s];

	return "";
}

double Determinant(double *a, int n)
{
	int i, j, j1, j2;
	double det = 0;
	double *m;

	if (n < 1) {

	}
	else if (n == 1) {
		det = a[0];
	}
	else if (n == 2) {
		det = a[0] * a[3] - a[1] * a[2];
	}
	else {
		det = 0;
		for (j1 = 0; j1<n; j1++) {
			MALLOC(m, double, (n - 1)*(n - 1));
			for (i = 1; i<n; i++) {
				j2 = 0;
				for (j = 0; j<n; j++) {
					if (j == j1)
						continue;
					m[(i - 1)*(n-1)+j2] = a[i*n+j];
					j2++;
				}
			}
			det += pow(-1.0, 1.0 + j1 + 1.0) * a[j1] * Determinant(m, n - 1);
			FREE(m);
		}
	}
	return(det);
}

inline int Factorial(int x) {
	return (x == 0 ? 1 : (x == 1 ? x : x * Factorial(x - 1)));
}

int takeElement(int k, int n, int last, int *offsetBlocks, const int kInitial, int *C, int *CBlock) {
	int t,i;
	if (k == 0) { 
		for (i = 0; i < kInitial; i++) {
			C[(kInitial - i - 1) + offsetBlocks[0] * kInitial] = CBlock[i];
		}
		return 0;	
	};
	for (t = (last+1); t <= (n-k+1); t++){ 
		if (last + 1 - t < 0) { offsetBlocks[0] = offsetBlocks[0] + 1; }
		CBlock[k - 1] = t;
		takeElement(k - 1, n, t, offsetBlocks, kInitial, C, CBlock);
	}
}

int calculateElementarySymmetricPolynomials(int n, double *sigma, double *A) {
	int i, size_C, k, s, kInitial, nbrBlocks,ii,jj;
	int *C, *CBlock, *offsetBlocks;
	double det;
	double *M;
	for (i = 0; i < n; i++) {
		k = i+1;
		kInitial = k;
		nbrBlocks = static_cast<int>(Factorial(n) / (Factorial(k)*Factorial(n - k)));
		size_C = nbrBlocks * k;
		MALLOC(C, int, size_C);
		MALLOC(CBlock, int, k);
		MALLOC(offsetBlocks, int, 1);
		offsetBlocks[0] = 0;
		takeElement(k, n, 0, offsetBlocks, kInitial, C, CBlock);
		FREE(CBlock);
		FREE(offsetBlocks);

		MALLOC(M, double, k*k);
		sigma[i] = 0;
		for (s = 0; s < nbrBlocks; s++) {
			for (ii = 0; ii < k; ii++) {
				for (jj = 0; jj < k; jj++) {
					M[ii*k + jj] = A[(C[s*k + ii]-1) * n + (C[s*k + jj]-1)];
				}
			}
			det = Determinant(M, k);
			sigma[i] += det;
		}
		FREE(M);
	}
}

int calculateSymmetricPolynomialsNth(int maxg, int n, double *sigma, double *B, bool withAlternateFunction) {
	int g,i;
	double sum;
	for (g = 0; g <= maxg; g++) {
		if (g <= (n - 2)) {
			B[g] = 0.0;
		}
		else if (g == (n - 1)) {
			B[g] = 1.0;
		}
		else {
			sum = 0;
			for (i = 1; i <= n; i++) {
				if (withAlternateFunction) {
					sum += pow(-1.0, ((static_cast<double>(i)) - 1.0)) * sigma[i - 1] * B[g - i];
				}
				else {
					sum += sigma[i] * B[g - i];
				}
			}
			B[g] = sum;
		}
	}
	return 0;
}

int squareMatrixCopy(int n, double *A, double *B, int offsetA, int offsetB) {
	int nn;
	for (nn = 0; nn < n*n; nn++) { 
		B[nn + offsetB] = A[nn + offsetA];
	}
	return 0;
}

int matrixByMatrixMultiplication(int n, double *A, double *B, double *C, int offsetA, int offsetB, int offsetC) {
	int i, j, k;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			C[i*n + j + offsetC] = 0;
			for (k = 0; k < n; k++) {
				C[i*n + j + offsetC] += A[i*n + k + offsetA] * B[k*n + j + offsetB];
			}
		}
	}
	return 0;
}

int matrixPowerRoutine(int n, double *A, double *sigma, double *B, int power, double *X) {
	double multipler;
	int pw, i, j, nn, l, g;		
	for (int pw = 0; pw <= power; pw++)
	{
		if (pw == 0)
		{
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					if (i == j) {
						X[i*n + j] = 1;
					}
				}
			}
		}
		else if (pw == 1) {
			for (nn = 0; nn < n*n; nn++) {
				X[n*n + nn] = A[nn];
			}
		}
		else if (pw <= n + 1) {
			matrixByMatrixMultiplication(n, A, X, X, 0, (pw - 1)*n*n, (pw)*n*n);
		}
		else {
			for (l = 0; l < n; l++) {
				multipler = 0;
				for (g = 0; g <= l; g++) {
					multipler += pow(-1.0, ((static_cast<double>(n - l + g)) - 1.0)) * sigma[n - l + g - 1] * B[pw - 1 - g];
				}
				for (nn = 0; nn < n*n; nn++) {
					X[(pw)*n*n + nn] += X[(l)*n*n + nn] * multipler;
				}
			}
		}
	}
	return 0;
}

int matrixPower(double *A, int n, int power) {
	double *sigma, *B, *X;
	int nn;
	/* sigma array (elementary symmetric polynomials) calculation (start) */
	MALLOC(sigma, double, n);
	calculateElementarySymmetricPolynomials(n, sigma, A);
	/* sigma array (elementary symmetric polynomials) calculation (finish) */

	/* Symmetric polynomials n'th order calculation (start) */
	MALLOC(B, double, power + 1);
	calculateSymmetricPolynomialsNth(power, n, sigma, B, true);
	/* Symmetric polynomials n'th order calculation (finish) */

	/* Result matrix calculation (start) */
	MALLOC(X, double, (power + 1)*n*n);
	for (nn = 0; nn < (power + 1)*n*n; nn++) { X[nn] = 0; }
	matrixPowerRoutine(n, A, sigma, B, power, X);
	FREE(sigma);
	FREE(B);
	squareMatrixCopy(n, X, A, power*n*n, 0);
	FREE(X);
	/* Result matrix calculation (finish) */
}

double matrixTrace(double *X, int n, int offsetX) {
	int i;
	double tr;
	tr = 0;
	for (i = 0; i < n; i++) {
		tr += X[i*n + i + offsetX];
	}
	return tr;
}

int calculateTraces(double *X, int n, double *S) {
	int g;
	for (g = 0; g <= n; g++) {
		S[g] = matrixTrace(X, n, g*n*n);
	}
	return 0;
}

int calculateCoefficientsP(double *S, int n, double *P) {
	int g,i;
	double sum;
	P[0] = 0;
	for (g = 1; g <= n; g++) {
		sum = S[g];
		for (i = 1; i < g; i++) {
			sum -= P[i] * S[g - i];
		}
		P[g] = sum/g;
	}
	return 0;
}

int calculateFactorialArray(double *F, int n) {
	int i;
	for (i = 0; i <= n; i++) {
		F[i] = Factorial(i);
	}
	return 0;
}

int calculateSums(int n, int J, double *B, double *F, double *SUM) {
	int g,j;
	for (g = 0; g < n; g++) {
		SUM[g] = 0;
		for (j = n; j <= J; j++) {
			SUM[g] += B[j - 1 - g] / F[j];
		}
	}
	return 0;
}

int calculateMultipler(int n, double *P, double *SUM, double *F, double *M) {
	int i, l;
	for (i = 0; i < n; i++) {
		M[i] = 1 / F[i];
		for (l = 0; l <= i; l++) {
			M[i] += P[n - l] * SUM[i - l];
		}
	}
	return 0;
}

int matrixExponentialRoutine(int n, double *X, double *M, double *EXP) {
	int l,nn;
	for (nn = 0; nn < n*n; nn++) {
		EXP[nn] = 0;
	}
	for (l = 0; l < n; l++) {
		for (nn = 0; nn < n*n; nn++) {
			EXP[nn] += X[nn + l*n*n] * M[l];
		}
	}
	return 0;
}

int matrixExponential(double *A, int n, int J) {
	double *X, *B, *sigma, *S, *P, *F, *SUM, *M, *EXP;
	int nn;

	/* Calculate matrix powers (from 0 to n) (start) */
	MALLOC(X, double, (n + 1)*n*n);
	for (nn = 0; nn < (n + 1)*n*n; nn++) { X[nn] = 0; }
	MALLOC(B, double, 1); //just initialization for matrixPowerRoutine (not used this time)
	MALLOC(sigma, double, 1); //just initialization for matrixPowerRoutine (not used this time)
	matrixPowerRoutine(n, A, sigma, B, n, X);
	FREE(B);
	FREE(sigma);
	/* Calculate matrix powers (from 0 to n) (finish) */

	/* Calculate traces (start) */
	MALLOC(S, double, n + 1);
	calculateTraces(X, n, S);
	/* Calculate traces (finish) */

	/* Calculate P array (start) */
	MALLOC(P, double, n + 1);
	calculateCoefficientsP(S, n, P);
	FREE(S);
	/* Calculate P array (finish) */

	/* Symmetric polynomials n'th order calculation (start) */
	MALLOC(B, double, J);
	calculateSymmetricPolynomialsNth(J-1, n, P, B, false);
	/* Symmetric polynomials n'th order calculation (finish) */

	/* Calculate Factorials (start) */
	MALLOC(F, double, J+1);
	calculateFactorialArray(F, J);
	/* Calculate Factorials (finish) */

	/* Calculate Sums (start) */
	MALLOC(SUM, double, n);
	calculateSums(n,J,B,F,SUM);
	FREE(B);
	/* Calculate Sums (finish) */

	/* Calculate Multipler (start) */
	MALLOC(M, double, n);
	calculateMultipler(n,P,SUM,F,M);
	FREE(P);
	FREE(F);
	FREE(SUM);
	/* Calculate Multipler (finish) */

	/* Calculate Matrix Exponential (start) */
	MALLOC(EXP, double, n*n);
	matrixExponentialRoutine(n, X, M, EXP);
	squareMatrixCopy(n, EXP, A, 0, 0);
	FREE(M);
	FREE(X);
	FREE(EXP);
	/* Calculate Matrix Exponential (finish) */
}

int main(int argc, char * argv[])
{
	int n = INP_MATRIX_RANK, power = POWER, job = MATRIX_EXPONENT_JOB, exp_precision = EXPONENT_PRECISION; //default values
	
	/* Program arguments parsing (start) */
	CLParser cmd_line(argc, argv, true);
	std::string temp;
	std::string filename;
	std::string out_filename;
	temp = cmd_line.get_arg("-n");
	n = atoi(temp.c_str());
	cout << "Input matrix rank is: " << n << "\n";
	temp = cmd_line.get_arg("-p");
	power = atoi(temp.c_str());
	cout << "Specified matrix power is: " << power << "\n";
	temp = cmd_line.get_arg("-e");
	exp_precision = atoi(temp.c_str());
	cout << "Specified exponential precision is: " << exp_precision << "\n";
	temp = cmd_line.get_arg("-j");
	job = atoi(temp.c_str());
	if (job == 0) {
		cout << "Specified job is: matrix power \n";
	}
	else if (job == 1) {
		cout << "Specified job is: matrix exponential \n";
	}
	filename = cmd_line.get_arg("-f");
	cout << "Specified input file name is: " << filename.c_str() << "\n";
	out_filename = cmd_line.get_arg("-o");
	cout << "Specified output file name is: " << out_filename.c_str() << "\n";
	/* Program arguments parsing (finish) */

	/* Initial matrix reading (start) */
	double *A;
	MALLOC(A, double, n*n);
	if (read_matrix(n, n, A, filename.c_str()) != 0) {
		printf("input matrix contains errors");
		exit(-1);
	};
	/* Initial matrix reading (finish) */

	/* Execute program (depending on job) (start)*/
	if (job == MATRIX_POWER_JOB) {
		matrixPower(A, n, power);
	}
	else if (job == MATRIX_EXPONENT_JOB) {
		matrixExponential(A, n, exp_precision);
	}
	/* Execute program (depending on job) (finish)*/

	/* Save results to file (start) */
	ofstream file;
	file.open(out_filename.c_str());
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			file << i << " " << j << " " << A[i*n+j] << "\n";
		}
	}
	file.close();
	/* Save results to file (finish) */

	return 0;
}

