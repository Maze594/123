#include <math.h>
#include <complex>
#include <iostream>
#include <fstream>
using namespace std;

#define PI 3.141592653589793

const int N = 625;
// 
const int N2 = 25;
const int N1 = 25;
// 
const int N4 = 5;
const int N3 = 5;

// подсчет числа операций
int m = 0;
int a = 0;
// входной сигнал
complex <double> X[N];
complex <double>* DFT(complex<double> B[], int count);
void BinReverse(complex<double> A[], int N);
void CooleyTukey5_4(complex<double> B[], bool INV);
void CooleyTukey32_20(complex<double> B[], bool INV);

int main()
{
	setlocale(LC_ALL, "Russian");
	string fileName = "file1.txt";
	fstream file(fileName);
	file.open(fileName, std::ios::app);

	int fd = 1200;
	complex<double> *X_2 = new complex<double>[N];

	// входной сигнал
	for (int i = 0; i < N; i++)
		X[i] = sin(2 * PI * 120 / fd * i) + sin(2 * PI * 240 / fd * i);
	cout << "--Входной сигнал Х: сохранено в file.txt" << endl;
	file << "--Входной сигнал Х:" << endl;
	for (int i = 0; i < N; i++)
		file << X[i] << "  ";
	cout << endl;
	file << endl << endl;

	// ДПФ по определению 
	X_2 = DFT(X, N);
	cout << "--ДПФ по определению: сохранено в file.txt" << endl;
	file << "--ДПФ по определению:" << endl;
	for (int i = 0; i < N; i++)
		file << X_2[i].real() << "+" << X_2[i].imag() << ",";
	file << endl << endl;
	cout << "\tКоличество сложений:  " << a - N << endl;
	cout << "\tКоличество умножений:  " << m << endl << endl;

	// Кули-Тьюки
	cout << "--Кули-Тьюки: сохранено в file.txt" << endl;
	file << "--Кули-Тьюки:" << endl;
	m = 0;
	a = 0;
	CooleyTukey32_20(X, false);
	for (int i = 0; i < N; i++)
		file << X[i].real() << "+" << X[i].imag() << ",";
	file << endl << endl;
	cout << "\tКоличество сложений:  " << a - 4 * N << endl;
	cout << "\tКоличество умножений:  " << m << endl << endl;

	// Обратное БПФ
	cout << "--Обратное БПФ: сохранено в file.txt" << endl;
	file << "--Обратное БПФ:" << endl;
	m = 0;
	a = 0;
	CooleyTukey32_20(X, true);
	for (int i = 0; i < N; i++)
		file << X[i] << "  ";
	file << endl;
	cout << "\tКоличество сложений:  " << a - 4 * N << endl;
	cout << "\tКоличество умножений:  " << m << endl << endl;

	file.close();
	system("pause");
	return 0;
}

// ДПФ по определению
complex <double>* DFT(complex<double> B[], int count)
{
	complex<double> *A = new complex<double>[count];
	complex<double> W;

	for (int k = 0; k < count; k++)
	{
		A[k] = 0;
		for (int n = 0; n < count; n++)
		{
			W.real(cos(2 * PI / count * k*n));
			W.imag((-1)*sin(2 * PI / count * k*n));
			A[k] = A[k] + B[n] * W;
			a++;
			m++;
		}
	}
	return A;
}

// двоичное реверсирование
void BinReverse(complex<double> A[], int N)
{
	complex <double> T;
	int k, j, NV2;
	NV2 = (int)round((double)N / 2);
	j = 1;

	for (int i = 1; i < N; i++)
	{
		if (i < j)
		{
			T = A[j - 1];
			A[j - 1] = A[i - 1];
			A[i - 1] = T;
		}
		k = NV2;
		while (k < j)
		{
			j -= k;
			k = (int)round((double)k / 2);
		}
		j += k;
	}
}

// Кули-Тьюки для 25
void CooleyTukey5_4(complex<double> B[], bool INV)
{
	// для преобразования в двумерный
	complex<double> Y[N4][N3], W;
	// временное хранение данных для N4
	complex<double> *Tmp_N4 = new complex<double>[N4];
	// временное хранение данных для N3
	complex<double> *Tmp_N3 = new complex<double>[N3];

	// обратное БПФ
	if (INV)
	{
		for (int i = 0; i < N; i++)
			X[i] = conj(X[i]);
	}

	// преобразование одномерный массив в двумерный
	for (int i = 0; i < N3; i++)
		for (int j = 0; j < N4; j++)
			Y[j][i] = X[i + j * N3];

	// 5 точечное ДПФ по определению
	for (int i = 0; i < N3; i++)
	{
		for (int j = 0; j < N4; j++)
			Tmp_N4[j] = Y[j][i];
		Tmp_N4 = DFT(Tmp_N4, N4);
		for (int j = 0; j < N4; j++)
			Y[j][i] = Tmp_N4[j];
	}

	// умножение на поворачивающий множитель
	for (int i = 0; i < N4; i++)
	{
		for (int j = 0; j < N3; j++)
		{
			W.real(cos(2 * PI / N * i*j));
			W.imag((-1) * sin(2 * PI / N * i*j));
			Y[i][j] = Y[i][j] * W;
		}
	}
	m += N4 * N3;

	// 5 точечное ДПФ по определению
	for (int i = 0; i < N4; i++)
	{
		for (int j = 0; j < N3; j++)
			Tmp_N3[j] = Y[j][i];
		Tmp_N3 = DFT(Tmp_N3, N3);
		for (int j = 0; j < N3; j++)
			Y[j][i] = Tmp_N3[j];
	}

	// преобразовать двумерный массив в одномерный
	for (int i = 0; i < N3; i++)
		for (int j = 0; j < N4; j++)
			X[i*N4 + j] = Y[j][i];

	if (INV)
	{
		for (int i = 0; i < N; i++)
			X[i] /= N;
	}
}

// Кули-Тьюки для 625
void CooleyTukey32_20(complex<double> B[], bool INV)
{
	// для преобразования в двумерный
	complex<double> Y[N2][N1], W;
	// временное хранение данных для N2
	complex<double> *Tmp_N2 = new complex<double>[N2];
	// временное хранение данных для N1
	complex<double> *Tmp_N1 = new complex<double>[N1];

	// обратное БПФ
	if (INV)
	{
		for (int i = 0; i < N; i++)
			X[i] = conj(X[i]);
	}

	// преобразование одномерный массив в двумерный
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			Y[j][i] = X[i + j * N1];

	m += N1 * N2;


	for (int i = 0; i < N1; i++)
	{
		for (int j = 0; j < N2; j++)
			Tmp_N2[j] = Y[j][i];
		CooleyTukey5_4(Tmp_N2, INV);
		for (int j = 0; j < N2; j++)
			Y[j][i] = Tmp_N2[j];
	}

	for (int i = 0; i < N2; i++)
	{
		for (int j = 0; j < N1; j++)
			Tmp_N1[j] = Y[j][i];
		CooleyTukey5_4(Tmp_N1, INV);
		for (int j = 0; j < N1; j++)
			Y[j][i] = Tmp_N1[j];
	}

	// преобразовать двумерный массив в одномерный
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			X[i*N2 + j] = Y[j][i];

	// обратное БПФ
	if (INV)
	{
		for (int i = 0; i < N; i++)
			X[i] = conj(X[i]);
	}
}