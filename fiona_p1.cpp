// Computational Finance HW1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cstring>
#include <time.h>  
#include <stdlib.h>  
#include <iomanip>  
using namespace std;


long long* LGM_generator(long long *seq, int size)
{
	long long m = pow(2, 31) - 1;
	long long a = pow(7, 5);
	long b = 0;
	for (int i = 1; i < size; i++)
	{
		seq[i] = (a*(seq[i - 1]) + b) % m;
		//	cout << seq[i] << endl;
	}
	return seq;
}

void Compare_BM_And_PM(int SourceSize, int TargetSize, int ExecutionTimes)
{
	int initial_seed = 1;
	long m = pow(2, 31) - 1;
	long long*p1 = new long long[SourceSize] {initial_seed};
	long long* uniform_pointer_normal_1;
	uniform_pointer_normal_1 = LGM_generator(p1, SourceSize);
	double *u1 = new double[SourceSize];
	for (int i = 0; i < SourceSize; i++)
	{
		u1[i] = double(*(uniform_pointer_normal_1++)) / double(m);
	}

	clock_t clockBegin_BM, clockEnd_BM, Sum_clockDuration_BM;
	clock_t clockBegin_PM, clockEnd_PM, Sum_clockDuration_PM;
	Sum_clockDuration_BM = 0;
	Sum_clockDuration_PM = 0;

	for (int k = 0; k < ExecutionTimes; k++)
	{
		double* normal_BM_test1 = new double[SourceSize];
		int num_BM_test1 = 0;
		clockBegin_BM = clock();
		for (int i = 0; i < SourceSize; i++)
		{
			normal_BM_test1[i] = sqrt(-2 * log(u1[i]))*cos(2 * 3.14 * u1[i + 1]);
			normal_BM_test1[i + 1] = sqrt(-2 * log(u1[i]))*sin(2 * 3.14 * u1[i + 1]);
			//	cout << normal_BM[i] << normal_BM[i + 1] << endl;
			num_BM_test1 += 2;
			if (num_BM_test1 == TargetSize)
			{
				clockEnd_BM = clock();
				delete[] normal_BM_test1;
				break;
			}
			i++;
		}

		Sum_clockDuration_BM = clockEnd_BM - clockBegin_BM;

		double* normal_PM_test1 = new double[TargetSize];
		int num_PM_test1 = 0;
		clockBegin_PM = clock();
		for (int i = 0; i < SourceSize; i++)
		{
			double v1 = 2 * u1[i] - 1;
			double v2 = 2 * u1[i + 1] - 1;
			double w = pow(v1, 2) + pow(v2, 2);
			if (w <= 1)
			{
				normal_PM_test1[num_PM_test1] = v1 * sqrt(-2 * log(w) / w);
				normal_PM_test1[num_PM_test1 + 1] = v2 * sqrt(-2 * log(w) / w);
			}
			//	cout << normal_BM[i] << normal_BM[i + 1] << endl;
			num_PM_test1 += 2;
			if (num_PM_test1 == TargetSize)
			{
				clockEnd_PM = clock();
				//cout << "size of PM has reached " << TargetSize << endl;
				break;
			}
			i++;
		}

		delete[] normal_PM_test1;
		Sum_clockDuration_PM += (clockEnd_PM - clockBegin_PM);
	}

//	cout << "source size: " << SourceSize << endl;
//	cout << "target size: " << TargetSize << endl;
	cout << "execution times: " << ExecutionTimes << endl;
	cout << "average execution time: " << endl;
	cout << "BM: " << (double)Sum_clockDuration_BM / ExecutionTimes << endl;
	cout << "PM: " << (double)Sum_clockDuration_PM / ExecutionTimes << endl;
	cout << endl;

	delete[] p1;
	delete[] u1;
}


int main()
{
	// Problem 1 (1)
	const int size = 10000;
	int initial_seed = 1;
	long long* seq_1 = new long long[size] { initial_seed };
	long long* uniform_pointer = new long long[size];
	uniform_pointer = LGM_generator(seq_1, size);
	//cout << uniform_pointer;
	double* uniform = new double[size];
	long m = pow(2, 31) - 1;
	double sum_LGM = 0;
	double mean_LGM = 0;
	for (int i = 0; i < size; i++)
	{
		//	cout << *uniform_pointer << endl;

		uniform[i] = double(*(uniform_pointer++)) / double(m);
		sum_LGM += uniform[i];
		//	cout << uniform[i]<<endl;
	}
	//uniform_pointer -= size;
	//cout << sum;
	mean_LGM = sum_LGM / size;
	double std_LGM = 0;
	for (int i = 0; i < size; i++)
	{
		std_LGM += pow((uniform[i] - mean_LGM), 2);
	}
	std_LGM = pow((std_LGM / size), 0.5);
	cout << "The mean of the random numbers in LGM method is:" << mean_LGM << ", and the standard deviation is :" << std_LGM << endl;

	//Problem 1	(2)
	double a = 0;
	double b = 1;
	double* target = new double[size];
	srand(1);
	double sum_built_in = 0;
	double mean_built_in = 0;
	double std_built_in = 0;
	for (int i = 0; i < size; i++)
	{
		target[i] = double(rand()) / RAND_MAX * (b - a) + a;
		//	cout << target[i] << endl;
		sum_built_in += target[i];
	}
	mean_built_in = sum_built_in / size;
	for (int i = 0; i < size; i++)
	{
		std_built_in += pow((target[i] - mean_built_in), 2);
	}
	std_built_in = pow((std_built_in / size), 0.5);
	cout << "The mean of the random numbers in built in function is:" << mean_built_in << ", and the standard deviation is :" << std_built_in << endl;

	//Problem 2(1)
	double* GDD = new double[size];
	for (int i = 0;i < size;i++)
	{
		if (uniform[i] < 0.3)
			GDD[i] = -1;
		else if (uniform[i] < 0.65)
			GDD[i] = 0;
		else if (uniform[i] < 0.85)
			GDD[i] = 1;
		else GDD[i] = 2;
		//	cout << GDD[i];
	}

	//Problem 3(1)

	const int size_bernoulli = 44000;
	long long* seq_3 = new long long[size_bernoulli] { initial_seed };
	long long* uniform_pointer_3 = new long long[size_bernoulli];
	uniform_pointer_3 = LGM_generator(seq_3, size_bernoulli);
	double* uniform_size3 = new double[size_bernoulli];
	for (int i = 0; i < size_bernoulli; i++)
	{
		uniform_size3[i] = double(*(uniform_pointer_3++)) / double(m);
	}
	int bernoulli[size_bernoulli];
	double p = 0.64;
	for (int i = 0; i < size_bernoulli; i++)
	{
		if (uniform_size3[i] < p)
			bernoulli[i] = 1;
		else bernoulli[i] = 0;
	}

	int n = 44;
	const int size_binomial = 1000;
	int* binomial = new int[size_binomial] { 0 };
	for (int i = 0; i < size_binomial; i++)
	{
		for (int j = 0;j < n;j++)
		{
			binomial[i] += bernoulli[i*n + j];
		}
		//	cout << binomial[i]<<endl;
	}

	// Problem 4(1)
	double* exponential = new double[size];
	double lamda = 1.5;
	int lar_1 = 0;
	int lar_4 = 0;
	for (int i = 0;i < size;i++)
	{
		exponential[i] = -lamda * log(1 - uniform[i]);
		if (exponential[i] >= 4)
		{
			lar_4 += 1;
			lar_1 += 1;
		}
		else if (exponential[i] >= 1)
			lar_1 += 1;
	}
	double probability_1 = double(lar_1) / double(size);
	double probability_4 = double(lar_4) / double(size);
	cout << probability_1 << endl;
	cout << probability_4 << endl;

	//Problem 5(1)
	const int size_5 = 5000;
	long long* seq_5 = new long long[size_5] { initial_seed };
	long long* uniform_pointer_5;
	uniform_pointer_5 = LGM_generator(seq_5, size_5);
	double* uniform_size5 = new double[size_5];
	for (int i = 0; i < size_5; i++)
	{
		uniform_size5[i] = double(*(uniform_pointer_5++)) / double(m);
	}
	//Problem 5(2)
	double* normal_BM = new double[size_5];
	for (int i = 0;i < size_5;i++)
	{
		normal_BM[i] = sqrt(-2 * log(uniform_size5[i]))*cos(2 * 3.14 * uniform_size5[i + 1]);
		normal_BM[i + 1] = sqrt(-2 * log(uniform_size5[i]))*sin(2 * 3.14 * uniform_size5[i + 1]);
		//	cout << normal_BM[i] << normal_BM[i + 1] << endl;
		i++;
	}
	//Problem 5(3)
	double sum_normal_BM = 0;
	double mean_normal_BM = 0;
	double std_normal_BM = 0;
	for (int i = 0;i < size_5;i++)
	{
		sum_normal_BM += normal_BM[i];
	}
	mean_normal_BM = sum_normal_BM / size_5;
	for (int i = 0;i < size_5;i++)
	{
		std_normal_BM += pow((normal_BM[i] - mean_normal_BM), 2);
	}
	std_normal_BM = pow((std_normal_BM / size_5), 0.5);

	//Problem 5(4)
	vector<double> normal_PM;
	for (int i = 0;i < size_5;i++)
	{
		double v1 = 2 * uniform_size5[i] - 1;
		double v2 = 2 * uniform_size5[i + 1] - 1;
		double w = pow(v1, 2) + pow(v2, 2);
		if (w <= 1)
		{
			normal_PM.push_back(v1*sqrt(-2 * log(w) / w));
			normal_PM.push_back(v2*sqrt(-2 * log(w) / w));
		}
		i++;
	}

	//Problem 5(5)
	//	vector<double>::iterator it; //define an iteration 
	double sum_normal_PM = 0;
	double mean_normal_PM = 0;
	double std_normal_PM = 0;
	//	for (it = normal_PM.begin();it != normal_PM.end();it++)  //use iteration to traverse
	for (int i = 0; i<normal_PM.size();i++)
	{
		//	cout << *it << " ";
		sum_normal_PM += normal_PM[i];
	}
	mean_normal_PM = sum_normal_PM / normal_PM.size();
	for (int i = 0;i<normal_PM.size();i++)
	{
		std_normal_PM += pow((normal_PM[i] - mean_normal_PM), 2);
	}
	std_normal_PM = pow((std_normal_PM / normal_PM.size()), 0.5);
	cout << mean_normal_PM << "    " << std_normal_PM << endl;

	//Problem5(6)
	Compare_BM_And_PM(10000, 5000, 200);
	Compare_BM_And_PM(20000, 10000, 200);
	Compare_BM_And_PM(30000, 20000, 200);

	//release all the pointers
	//1-1
	delete[] uniform;
	delete[] seq_1;
	//1-2
	delete[] target;
	//2-1
	delete[] GDD;
	//3-1
	delete[] seq_3;
	delete[] uniform_size3;
	delete[] binomial;
	//4-1
	delete[] exponential;
	//5-1
	delete[] seq_5;
	delete[] uniform_size5;
	//5-2
	delete[] normal_BM;
	//5-6
	/*delete[] p1;
	delete[] u1;
	delete[] normal_BM_test1;*/

	system("pause");
	return 0;
}