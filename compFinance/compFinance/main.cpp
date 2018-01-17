//
//  main.cpp
//  compFinance
//
//  Created by Samuel Tu on 1/13/18.
//  Copyright © 2018 Samuel_Tu. All rights reserved.
//  MGMTMFE405.2@gmail.com
//

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <fstream>
#include <ctime>

#define PI 3.14159265
using namespace std;

/**
 *         Usage
 *  This C++ code implements project 1 for the computational finance class.
 *  The output includes all the mean, standard deviation and time calculated.
 *  For histogram plotting, I have code to output the random variables generated
 *  to several files and plot the histogram with R. The location where the files are
 *  located is my local directory. Therefore, if you want to see the files when running
 *  code, please change the addresses in the code accordingly.
 **/

// generate uniform random variables from 0 -> m
vector<long> uniform(int seed, int size, int m, int a, int b){
    vector<long> array(size);
    array[0] = seed;
    for (int i = 1; i < size; i++) {
        array[i] = (a * array[i-1] + b) % m;
    }
    return array;
}

// output File function
void outputFile (vector<double> v, string s) {
    ofstream outputFile;
    outputFile.open(s);
    
    for (int i = 0; i < v.size(); i++) {
        //cout << v[i] << endl;
        outputFile << v[i] << "\n";
    }
    outputFile.close();
}

// calculate mean
double mean (vector<double> v) {
    double sum = 0.0;
    for (int i = 0; i < v.size(); i++) {
        sum += v[i];
    }
    return sum/v.size();
}

// calculate standard deviation
double sd(vector<double> v, double mean) {
    double sum = 0.0;
    for (int i = 0; i < v.size(); i++) {
        sum += (v[i] - mean) * (v[i] - mean);
    }
    
    return sqrt(sum/(v.size() - 1));
}

// general output function to print
void calculation (vector<double> vb) {
    double avg = mean(vb);
    double std = sd(vb, avg);
    cout << "The Empirical mean = " << avg << endl;
    cout << "The Empirical std  = " << std << endl;
}

// Question 1
vector<double> q1 (int flag, int size) {
    // set parameters
    int a = (int)pow(7.0, 5.0);
    int b = 0;
    int m = (int)pow(2, 31) - 1;
    //int size = 10000;
    int seed = 1;
    
    // generate uniform random using LGM
    vector<long> v = uniform(seed, size, m, a, b);
    
    // generate uniform from 0 -> 1 by dividing
    vector<double> vb(size);
    for (int i = 0; i < size; i++) {
        vb[i] = (double)v[i]/(double)m;
    }
    
    // since I will call this function in multiple places, where I might not want to print the following
    // I have a flag variable to control where to print.
    if (flag == 1) {
        cout << "************* Question 1 ***********" << endl;
        calculation(vb);
        
        // generating uniform random using buildin function
        vector<double> buildinUni(size);
        for (int i = 0; i < size; i++) {
            buildinUni[i] = ((double) rand() / (RAND_MAX));
        }
        double buildinAvg = mean(buildinUni);
        double buildinStd = sd(buildinUni, buildinAvg);
        cout << "Mean using buildin = " << buildinAvg << endl;
        cout << "std  using buildin = " << buildinStd << endl;
        // The empirical mean and std are very close to the ones that generates using buildin function.
    }
    
    return vb;
    // two ways of calculating mean and std. I am just testing those.
    //double sum = accumulate(vb.begin(), vb.end(), 0.0);
    //double avg2 = sum / vb.size();
    
    //double sq_sum = std::inner_product(vb.begin(), vb.end(), vb.begin(), 0.0);
    //double stdev = sqrt(sq_sum / v.size() - avg2 * avg2);
    
}

// Question 2
void q2 () {
    cout << endl << "************* Question 2 ***********" << endl;
    
    vector<double> v = q1(0, 10000); // generate 10,000 uniform using (a)
    vector<double> result;
    
    for (int i = 0; i < v.size(); i++) {
        if (v[i] >= 0 && v[i] < 0.3) {
            result.push_back(-1);
        } else if (v[i] >= 0.3 && v[i] < 0.65) {
            result.push_back(0);
        } else if (v[i] >= 0.65 && v[i] < 0.85) {
            result.push_back(1);
        } else {
            result.push_back(2);
        }
    } // for loop
    
    // address to create csv file
    string address = "/Users/samueltu/Desktop/Winter_2018/405_ComputFin/Week1/q2.csv";
    outputFile(result, address);
    calculation(result);
    
}

// getting sum to sum up bernulli
double sum (vector<double> v) {
    double sum = 0;
    for (int i = 0; i < v.size(); i++) {
        sum += v[i];
    }
    return sum;
}

// Question 3
void q3 () {
    cout << endl << "************* Question 3 ***********" << endl;
    int n = 44;
    double p = 0.64;
    int size = 1000;
    int count = 0;
    
    vector<double> result;
    vector<double> tmp = q1(0, n * size);
    
    // sum up bernulli
    for (int i = 0; i < size * n; i += 44) {
        vector<double> res;
        for (int j = i; j < i+n; j++) {
            if (tmp[j] < p) {
                res.push_back(1);
            } else {
                res.push_back(0);
            }
        } // inner for loop
        double temp_num = sum(res);
        result.push_back(temp_num);
        res.clear();
        // compute probability
        if (temp_num >= 40.0) {
            count++;
        }
    } // outter for loop
    
    // outputting files
    string address = "/Users/samueltu/Desktop/Winter_2018/405_ComputFin/Week1/q3.csv";
    outputFile(result, address);
    
    // probability that P(X >= 40) = 0
    // Theoretical value is 4.823664e-05, which is approx to 0
    double prob = (double)count/(double)size;
    cout << "count = " << count << endl;
    cout << "probability that P(X >= 40) = " << prob << endl;
    
}

// Question 4
void q4 () {
    cout << endl << "************* Question 4 ***********" << endl;
    // (a)
    int size = 10000;
    vector<double> v = q1(0, size);
    double lamda = 1.5;
    int count1 = 0;
    int count2 = 0;
    vector<double> res;
    
    // constructing vector and counting probabilities for two conditions.
    for (int i = 0; i < v.size(); i++) {
        double temp = -lamda * log(v[i]);
        res.push_back(temp);
        if (temp >= 1) {
            count1++;
        }
        if (temp >= 4) {
            count2++;
        }
    } // for loop

    string address = "/Users/samueltu/Desktop/Winter_2018/405_ComputFin/Week1/q4.csv";
    outputFile(res, address);
    
    // (b) calculate propability
    double prob1 = (double)count1/(double) size;
    double prob2 = (double)count2/(double) size;
    cout << "P(X >= 1) = " << prob1 << endl;
    cout << "P(X >= 4) = " << prob2 << endl;
    
    // (c) calculate mean and std
    calculation(res);
}

// Question 5
void q5 () {
    cout << endl << "************* Question 5 ***********" << endl;
    clock_t t1, t2;
    // (a)
    int size = 50000; // number of random variables
    
    // generate more for the pm method. choose 0.3 because the probability for failure are pi/4
    // therefore, I chose 0.3 > 0.21 (1-pi/4) to ensure there is enough number to pick.
    vector<double> v1 = q1(0, size * (1 + 0.3));
    
    // (b) Box-Muller Method
    vector<double> res_bm;
    
    // get 2 consecutive uniform random as U1 and U2
    t1 = clock();
    for (int i = 0; i < size; i += 2) {
        double z1 = sqrt(-2*log(v1[i])) * cos(2*PI*v1[i+1]);
        double z2 = sqrt(-2*log(v1[i])) * sin(2*PI*v1[i+1]);
        res_bm.push_back(z1);
        res_bm.push_back(z2);
    }
    t1 = clock() - t1;
    
    // (c) mean and sd
    cout << "Box-Muller Method:" << endl;
    calculation(res_bm);
    
    // (d) Polar-Marsaglia Method
    vector<double> res_pm;
    int size2 = (int)v1.size();

    t2 = clock();
    for (int i = 0; i < size2; i += 2) {
        // break when it reach "size" number of variables.
        if (res_pm.size() == size) break;
        double v_1 = 2 * v1[i] - 1;
        double v_2 = 2 * v1[i+1] - 1;
        double w = v_1 * v_1 + v_2 * v_2;
        
        if (w <= 1) {
            double z_1 = v_1 * sqrt((-2*log(w))/w);
            double z_2 = v_2 * sqrt((-2*log(w))/w);
            res_pm.push_back(z_1);
            res_pm.push_back(z_2);
        }
       
    }
    t2 = clock() - t2;
    
    // (e) mean and std
    cout << endl;
    cout << "Polar-Marsaglia Method:" << endl;
    calculation(res_pm);
    
    // (f)
    cout << endl;
    cout << "size = " << size << endl;
    cout << "Box-Muller time: " << ((float)t1/CLOCKS_PER_SEC) * 1000 << " ms" << endl;
    cout << "Polar-Marsaglia time: " << ((float)t2/CLOCKS_PER_SEC) * 1000 << " ms" << endl;
    // box-muller is slower because the sin and cos calculation.
}

int main(int argc, const char * argv[]) {
    /* Project 1*/
    
    q1(1, 10000);   // Question 1
    q2();           // Question 2
    q3();           // Question 3
    q4();           // Question 4
    q5();           // Question 5
    return 0;
}
