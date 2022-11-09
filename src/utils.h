#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <bits/stdc++.h>
namespace std{

struct myFile{
    string filename;
    int row;
    int col;
};


class utils
{
public:

    /**
     * functions to load files
     */
    static vector<vector<double> > loadMatrix(string file_name, unsigned int row, unsigned int col, int row_num, int col_num);
    static vector<double> loadVector(string file_name, unsigned int l);
    static void loadArray(string file_name, unsigned int l, double* pArray);
    static vector<vector<double>>  readMatrixFromDataFile(string file_name);
    static vector<vector<double>> matrixReshape(vector<vector<double> > &nums, int r, int c);
    static vector<vector<double >> readTabSepFile(string filename);
    static vector<double> loadVectorFromDataFile(string file_name);
    static vector<double> rowSumsMatrix(vector< vector<double> > data);
    static double meanVector(vector<double> v);
    static double computePctile(vector<double> g, double p);
    static vector<double> columnSumMatrix(vector<vector<double>> v);
    static double columnSumRow(vector<vector<double>> v, int row_num);
    static vector<double> rowSumMatrix(vector<vector<double>> v);
    static double rowSumColumn(vector<vector<double>> v, int col_num);
    /**
     * functions to transport matrices
     */
     static void transpose(vector<vector<double> > &b);
     static vector<double> slicing(vector<double>& arr, int X, int Y);
     static void writeMatrix(vector<vector<double>> v,string filename, int N, int M);
    /**
     * function to print a vetor to a file (App append the file)
     */
    static void logVector(vector<double> x, string file_name);
    static void logVectorApp(vector<double> x, string file_name);

    /**
      * linear interpolation
      */
    static double interp_lin(vector<double> X, vector<double> Y, double x);
    static double BilinearInterpolation(double q11, double q12, double q21, double q22, double  x1, double  x2, double  y1, double  y2, double x, double y);
    static int binarySearch(vector<double> v, double x);

    /**
      * conversion from gallon to cubic feet
      */
    static double gallonToCubicFeet(double x);

    /**
      * conversion from inches to feet
      */
    static double inchesToFeet(double x);

    /**
      * conversion from cubic feet to cubic meters
      */
    static double cubicFeetToCubicMeters(double x);

    /**
      * conversion from feet to meters
      */
    static double feetToMeters(double x);

    /**
      * conversion from acre to squared feet
      */
    static double acreToSquaredFeet(double x);

    /**
      * conversion from acre-feet to cubic feet
      */
    static double acreFeetToCubicFeet(double x);

    /**
      * conversion from cubic feet to acre-feet
      */
    static double cubicFeetToAcreFeet(double x);

    /**
      * Basic operations on vector
      */
    static double computeSum(vector<double> g);
    static double computeMax(vector<double> g);
    static double computeMin(vector<double> g);
    static double computeMean(vector<double> g);
    static double computeVariance(vector<double> g);

    /**
      * Normalization and De-normalization
      */
    static vector<double> normalizeVector( vector<double> X, vector<double> m, vector<double> M );
    static vector<double> deNormalizeVector( vector<double> X, vector<double> m, vector<double> M );


    /**
      * Standardization and De-standardization
      */
    static vector<double> standardizeVector( vector<double> X, vector<double> m, vector<double> s );
    static vector<double> deStandardizeVector( vector<double> X, vector<double> m, vector<double> s );

    /**
     * Random number from uniform distribution
     */
    static double generateRandomUnif(double lower_bound, double upper_bound);

};
}

#endif // UTILS_H
