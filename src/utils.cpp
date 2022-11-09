#include "utils.h"
#include <math.h>
#include <cmath>
#include <limits>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>
#include <stdlib.h>
#include <bits/stdc++.h>
using namespace std;

vector<vector<double> > utils::loadMatrix(string file_name, unsigned int row, unsigned int col, int row_num, int col_num){

    double data;
    vector<vector<double> > output;
    vector<double> temp1;
    ifstream input(file_name.c_str(), ifstream::in);
    if (input.is_open()){   
        for(unsigned int i = 0; i<row+row_num; i++){
            if(i<row_num)
            {
                    for(unsigned int j=col_num; j<col+col_num; j++){
                        input >> data;
                    }
            }
            else
            {
                    for(unsigned int j=col_num; j<col+col_num; j++){
                        input >> data;
                        if(i>=row_num)
                        temp1.push_back(data);
                    }
                    output.push_back(temp1);
                    temp1.clear();
            
            }
            
        }
        input.close();
    }
    else cout << "Unable to open file";

    return output;
}

vector<double> utils::loadVector(string file_name, unsigned int l){

    double data;
    vector<double> output;
    ifstream input(file_name.c_str(), ifstream::in);
    if (input.is_open()){
        for(unsigned int i = 0; i<l; i++){
                input >> data;
                output.push_back(data);
        }
        input.close();
    }
    else cout << "Unable to open file1";

    return output;
}

vector<double> utils::loadVectorFromDataFile(string file_name){

        
        ifstream inputFile(file_name.c_str(), ifstream::in);
        vector<double> output;
        // test file open   
        if (inputFile) {        
                    double value;

                    // read the elements in the file into a vector  
                    while ( inputFile >> value ) {
                        output.push_back(value);
                          }
                       }
        else cout << "Unable to open file1";
        
    return output;
}

vector<vector<double>>  utils::readMatrixFromDataFile(string filename)
{

    std::vector<std::vector<double> > data;

    std::ifstream file(filename);
    std::string line;
    //Read one line at a time into the variable line:
    while(std::getline(file, line))
    {
    	std::vector<double> lineData;
    	std::stringstream lineStream(line);
    	double value;
        while(lineStream >> value)
        {
           lineData.push_back(value);
        }                                                                                                                             data.push_back(lineData);
                                                                                                                                   }
   return data;
}
vector<vector<double>> utils::matrixReshape(vector<vector<double> > &nums, int r, int c) {
        int x = nums.size();
        int y = nums[0].size();
        if (x * y != r * c) {
            return nums;
        }
        queue<double> data;
        for (int j=0;j<y;j++) {
            for (int i=0;i<x;i++) {
                data.push(nums[i][j]);
            }
        }
        vector<vector<double>> result(r, vector<double>(c));
        for (int col = 0; col < c; ++ col)
        {
        for (int row = 0; row < r; ++ row) 
             {                
                result[row][col] = data.front();
                data.pop();
            }
        }        
        return result;
    }

vector<double> rowSumsMatrix(vector< vector<double> > data)
{
          vector<double> rowmeans( data.size() );
          auto Sum = [](std::vector<double> const& d) { return std::accumulate(d.begin(), d.end(), 0.0); };
          std::transform(data.begin(), data.end(), rowmeans.begin(), Sum);
          return rowmeans;

}
double utils::computePctile(vector<double> g, double p){
    sort(g.begin(), g.end());
    int size = g.size();
    double pct = g[int(p*size)];
    return pct;
}

double meanVector(vector<double> v)
{

          double sum = accumulate(v.begin(),v.end(),decltype(v)::value_type(0));
          return sum/v.size();

}
void utils::transpose(vector<vector<double> > &b)
{
    if (b.size() == 0)
        return;

    vector<vector<double> > trans_vec(b[0].size(), vector<double>());

    for (int i = 0; i < b.size(); i++)
    {
        for (int j = 0; j < b[i].size(); j++)
        {
            trans_vec[j].push_back(b[i][j]);
        }
    }

    b = trans_vec;    // <--- reassign here
}

int utils::binarySearch(vector<double> v, double x)
{
	return lower_bound(v.begin(), v.end(), x)- v.begin();

}


double utils::BilinearInterpolation(double q11, double q12, double q21, double q22, double  x1, double  x2, double  y1, double  y2, double x, double y) 
{
    double x2x1, y2y1, x2x, y2y, yy1, xx1;
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    return 1.0 / (x2x1 * y2y1) * (
        q11 * x2x * y2y +
        q21 * xx1 * y2y +
        q12 * x2x * yy1 +
        q22 * xx1 * yy1
    );
}

vector<vector<double >> utils::readTabSepFile(string filename)
{

        std::fstream in(filename);
        std::vector<vector<double >>v;
        std::string line;
        int i = 0;
        while (std::getline(in, line))
        {
                double value;
                std::stringstream ss(line);
                v.push_back(std::vector<double>());
                while (ss >> value)
                {
                         v[i].push_back(value);
                }
                 ++i;
        }

        return v;


}


void utils::writeMatrix(vector<vector<double>> v,string filename, int N, int M)
{
        transpose(v);
	//cout<<"first transposed"<<endl;
        ofstream outputfile;
        outputfile.open(filename);
        for(int i=0;i<N;i++)
        {
             //cout<<"i is"<<i<<endl;
	     for(int j=0;j<M;j++)
             {
                 //cout<<"j is"<<j<<endl;
		 outputfile<<to_string(v[i][j])+"\t";
              }
           outputfile<<'\n';
        }
        transpose(v);
	//cout<<"second transposed"<<endl;
}
vector<double> utils::slicing(vector<double>& arr, 
                    int X, int Y) 
{ 
  
    // Starting and Ending iterators 
    auto start = arr.begin() + X; 
    auto end = arr.begin() + Y + 1; 
  
    // To store the sliced vector 
    vector<double> result(Y - X + 1); 
  
    // Copy vector using copy function() 
    copy(start, end, result.begin()); 
  
    // Return the final sliced vector 
    return result; 
} 
void utils::loadArray(string file_name, unsigned int l, double *pArray){

    double x;
    ifstream input(file_name.c_str(), ifstream::in);
    if (input.is_open()){
        for(unsigned int i = 0; i<l; i++){
            input >> x;
            pArray[i] = x;
        }
        input.close();
    }
    else cout << "Unable to open file1";

}

void utils::logVector(vector<double> x, string filename){

    ofstream logResult;
    logResult.open(filename.c_str(), ios::out);
    for(unsigned int i = 0; i < x.size(); i++){
        logResult << x[i] << endl;
    }
    logResult.close();
}


void utils::logVectorApp(vector<double> x, string filename){

    ofstream logResult;
    logResult.open(filename.c_str(), ios::app);
    for(unsigned int i = 0; i < x.size(); i++){
        logResult << x[i] << endl;
    }
    logResult.close();
}

vector<double> utils::columnSumMatrix(vector<vector<double>> v){

    vector<double> column_sum(v.size(),0);
    for(int i=0;i<v.size();i++)
    {
       double col_sum=0;
       for(int j=0;j<v[0].size();j++)
           col_sum+=v[i][j];
       column_sum[i]=col_sum;
    }
    return column_sum;
}

double utils::columnSumRow(vector<vector<double>> v, int row_num){
       double col_sum=0;
       for(int j=0;j<v[row_num].size();j++)
           col_sum+=v[row_num][j];
       return col_sum;
    
}

vector<double> utils::rowSumMatrix(vector<vector<double>> v){

    vector<double> row_sum(v[0].size(),0);
    for(int i=0;i<v[0].size();i++)
    {
       double row=0;
       for(int j=0;j<v.size();j++)
           row+=v[j][i];
       row_sum[i]=row;
    }
    return row_sum;
}

double utils::rowSumColumn(vector<vector<double>> v, int col_num){
       double row_sum=0;
       for(int j=0;j<v.size();j++)
           row_sum+=v[j][col_num];
       return row_sum;
    
}

double utils::interp_lin(vector<double> X, vector<double> Y, double x){

    int dim = X.size()-1;
    double y;

    // extreme cases (x<X(0) or x>X(end): extrapolation
    if(x <= X[0]){
        y = ( Y[1] - Y[0] ) / ( X[1] - X[0] )*( x - X[0] ) + Y[0] ;
        return y;
    }
    if(x >= X[dim]){
        y = Y[dim] + ( Y[dim] - Y[dim-1] ) / ( X[dim] - X[dim-1] ) * ( x - X[dim] );
        return y;
    }

    // otherwise
    // [ x - X(A) ] / [ X(B) - x ] = [ y - Y(A) ] / [ Y(B) - y ]
    // y = [ Y(B)*x - X(A)*Y(B) + X(B)*Y(A) - x*Y(A) ] / [ X(B) - X(A) ]
    double delta = 0.0;
    double min_d = numeric_limits<double>::max( );
    int j = -99;

    for(unsigned int i=0; i<X.size(); i++){
        if(X[i] == x){
            y = Y[i];
            return y;
        }
        delta = abs( X[i] - x ) ;
        if(delta < min_d){
            min_d = delta ;
            j = i;
        }
    }
    int k;
    if(X[j] < x){
        k = j;
    }else{
        k = j-1;
    }

    double a = (Y[k+1] - Y[k]) / (X[k+1] - X[k]) ;
    double b = Y[k] - a*X[k];
    y = a*x + b;

    return y;

}


double utils::gallonToCubicFeet(double x){

    double conv = 0.13368 ; // 1 gallon = 0.13368 cf
    return x*conv;
}

double utils::inchesToFeet(double x){

    double conv = 0.08333 ; // 1 inch = 0.08333 ft
    return x*conv;
}


double utils::cubicFeetToCubicMeters(double x){

    double conv = 0.0283 ; // 1 cf = 0.0283 m3
    return x*conv;
}


double utils::feetToMeters(double x){

    double conv = 0.3048 ; // 1 ft = 0.3048 m
    return x*conv;
}

double utils::acreToSquaredFeet(double x){
    double conv = 43560 ; // 1 acre = 43560 feet2
    return x*conv;
}

double utils::acreFeetToCubicFeet(double x){
    double conv = 43560 ; // 1 acre-feet = 43560 feet3
    return x*conv;
}

double utils::cubicFeetToAcreFeet(double x){
    double conv = 43560 ; // 1 acre = 43560 feet2
    return x/conv;
}



double utils::computeSum(vector<double> g){


    double z = 0.0;
    for(unsigned int i=0; i<g.size(); i++){
        z = z + g[i];
    }
    //double z = std::accumulate( g.begin(), g.end(), 0.0 ) ;
    return z;
}

double utils::computeMax(vector<double> g){
    double m = -1*numeric_limits<double>::max( );
    for(unsigned int i=0; i<g.size(); i++){
        if(g[i]>m){
            m = g[i];
        }
    }
    return m;
}

double utils::computeMin(vector<double> g){
    double m = numeric_limits<double>::max( );
    for(unsigned int i=0; i<g.size(); i++){
        if(g[i]<m){
            m = g[i];
        }
    }
    return m;
}

double utils::computeMean(vector<double> g){
    double z = computeSum(g)/g.size();
    return z;
}

double utils::computeVariance(vector<double> g){
    double v = 0.0;
    double M = computeMean(g);
    for(unsigned int i=0; i<g.size(); i++){
        v += ( g[i]-M )*( g[i]-M );
    }
    return v/g.size();
}


vector<double> utils::normalizeVector(vector<double> X, vector<double> m, vector<double> M){

    vector<double> Y;
    double z;
    for(unsigned int i=0; i<X.size(); i++){
        z = ( X[i] - m[i] ) / ( M[i] - m[i] );
        Y.push_back( z );
    }
    return Y;
}

vector<double> utils::deNormalizeVector(vector<double> X, vector<double> m, vector<double> M){

    vector<double> Y;
    double z;
    for(unsigned int i=0; i<X.size(); i++){
        z = X[i]*( M[i] - m[i] ) + m[i] ;
        Y.push_back( z );
    }
    return Y;

}


vector<double> utils::standardizeVector(vector<double> X, vector<double> m, vector<double> s){

    vector<double> Y;
    double z;
    for(unsigned int i=0; i<X.size(); i++){
        z = ( X[i] - m[i] ) / ( s[i] );
        Y.push_back( z );
    }
    return Y;
}

vector<double> utils::deStandardizeVector(vector<double> X, vector<double> m, vector<double> s){

    vector<double> Y;
    double z;
    for(unsigned int i=0; i<X.size(); i++){
        z = X[i]*s[i] + m[i] ;
        Y.push_back( z );
    }
    return Y;

}

double utils::generateRandomUnif(double lower_bound, double upper_bound){
    
    double n = ((double)rand() * (upper_bound - lower_bound)) / (double)RAND_MAX + lower_bound;
    return n;
    
}
