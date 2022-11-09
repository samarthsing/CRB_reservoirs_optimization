/*
 * param_function.cpp
 *
 *  Created on: 31/oct/2013
 *      Author: MatteoG
 */


#include "param_function.h"

using namespace std;

param_function::param_function()
{
}

param_function::~param_function()
{
}

unsigned int param_function::getInputNumber(){
    return M;
}

unsigned int param_function::getOutputNumber(){
    return K;
}

vector<double> param_function::get_StdOutput(vector<double> pInput){

    vector<double> x, y, z;

    x = utils::standardizeVector( pInput, input_mean, input_std );
    z = get_output( x ) ;
    y = utils::deStandardizeVector( z, output_mean, output_std );
    return y ;
}

vector<double> param_function::get_NormOutput(vector<double> pInput){

    vector<double> x, y, z;
    /*
    cout<<"pInput.size() "<<pInput.size()<<endl;
    for(auto zz: pInput)
        cout<<zz<<" "<<endl;
    
    cout<<"input_min.size() "<<input_min.size()<<endl;
    for(auto zz: input_min)
        cout<<zz<<" "<<endl;

    cout<<"input_max.size() "<<input_max.size()<<endl;
    for(auto zz: input_max)
        cout<<zz<<" "<<endl;
    cout<<"output_min.size() "<<output_min.size()<<endl;
    for(auto zz: output_min)
        cout<<zz<<" "<<endl;
    cout<<"output_max.size() "<<output_max.size()<<endl;
    for(auto zz: output_max)
        cout<<zz<<" "<<endl;

    */

    

    

    x = utils::normalizeVector( pInput, input_min, input_max );
    /*
    cout<<"x.size() "<<x.size()<<endl;
    for(auto zz: x)
        cout<<zz<<" "<<endl;*/
    z = get_output( x ) ;
    /*
    cout<<"z.size() "<<z.size()<<endl;
    for(auto zz: z)
        cout<<zz<<" "<<endl;*/
    
    y = utils::deNormalizeVector( z, output_min, output_max );
    /*
    cout<<"y.size() "<<y.size()<<endl;
    for(auto zz: y)
        cout<<zz<<" "<<endl;*/
    return y ;
}

vector<double> param_function::get_NormOutputMaster(vector<double> pInput,unsigned int rank_no,unsigned int master_no){

    vector<double> x, y, z;
    
    cout<<"pInput.size() "<<pInput.size()<<endl;
    for(auto zz: pInput)
        cout<<zz<<" "<<rank_no<<" "<<master_no<<endl;
    
    cout<<"input_min.size() "<<input_min.size()<<endl;
    for(auto zz: input_min)
        cout<<zz<<" "<<rank_no<<" "<<master_no<<endl;

    cout<<"input_max.size() "<<input_max.size()<<endl;
    for(auto zz: input_max)
        cout<<zz<<" "<<rank_no<<" "<<master_no<<endl;
    cout<<"output_min.size() "<<output_min.size()<<endl;
    for(auto zz: output_min)
        cout<<zz<<" "<<rank_no<<" "<<master_no<<endl;
    cout<<"output_max.size() "<<output_max.size()<<endl;
    for(auto zz: output_max)
        cout<<zz<<" "<<rank_no<<" "<<master_no<<endl;

    

    

    

    x = utils::normalizeVector( pInput, input_min, input_max );
    cout<<"x.size() "<<x.size()<<endl;
    for(auto zz: x)
        cout<<zz<<" "<<rank_no<<" "<<master_no<<endl;
    z = get_output( x ) ;
    cout<<"z.size() "<<z.size()<<endl;
    for(auto zz: z)
        cout<<zz<<" "<<rank_no<<" "<<master_no<<endl;
    
    y = utils::deNormalizeVector( z, output_min, output_max );
    cout<<"y.size() "<<y.size()<<endl;
    for(auto zz: y)
        cout<<zz<<" "<<rank_no<<" "<<master_no<<endl;
    return y ;
}

void param_function::setMaxInput(vector<double> pV){

    for(unsigned int i=0; i<pV.size(); i++){
        input_max.push_back( pV[i] );
    }
}

void param_function::setMaxOutput(vector<double> pV){

    for(unsigned int i=0; i<pV.size(); i++){
        output_max.push_back( pV[i] );
    }
}

void param_function::setMinInput(vector<double> pV){

    for(unsigned int i=0; i<pV.size(); i++){
        input_min.push_back( pV[i] );
    }
}

void param_function::setMinOutput(vector<double> pV){

    for(unsigned int i=0; i<pV.size(); i++){
        output_min.push_back( pV[i] );
    }
}
/*
vector<double> param_function::getMaxInput(){
    return input_max;
}

vector<double> param_function::getMaxOutput(){
    return output_max;
}

vector<double> param_function::getMinInput(){
    return input_min;
}

vector<double> param_function::getMinOutput(){
    return output_min;
}
*/

void param_function::setMeanInput(vector<double> pV){

    for(unsigned int i=0; i<pV.size(); i++){
        input_mean.push_back( pV[i] );
    }
}

void param_function::setMeanOutput(vector<double> pV){

    for(unsigned int i=0; i<pV.size(); i++){
        output_mean.push_back( pV[i] );
    }
}

void param_function::setStdInput(vector<double> pV){

    for(unsigned int i=0; i<pV.size(); i++){
        input_std.push_back( pV[i] );
    }
}

void param_function::setStdOutput(vector<double> pV){

    for(unsigned int i=0; i<pV.size(); i++){
        output_std.push_back( pV[i] );
    }
}

/*
vector<double> param_function::getMeanInput(){
    return input_mean;
}

vector<double> param_function::getMeanOutput(){
    return output_mean;
}

vector<double> param_function::getStdInput(){
    return input_std;
}

vector<double> param_function::getStdOutput(){
    return output_std;
}
*/




