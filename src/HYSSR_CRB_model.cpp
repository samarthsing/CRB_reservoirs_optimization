#include "utils.h"
#include "utils.cpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <limits>
#include <map>
#include <fstream>
#include <numeric>
#include <stdlib.h>
#include <array>
#include <sstream>
#include <bits/stdc++.h>
#include <cmath>
#include <functional>
#include <numeric>
#include "param_function.h"
#include "rbf.h"
#include "param_function.cpp"
#include "rbf.cpp"
#include "ncRBF.cpp"
#include "ncRBF.h"
using namespace std;

class MyClass {
public:

const int k=1000;
const double EPSILON=0.00000001;
const double pi=M_PI;
const int M=47;
const int M_comp=55;
const int M_H=10;
const int M_HR=15;
int sim_years=100;
int sim_days=sim_years*365;
int no_years=sim_years+2;
int bpa_rev_sim_years=sim_years-1;
int bpa_rev_sim_days=sim_days-365;
int Nobj=1;
bool write_output=true;
bool CompareDoubles2 (double A, double B)
{
   double diff = A - B;
   return (diff < EPSILON) && (-diff < EPSILON);
}

double cfs_kaf_conversion_factor = 0.001987;
double powerhouse_conversion_factor= 14.27641;
double density_water=1000;
double gravity=9.81;
double head_metre_conversion=0.3048;
double kaf_ksfd_conversion=0.504165;
//int start_row=173010;
int start_row=243;
int start_col=0;
vector<vector<double> > modelled_discharge=utils::loadMatrix("../data/canadian_dams/discharge_48.txt",sim_days+1,M,365-1,start_col);//water year
vector<vector<double> > modelled_storage=utils::loadMatrix("../data/canadian_dams/storage_48.txt",sim_days+1,M,365-1,start_col);//water year
vector<vector<double> > modelled_generation=utils::loadMatrix("../data/canadian_dams/generation_48.txt",sim_days+1,M,365-1,start_col);//water year
vector<vector<double> > modelled_powerflow=utils::loadMatrix("../data/canadian_dams/powerflow_48.txt",sim_days+1,M,365-1,start_col);//water year
vector<vector<double> > inflow=utils::loadMatrix("../data/canadian_dams/inflow_candians.txt",sim_days+1,M_comp,start_row+365-1,start_col);//jan 1
vector<vector<double> > maxhydro_matrix=utils::readMatrixFromDataFile("../data/HYSSR/max_hydropower.txt");
vector<vector<double> > maxstorage_matrix=utils::readMatrixFromDataFile("../data/HYSSR/max_storage.txt");
vector<double> power_max;
vector<double> s_max; 
vector<vector<double> > StorEq=utils::readMatrixFromDataFile("../data/HYSSR/storage_elevation.txt");
vector<vector<double> > TailEq=utils::readMatrixFromDataFile("../data/HYSSR/tailwater_elevation.txt");
vector<vector<double> > power_constants=utils::readMatrixFromDataFile("../data/HYSSR/generations_constant.txt");
vector<double> hydraulic_head;
vector<double> powerhouse_efficiency;
vector<double> spokane_flow;
//data structures to model HYSSR
vector<vector<double> > discharge=modelled_discharge;
//vector<vector<double> > storage(M,vector<double> (sim_days+1,0));
vector<vector<double> > storage=modelled_storage;
vector<vector<double> > total_inflow=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > generation=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > spill=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > powerflow=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > tailwater=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > stage_hyssr=vector<vector<double> >(M_H,vector<double> (sim_days+1,0));
vector<vector<double> > tailwater_hyssr=vector<vector<double> >(M_H,vector<double> (sim_days+1,0));
vector<vector<double> > heads=vector<vector<double> >(M_H,vector<double> (sim_days+1,0));
vector<vector<double> > rbfs_releases=vector<vector<double> >(3,vector<double> (sim_days,0));
vector<vector<double> > rbfs_inputs=vector<vector<double> >(8,vector<double> (sim_days,0));

//Data structures to model GCL(storage dam)
vector<vector<double> > calender=utils::readMatrixFromDataFile("../data/HYSSR/GCL/calender.txt");
vector<vector<double> > GCL_CRC=utils::readMatrixFromDataFile("../data/HYSSR/GCL/GCL_CRC_hm3.txt");
vector<vector<double> > GCL_ARC=utils::readMatrixFromDataFile("../data/HYSSR/GCL/GCL_ARC_hm3.txt");
vector<vector<double> > GCL_fc=utils::readMatrixFromDataFile("../data/HYSSR/GCL/GCL_fc.txt");
vector<vector<double> > c=vector<vector<double> >(4,vector<double> (0,0));
vector<double> years;
vector<double> months;
vector<double> julians;
vector<double> days;
vector<double> TDA_unreg;
vector<double> ICFs;

//string names[M_H]={"GCL","CHJ","LWG","LGS","LMN","IHR","MCN","JDA","TDA","BON"};


double RC_conversion_factor=1.23;

//Data Structures to form the network
set<int> shoal_reservoirs={9,28,35,36,37,38,39,40,41,42};
set<int> hyssr_reservoirs={9,28,29,30,31,32,33,35,36,37,38,39,40,41,42};
set<int> opt_reservoirs={0,1,3};
set<int> run_of_river_reservoirs={39,40,41};
vector<int> reservoirs_network ={0,15,1,3,4,17};

map<int, vector<pair<int, int>>> flow_route;
map<int, int> opt_index;
map<int, int> sim_index;
map<int,int> inverse_opt;
map<int,int> inverse_sim;

//vector<vector<double> > extra_generation = utils::loadMatrix("../data/cons_opt/Datafiles/extra_gen.txt",sim_days,1,start_row,start_col);//jan 1
vector<vector<double> > storage_opt=vector<vector<double> >(sim_days+1,vector<double> (opt_reservoirs.size(),0));
vector<vector<double> > generation_sim=vector<vector<double> >(sim_days,vector<double> (reservoirs_network.size(),0));

bool do_transpose=true;


//global declaration of rbfs
vector<pFunction_param> p_param=vector<pFunction_param>(1);
vector<param_function*> mPolicy=vector<param_function*>(1);


vector<vector<double> > min_spills=utils::readMatrixFromDataFile("../data/HYSSR/min_constraints/min_spills_opt.txt");
vector<vector<double> > max_spills=utils::readMatrixFromDataFile("../data/HYSSR/min_constraints/max_spills_opt.txt");
vector<vector<double> > min_powerh=utils::readMatrixFromDataFile("../data/HYSSR/min_constraints/min_ph_opt.txt");
vector<double> max_powerh=utils::loadVectorFromDataFile("../data/HYSSR/min_constraints/max_ph_opt.txt");
vector<double> max_discharge_gcl=utils::loadVectorFromDataFile("../data/cons_opt/Datafiles/max_GCL_discharge.txt");
//Revenue model

void makeFlowRoute()
{
        
                        //spokane flow
                        flow_route[23].push_back(make_pair(8,0));//Upper falls to Post falls
                        flow_route[24].push_back(make_pair(23,0));
                        flow_route[25].push_back(make_pair(24,0));
                        flow_route[26].push_back(make_pair(25,0));
                        flow_route[27].push_back(make_pair(26,0));

                        //upper columbia
                        flow_route[15].push_back(make_pair(0,0));//RVC to Mica
                        flow_route[1].push_back(make_pair(15,0));

                        //kootenai
                        flow_route[46].push_back(make_pair(2,0));//Bonn's Ferry to Libby
                        flow_route[4].push_back(make_pair(46,0));
                        flow_route[4].push_back(make_pair(3,0));
                        flow_route[17].push_back(make_pair(4,0));

                        //flathead 
                        flow_route[6].push_back(make_pair(5,0));//Kerr to Hungry Horse
                        flow_route[6].push_back(make_pair(48,0));
                        flow_route[18].push_back(make_pair(6,0));
                        flow_route[13].push_back(make_pair(18,1));
                        flow_route[19].push_back(make_pair(13,0));
                        flow_route[7].push_back(make_pair(45,0));
                        flow_route[7].push_back(make_pair(19,0));
                        flow_route[20].push_back(make_pair(7,1));
                        flow_route[21].push_back(make_pair(20,0));
                        flow_route[16].push_back(make_pair(21,0));
                        flow_route[22].push_back(make_pair(16,0));


                        //mid columbia
                        flow_route[9].push_back(make_pair(1,1));
                        flow_route[9].push_back(make_pair(17,1));
                        flow_route[9].push_back(make_pair(22,1));
                        flow_route[28].push_back(make_pair(9,0));
                        flow_route[29].push_back(make_pair(28,0));
                        flow_route[30].push_back(make_pair(29,0));
                        flow_route[30].push_back(make_pair(10,0));
                        flow_route[31].push_back(make_pair(30,0));
                        flow_route[32].push_back(make_pair(31,0));
                        flow_route[33].push_back(make_pair(32,0));

                        //snake river
                        flow_route[34].push_back(make_pair(11,0));//Oxbow to Borwnlee
                        flow_route[43].push_back(make_pair(34,0));
                        flow_route[35].push_back(make_pair(12,1));//LWG to Dworshak
                        flow_route[35].push_back(make_pair(43,1));
                        flow_route[35].push_back(make_pair(49,0));
                        flow_route[35].push_back(make_pair(50,0));
                        flow_route[35].push_back(make_pair(51,0));
                        flow_route[35].push_back(make_pair(52,0));
                        flow_route[35].push_back(make_pair(53,0));
                        flow_route[36].push_back(make_pair(35,0));
                        flow_route[37].push_back(make_pair(36,0));
                        flow_route[38].push_back(make_pair(37,0));

                        //lower columbia
                        flow_route[44].push_back(make_pair(14,0));//PLT to RBU
                        flow_route[39].push_back(make_pair(38,0));
                        flow_route[39].push_back(make_pair(33,0));
                        flow_route[39].push_back(make_pair(54,0));
                        flow_route[40].push_back(make_pair(39,0));
                        flow_route[41].push_back(make_pair(40,0));
                        flow_route[41].push_back(make_pair(44,0));
                        flow_route[42].push_back(make_pair(41,0));
                        do_transpose=false;
                        //cout<<"FLOW ROUTE "<<to_string(rank_no)<<" master "<<to_string(master_no)<<endl;
}



void makeOPTIndexing()
{
     opt_index[0]=0;
     opt_index[1]=1;
     opt_index[3]=2;
    
}

void makeinverseOPTIndexing()
{

      inverse_opt[0]=0;
      inverse_opt[1]=1;
      inverse_opt[2]=3;
}
void makeSimIndexing()
{

     sim_index[0]=0;
     sim_index[15]=1;
     sim_index[1]=2;
     sim_index[3]=3;
     sim_index[4]=4;
     sim_index[17]=5;
    
}

void makeinverseSimIndexing()
{

      inverse_sim[0]=0;
      inverse_sim[1]=15;
      inverse_sim[2]=1;
      inverse_sim[3]=3;
      inverse_sim[4]=4;
      inverse_sim[5]=17;
}


//vector<vector<double> > generations_constant=utils::readMatrixFromDataFile("../data/HYSSR/generations_constant.txt");
//vector<vector<double> > max_hydropower=utils::readMatrixFromDataFile("../data/HYSSR/max_hydropower.txt");

double stageHYSSR(int i, double res_stor)
{

            double res_forebay = StorEq[i][2]*pow(res_stor,4) + StorEq[i][3]*pow(res_stor,3) + StorEq[i][4]*pow(res_stor,2) + StorEq[i][5]*res_stor + StorEq[i][6];
            return res_forebay;


}

double tailwaterHYSSR(int i, double res_discharge)
{

            double res_tailwater = (TailEq[i][2]*pow(res_discharge,4) + TailEq[i][3]*pow(res_discharge,3) + TailEq[i][4]*pow(res_discharge,2) + TailEq[i][5]*res_discharge+  TailEq[i][6]);
            return res_tailwater;


}

double generationHYSSR(double res_head, double res_powerflow, int res_num)
{
    double bias_correction=0.94;
    return max(0.0,((density_water*res_head*res_powerflow*powerhouse_conversion_factor*gravity*powerhouse_efficiency[res_num])/1000000)*24)*bias_correction;
}



double interpolateStageTDA(double _stor)
{
	      
        std::ostringstream stream;
        stream << "../data/StageStorageAF/TDA_StageStorage.txt";
        std::string file_path = stream.str();
        vector<vector<double>> vec=utils::readMatrixFromDataFile(file_path);
        utils::transpose(vec);
        vector<double> _stage=vec[0];
        vector<double> _storage=vec[1];
  	return utils::interp_lin(_storage,_stage,_stor);
 	      
}



void model_storage_res(int res_no, int t,vector<double> uu)
        {
                
                        int idt=t+1;
                        
                        double s_ti=storage[res_no][t-1]+inflow[res_no][t]*cfs_kaf_conversion_factor;
                        vector<pair<int,int>> connections=flow_route[res_no];
                        for(pair<int,int> con:connections)
                        {  
                                
                                int con_num=con.first;
                                int con_time=con.second;
                                if(con_num>47)
                                {
                                        s_ti += inflow[con_num][t-con_time]*cfs_kaf_conversion_factor;
                                }
                                else
                                {
                                        s_ti += discharge[con_num][t-con_time]*cfs_kaf_conversion_factor;
                                } 
                                
                        }

                        //total_inflow[res_no][t]=(s_ti-storage[res_no][t-1])*(1/cfs_kaf_conversion_factor);
                    
                        
                        double min_discharge=max(0.0,s_ti-s_max[res_no]);
                        double max_discharge=s_ti;
                        discharge[res_no][t] = min( max_discharge , max( min_discharge , uu[opt_index[res_no]]*cfs_kaf_conversion_factor) );
                        powerflow[res_no][t] = min(discharge[res_no][t],power_max[res_no]);
                        spill[res_no][t]    = max(0.0,discharge[res_no][t]- powerflow[res_no][t]);
                        storage[res_no][t]   = (s_ti-discharge[res_no][t]);
                        discharge[res_no][t] = discharge[res_no][t]*(1/cfs_kaf_conversion_factor);
                        powerflow[res_no][t] *= (1/cfs_kaf_conversion_factor);
                        spill[res_no][t]*=(1/cfs_kaf_conversion_factor);


                /*
                       if (t<5)
                       {

                                cout<<"res_no "<<res_no<<endl;
                                cout<<"s_ti "<<s_ti<<endl;
                                cout<<"s_max[res_no] "<<s_max[res_no]<<endl;
                                cout<<"max_discharge "<<max_discharge<<endl;
                                cout<<"power_max[res_no] "<<power_max[res_no]<<endl;
                                cout<<"discharge[res_no][t] "<<discharge[res_no][t]<<endl;
                                cout<<"spill[res_no][t] "<<spill[res_no][t]<<endl;
                                cout<<"powerflow[res_no][t] "<<powerflow[res_no][t]<<endl;
                                cout<<"uu[opt_index[res_no]] "<<uu[opt_index[res_no]]<<endl;
                                cout<<"uu[0] "<<uu[0]<<endl;
                                cout<<"uu[1] "<<uu[1]<<endl;
                                cout<<"uu[2] "<<uu[2]<<endl;
                                cout<<"powerflow[res_no][t]*cfs_kaf_conversion_factor "<<powerflow[res_no][t]*cfs_kaf_conversion_factor<<endl;


                       }
                */


                        calcGeneration(res_no,t);
                        generation_sim[t-1][sim_index[res_no]]=generation[res_no][t];

        }

        void model_run_of_river(int res_no, int t)
        {

                        int idt=t+1;
                        vector<pair<int,int>> connections=flow_route[res_no];
                        double d_ti=inflow[res_no][t]*cfs_kaf_conversion_factor;
                        double s_ti=storage[res_no][t-1]+inflow[res_no][t]*cfs_kaf_conversion_factor;
                        double min_discharge=0.0;
                        //cout<<"connections size is"<<connections.size();
                        for(pair<int,int> con:connections)
                        {

                                int con_num=con.first;
                                int con_time=con.second;
                                if(con_num>47)
                                {
                                        d_ti += inflow[con_num][t-con_time]*cfs_kaf_conversion_factor;
                                        s_ti += inflow[con_num][t-con_time]*cfs_kaf_conversion_factor;
                                }
                                else
                                {
                                        d_ti += discharge[con_num][t-con_time]*cfs_kaf_conversion_factor;
                                        s_ti += discharge[con_num][t-con_time]*cfs_kaf_conversion_factor;
                                }

                        }
                        discharge[res_no][t]=max(min_discharge,d_ti)*(1/cfs_kaf_conversion_factor);
                        storage[res_no][t]=(s_ti-d_ti);
                        powerflow[res_no][t] = min(d_ti,power_max[res_no])*(1/cfs_kaf_conversion_factor);
                        spill[res_no][t] = max(0.0,d_ti - power_max[res_no])*(1/cfs_kaf_conversion_factor);
                        calcGeneration(res_no,t);
                        generation_sim[t-1][sim_index[res_no]]=generation[res_no][t];

                       /* 
                       if (t<5)
                       {

                                cout<<"res_no "<<res_no<<endl;
                                cout<<"s_ti "<<s_ti<<endl;
                                cout<<"s_max[res_no] "<<s_max[res_no]<<endl;
                                cout<<"power_max[res_no] "<<power_max[res_no]<<endl;
                                cout<<"discharge[res_no][t] "<<discharge[res_no][t]<<endl;
                                cout<<"spill[res_no][t] "<<spill[res_no][t]<<endl;
                                cout<<"powerflow[res_no][t] "<<powerflow[res_no][t]<<endl;
                                cout<<"powerflow[res_no][t]*cfs_kaf_conversion_factor "<<powerflow[res_no][t]*cfs_kaf_conversion_factor<<endl;


                       }*/

        }



        void calcGeneration(int res_no, int t)
        {

                        generation[res_no][t] = generationHYSSR(hydraulic_head[res_no], powerflow[res_no][t]*cfs_kaf_conversion_factor,res_no);


        }

void fillSpokaneFlow()
        {

                        for(int t=0;t<sim_days+2;t++)
                        {
                                spokane_flow.push_back(inflow[8][t]+inflow[23][t]+inflow[24][t]+inflow[25][t]+inflow[26][t]+inflow[27][t]);
                        }

        }


void simulate(unsigned int rank_no,unsigned int master_no)
{

        
        
        

                
                       for(int t=1;t<=sim_days;t++)
                       { 
                            vector<double> input; 
                            vector<double> uu;
                            int julian = (int)julians[t-1];

                            input.push_back( storage[0][t-1] );
                            input.push_back( storage[1][t-1] );
                            input.push_back( storage[3][t-1] );
                            input.push_back( inflow[0][t-1] );
                            input.push_back( inflow[1][t-1] );
                            input.push_back( inflow[3][t-1] );
                            input.push_back( sin( 2*pi*julian/365));
                            input.push_back( cos( 2*pi*julian/365));
                            
                        /*
                            if (t<5)
                       {

                                
                                cout<<"printing input"<<endl;
                                cout<<input.size()<<endl;
                                for (auto res : input)
                                cout<<res<<" "<<endl;
                                 cout<<"end printing input"<<endl;
                                cout<<"storage[0][t-1] "<<storage[0][t-1]<<endl;
                                cout<<"storage[1][t-1]"<<storage[1][t-1]<<endl;
                                cout<<"storage[3][t-1] "<<storage[3][t-1]<<endl;
                                cout<<"inflow[0][t-1] "<<inflow[0][t-1]<<endl;
                                cout<<"inflow[1][t-1] "<<inflow[1][t-1]<<endl;
                                cout<<"inflow[3][t-1] "<<inflow[3][t-1]<<endl;
                                cout<<"sin( 2*pi*julian/365) "<<sin( 2*pi*julian/365)<<endl;
                                cout<<"cos( 2*pi*julian/365) "<<cos( 2*pi*julian/365)<<endl;


                       }*/

                            rbfs_inputs[0][t-1]=storage[0][t-1];
                            rbfs_inputs[1][t-1]=storage[1][t-1];
                            rbfs_inputs[2][t-1]=storage[3][t-1];
                            rbfs_inputs[3][t-1]=inflow[0][t-1];
                            rbfs_inputs[4][t-1]=inflow[1][t-1];
                            rbfs_inputs[5][t-1]=inflow[3][t-1];
                            rbfs_inputs[6][t-1]=sin( 2*pi*julian/365);
                            rbfs_inputs[7][t-1]=cos( 2*pi*julian/365);

                            uu = mPolicy[0]->get_NormOutput(input);
                            /*
                                 if (t<5)
                           {
                               cout<<"printing uu"<<endl;
                               cout<<uu.size()<<endl;
                               cout<<"end printing uu"<<endl;
                            }
                            */
                            rbfs_releases[0][t-1]=uu[0];
                            rbfs_releases[1][t-1]=uu[1];
                            rbfs_releases[2][t-1]=uu[2];
                            for (auto res : reservoirs_network)
                            {
                            if(opt_reservoirs.find(res)!=opt_reservoirs.end())
                                model_storage_res(res,t,uu);
                            else
                                model_run_of_river(res,t);
                            }
                            
                
                         }


}



void initialize()
{
    fillSpokaneFlow();
    //total_inflow[9][0]=modelled_discharge[1][0]+modelled_discharge[17][0]+modelled_discharge[22][0]+spokane_flow[0]+inflow[9][1];
    //total_inflow[35][0]=modelled_discharge[12][0]+modelled_discharge[43][0]+inflow[35][1]+ inflow[49][1]+ inflow[50][1]+ inflow[51][1]+ inflow[52][1]+ inflow[53][1];
}

void initializeGCL()
{

        //Pre-processing the calendar data structures for modelling GCL
        utils::transpose(calender);
        for(int i=0;i<(sim_days+1);i++)
        {
                  for(int j=0;j<4;j++)
                  {
                          c[j].insert(c[j].end(), calender[j+1].begin(), calender[j+1].end());
                          if(j==3)
                          std::transform(c[j].end()-365, c[j].end(), c[j].end()-365,[=](double k) { return k+i; });
                          
                  }
        }
    
        
        months = c[0];
        days = c[1];
        julians = c[2];
        years = c[3];

        TDA_unreg=inflow[47];
        
        vector<vector<double> > ICF_file=utils::loadMatrix("../data/optimize_four/ICF_48.txt",sim_days/365,1,1,0);
        utils::transpose(ICF_file);
        ICFs=ICF_file[0];

}

void preprocess()
{

        
        utils::transpose(modelled_discharge);
        utils::transpose(modelled_storage);
        utils::transpose(modelled_generation);
        utils::transpose(modelled_powerflow);
        utils::transpose(inflow);
        utils::transpose(maxhydro_matrix);
        utils::transpose(maxstorage_matrix);
        utils::transpose(power_constants);
        utils::transpose(discharge);
        utils::transpose(storage);
        power_max=maxhydro_matrix[1];
        s_max=maxstorage_matrix[1];
        hydraulic_head=power_constants[1];
        powerhouse_efficiency=power_constants[2];
        initializeGCL();
        
}

void writeSimulationOutput(unsigned int rank_no,unsigned int master_no)
{

        utils::writeMatrix(generation,"../output_opt/policy"+to_string(rank_no)+"/generation_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(storage,"../output_opt/policy"+to_string(rank_no)+"/storage_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(discharge,"../output_opt/policy"+to_string(rank_no)+"/discharge_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(powerflow,"../output_opt/policy"+to_string(rank_no)+"/powerflow_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(spill,"../output_opt/policy"+to_string(rank_no)+"/spill_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(heads,"../output_opt/policy"+to_string(rank_no)+"/head_hyssr"+to_string(master_no)+".txt", sim_days+1, M_H);
        utils::writeMatrix(stage_hyssr,"../output_opt/policy"+to_string(rank_no)+"/stage_hyssr"+to_string(master_no)+".txt", sim_days+1, M_H);
        utils::writeMatrix(tailwater_hyssr,"../output_opt/policy"+to_string(rank_no)+"/tailwater_hyssr"+to_string(master_no)+".txt", sim_days+1, M_H);
        utils::writeMatrix(rbfs_releases,"../output_opt/policy"+to_string(rank_no)+"/rbfs_releases"+to_string(master_no)+".txt", sim_days, 3);
        utils::writeMatrix(rbfs_inputs,"../output_opt/policy"+to_string(rank_no)+"/rbfs_inputs"+to_string(master_no)+".txt", sim_days, 8);
        utils::writeMatrix(modelled_generation,"../output_opt/policy"+to_string(rank_no)+"/modelled_generation"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(modelled_discharge,"../output_opt/policy"+to_string(rank_no)+"/modelled_discharge"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(modelled_storage,"../output_opt/policy"+to_string(rank_no)+"/modelled_storage"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(modelled_powerflow,"../output_opt/policy"+to_string(rank_no)+"/modelled_powerflow"+to_string(master_no)+".txt", sim_days+1, M);
        //cout<<"write simulatio output "<<to_string(rank_no)+" "<<to_string(master_no)<<endl;
        utils::transpose(generation_sim);
        utils::writeMatrix(generation_sim,"../output_opt/policy"+to_string(rank_no)+"/generation_sim"+to_string(master_no)+".txt", sim_days, 6);
        utils::transpose(generation_sim);

}



vector<double> maxRenewablesGeneration()
{

        vector<double> renewables(366,0);
        for(int t=0;t<sim_days;t++)
        {
        
            int julian=(int)julians[t];
            renewables[julian]+=utils::columnSumRow(generation_sim,t);
        }
        std::transform(renewables.begin(), renewables.end(), renewables.begin(),[=](double i) { return i/sim_years; });
        renewables.erase (renewables.begin());
        return renewables;

}




vector<double> getObjectives(unsigned int rank_no,unsigned int master_no)
{
        simulate(rank_no,master_no);
        if(write_output)
        writeSimulationOutput(rank_no,master_no);
        vector<double> J;
        
        J.push_back(-1*meanVector(maxRenewablesGeneration()));
        
        return J;


}

void initializeRBFs(unsigned int rank_no,unsigned int master_no)
        {
                if(do_transpose)
                {
                preprocess();
                makeSimIndexing();
                makeinverseSimIndexing();
                makeOPTIndexing();
                makeinverseOPTIndexing();
                initialize();
                makeFlowRoute();
                //cout<<"INITIALIZE RBFS "<<to_string(rank_no)<<" "<<to_string(master_no)<<endl;
                }
                
                vector<vector<double> > min_max_inflow_rbfs_range=utils::readMatrixFromDataFile("../data/canadian_dams/inflows_rbfs_ranges.txt");
                //cout<<to_string(rank_no)+"inflow range size "<<to_string(master_no)<< "min_max_inflow_rbfs_range[0]"<< min_max_inflow_rbfs_range[5][1]<<endl;
                vector<vector<double> > min_max_releases_range=utils::readMatrixFromDataFile("../data/canadian_dams/discharge_rbfs_ranges.txt");
                for(int i=0;i<1;i++)
                {   
                p_param[i].tPolicy=1;
                p_param[i].policyInput=8;
                p_param[i].policyOutput=3;
                p_param[i].policyStr=12;
                param_function* mp1 = new rbf(p_param[i].policyInput,p_param[i].policyOutput,p_param[i].policyStr);
                //param_function* mp1 = new ncRBF(p_param[i].policyInput,p_param[i].policyOutput,p_param[i].policyStr);
                mPolicy[i] = mp1;
                }
                for(int i=0;i<1;i++)
                {   
                
                        
                        p_param[i].mIn.push_back(0);
                        p_param[i].mIn.push_back(0);
                        p_param[i].mIn.push_back(0);
                        p_param[i].mIn.push_back(min_max_inflow_rbfs_range[0][0]);
                        p_param[i].mIn.push_back(min_max_inflow_rbfs_range[1][0]);
                        p_param[i].mIn.push_back(min_max_inflow_rbfs_range[2][0]);
                        p_param[i].mIn.push_back(-1);
                        p_param[i].mIn.push_back(-1);
                        
                        p_param[i].MIn.push_back(s_max[0]);
                        p_param[i].MIn.push_back(s_max[1]);
                        p_param[i].MIn.push_back(s_max[3]);
                        p_param[i].MIn.push_back(min_max_inflow_rbfs_range[0][1]);
                        p_param[i].MIn.push_back(min_max_inflow_rbfs_range[1][1]);
                        p_param[i].MIn.push_back(min_max_inflow_rbfs_range[2][1]);
                        p_param[i].MIn.push_back(1);
                        p_param[i].MIn.push_back(1);

                        p_param[i].mOut.push_back(0);
                        p_param[i].mOut.push_back(0);
                        p_param[i].mOut.push_back(0);
                        p_param[i].MOut.push_back(min_max_releases_range[0][1]);
                        p_param[i].MOut.push_back(min_max_releases_range[1][1]);
                        p_param[i].MOut.push_back(min_max_releases_range[2][1]);
                        
                        mPolicy[i]->setMinInput(p_param[i].mIn); 
                        mPolicy[i]->setMinOutput(p_param[i].mOut);
                        mPolicy[i]->setMaxInput(p_param[i].MIn); 
                        mPolicy[i]->setMaxOutput(p_param[i].MOut);
                        
                }
                //cout<<"after initialize"<<endl;
        }


void evaluate(double* var, double* obj, double* constrs, unsigned int inp_rank_no,unsigned int inp_master_no){

        unsigned int rank_no= inp_rank_no,master_no=inp_master_no;
        int decsvars=204;
        initializeRBFs(rank_no,master_no);
        
        unsigned int p2 = 0;
        int getArrayLength;
        for(int k=0;k<1;k++)
                {
                
                double fullVar[p_param[k].policyStr*(2*p_param[k].policyInput+p_param[k].policyOutput)];
                //double fullVar[p_param[k].policyStr*(2*p_param[k].policyInput+p_param[k].policyOutput)];
                unsigned int p1 = 0;
                
                //ncRBF constant
                /*
                if (var[183]<0.5)
                {
                
                        for(unsigned int l = 0; l < p_param[k].policyOutput; l++)
                        {
                                fullVar[p1++] = 0.0;
                                p2++;
                        }
                
                }
                else
                {
                        for(unsigned int l = 0; l < p_param[k].policyOutput; l++)
                        {
                                fullVar[p1++] = var[p2++];
                        }
                }
                */
                /*
                for(unsigned int l = 0; l < p_param[k].policyOutput; l++)
                        {
                                fullVar[p1++] = 0.0;
                                p2++;
                        }

                */


                for(unsigned int i = 0; i < p_param[k].policyStr; i++){
                
                        for(unsigned int j = 0; j < p_param[k].policyInput-2; j++)
                        {
                                fullVar[p1++] = var[p2++];
                                fullVar[p1++] = var[p2++];
                        }
                        for(unsigned int j = 0; j < 2 ; j++)
                        {
                                // set centers to 0 and radii to 1 for sin and cos inputs
                                fullVar[p1++] = 0.0;
                                fullVar[p1++] = p2++;
                        }
                        for(unsigned int l = 0; l < p_param[k].policyOutput; l++)
                        {
                                fullVar[p1++] = var[p2++];
                        }
                        
                }
                
                /*
                int getArrayLength = sizeof(fullVar) / sizeof(double);
                cout <<"getArrayLength "<<getArrayLength<<endl;

                cout<<"printing fullVar"<<endl;
                for(int i=0;i<getArrayLength;i++)
                        cout<<fullVar[i]<<" "<<endl;
                cout<<"end printing fullVar"<<endl;
                */
                getArrayLength = sizeof(fullVar) / sizeof(double);
                vector<double> full_vars_vector(getArrayLength,0);
                    for(int i=0;i<getArrayLength;i++)
                         full_vars_vector[i]=fullVar[i];


                    if(write_output)
                {
                  vector<vector<double> >full_vars_matrix(1,vector<double> (getArrayLength,0));
                    full_vars_matrix[0]=full_vars_vector;
                    utils::writeMatrix(full_vars_matrix,"../output_opt/policy"+to_string(rank_no)+"/full_decs_vars"+to_string(master_no)+".txt", getArrayLength, 1);
                }

                mPolicy[k]->setParameters(fullVar);
                }
                
                
                vector<double> J ;
                J = getObjectives(rank_no,master_no);
                //cout<<"after get objectives"<<endl;
                for(unsigned int i=0; i<Nobj; i++){
                obj[i] = J[i];
                }
                

                if(write_output)
                {
                    vector<vector<double> >obj_matrix(1,vector<double> (Nobj,0));
                    obj_matrix[0]=J;
                    utils::writeMatrix(obj_matrix,"../output_opt/policy"+to_string(rank_no)+"/objectives"+to_string(master_no)+".txt", Nobj, 1);
                                    
                    
                    vector<double> vars_vector(decsvars,0);
                    for(int i=0;i<decsvars;i++)
                         vars_vector[i]=var[i];
                         
                    
                    
                    vector<vector<double> >vars_matrix(1,vector<double> (decsvars,0));
                    vars_matrix[0]=vars_vector;
                    utils::writeMatrix(vars_matrix,"../output_opt/policy"+to_string(rank_no)+"/decs_vars"+to_string(master_no)+".txt", decsvars, 1);
                }
                
                for(int k=0;k<1;k++)
                  mPolicy[k]->clearParameters();



        }


};




/*
void run_model()
{

         preprocess();
         makeHYSSRIndexing();
         makeinverseHYSSRIndexing();
         makeOPTIndexing();
         makeinverseOPTIndexing();
         makeFlowRoute();
         //initialize();
         getObjectives();

}
*/
/*
int main()
{
         
         preprocess();
         
         makeHYSSRIndexing();
         makeinverseHYSSRIndexing();
         makeOPTIndexing();
         makeinverseOPTIndexing();
         makeFlowRoute();
         //initialize();
         //writeSimulationOutput();
         //getObjectives();
         

         /*
         vector<double> arima_caiso_coeff=utils::loadVectorFromDataFile("../data/cons_opt/models/midc_arima_coef.txt");
         vector<vector<double> > residuals_matrix=utils::readMatrixFromDataFile("../data/cons_opt/models/MidC_cons_residuals_cpp.txt");
         cout<<arima_caiso_coeff.size()<<endl;
         cout<<residuals_matrix[0].size()<<endl;
         cout<<residuals_matrix.size()<<endl;
         cout<<storage.size()<<endl;
         cout<<storage[0].size()<<endl;
         for( int i=0;i<10;i++)
            cout<<storage[9][i]<<endl;

         
         return 0;
}
*/