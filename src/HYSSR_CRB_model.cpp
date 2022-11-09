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
const int M_H=14;
const int M_HR=30;
int sim_years=10;
int sim_days=sim_years*365;
int no_years=sim_years+2;
int bpa_rev_sim_years=sim_years-1;
int bpa_rev_sim_days=sim_days-365;
int Nobj=4;
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
double sdNegativeSumOpt;
vector<vector<double> > modelled_discharge=utils::loadMatrix("../data/optimize_four/discharge_48.txt",sim_days+1,M,365-1,start_col);//water year
//vector<vector<double> > modelled_discharge=utils::loadMatrix("../data/optimize_four/mod_discharge_sim_arrow.txt",sim_days+1,M,0,start_col);//arrow changed
//vector<vector<double> > modelled_discharge=utils::loadMatrix("../data/optimize_four/mod_discharge_sim_arrow_dun.txt",sim_days+1,M,0,start_col);//arrow changed

vector<vector<double> > modelled_storage=utils::loadMatrix("../data/optimize_four/storage_48.txt",sim_days+1,M,365-1,start_col);//water year
vector<vector<double> > modelled_generation=utils::loadMatrix("../data/optimize_four/generation_48.txt",sim_days+1,M+1,365-1,start_col);//water year
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

//Data structure to model the prices
vector<vector<double> > modelled_demands    = utils::loadMatrix("../data/optimize_four/demands_surrogate.txt",sim_days,15,start_row,start_col);//jan 1
vector<vector<double> > modelled_solar      = utils::loadMatrix("../data/optimize_four/solar_surrogate.txt",sim_days,1,start_row,start_col);//jan 1
vector<vector<double> > modelled_wind       = utils::loadMatrix("../data/optimize_four/wind_surrogate.txt",sim_days,1,start_row,start_col);//jan 1
vector<vector<double> > modelled_cahydro    = utils::loadMatrix("../data/optimize_four/CA_hydro_surrogate.txt",sim_days,1,start_row,start_col);//jan 1
vector<vector<double> > extra_generation = utils::loadMatrix("../data/optimize_four/extra_gen.txt",sim_days,1,start_row,start_col);//jan 1

map<int,int> PF_rates_months;
void makePF_rates_monthsIndexing()
{
     PF_rates_months[10]=0;
     PF_rates_months[11]=1;
     PF_rates_months[12]=2;
     PF_rates_months[1]=3;
     PF_rates_months[2]=4;
     PF_rates_months[3]=5;
     PF_rates_months[4]=6;
     PF_rates_months[5]=7;
     PF_rates_months[6]=8;
     PF_rates_months[7]=9;
     PF_rates_months[8]=10;
     PF_rates_months[9]=11;
    
}
double calculate_CRAC(double NR_, double tot_load)
  {
    double NR1,NR2,X;
    if(NR_ > 5*pow(10,6))
        {  
                if(NR_ > 100*pow(10,6))
                {
                    NR1=100*pow(10,6);
                    NR2=(NR_ - 100*pow(10,6))/2;
                }
        
                else
                {
                     NR1 = NR_;
                     NR2= 0;
                 }
                X=min((NR1+NR2)/(tot_load*24) ,  300*pow(10,6)/(tot_load*24));
        }
    else
        X=0;
    return X;
 }


vector<vector<double> > storage=modelled_storage;
vector<vector<double> > total_inflow=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > generation=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > spill=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > powerflow=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > tailwater=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > stage_hyssr=vector<vector<double> >(M_H,vector<double> (sim_days+1,0));
vector<vector<double> > tailwater_hyssr=vector<vector<double> >(M_H,vector<double> (sim_days+1,0));
vector<vector<double> > heads=vector<vector<double> >(M_H,vector<double> (sim_days+1,0));
vector<vector<double> > rbfs_releases=vector<vector<double> >(4,vector<double> (sim_days,0));
vector<vector<double> > rbfs_inputs=vector<vector<double> >(9,vector<double> (sim_days,0));

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
set<int> shoal_reservoirs={2,5,7,9,12,28,35,36,37,38,39,40,41,42};
set<int> hyssr_reservoirs={2,5,7,9,12,28,29,30,31,32,33,35,36,37,38,39,40,41,42,46,4,17,6,18,13,19,20,21,16,22};
set<int> spills_reservoirs={9,28,35,36,37,38,42};
set<int> opt_reservoirs={2,5,9,12};
set<int> run_of_river_reservoirs={39,40,41};
vector<int> reservoirs_network ={2,46,4,17,5,6,18,13,19,7,20,21,16,22,9,28,29,30,31,32,33,12,35,36,37,38,39,40,41,42};

map<int, vector<pair<int, int>>> flow_route;
map<int, int> opt_index;
map<int, int> sim_index;
map<int, int> hyssr_index;
map<int,int> inverse_opt;
map<int,int> inverse_sim;
map<int,int> inverse_hyssr;

//vector<vector<double> > extra_generation = utils::loadMatrix("../data/cons_opt/Datafiles/extra_gen.txt",sim_days,1,start_row,start_col);//jan 1
vector<vector<double> > storage_opt=vector<vector<double> >(sim_days+1,vector<double> (opt_reservoirs.size(),0));
vector<vector<double> > generation_sim=vector<vector<double> >(sim_days,vector<double> (reservoirs_network.size(),0));
vector<vector<double> > generation_fcrps=vector<vector<double> >(sim_days,vector<double> (hyssr_reservoirs.size(),0));
vector<vector<double> > generation_shoal=vector<vector<double> >(sim_days,vector<double> (shoal_reservoirs.size(),0));

bool do_transpose=true;


//global declaration of rbfs
vector<pFunction_param> p_param=vector<pFunction_param>(1);
vector<param_function*> mPolicy=vector<param_function*>(1);
double phi;


vector<vector<double> > min_spills=utils::readMatrixFromDataFile("../data/HYSSR/min_constraints/min_spills_opt.txt");
vector<vector<double> > max_spills=utils::readMatrixFromDataFile("../data/HYSSR/min_constraints/max_spills_opt.txt");
vector<vector<double> > min_powerh=utils::readMatrixFromDataFile("../data/HYSSR/min_constraints/min_ph_opt.txt");
vector<double> max_powerh=utils::loadVectorFromDataFile("../data/HYSSR/min_constraints/max_ph_opt.txt");
vector<double> max_discharge_alb=utils::loadVectorFromDataFile("../data/ALB/alb_max_discharge.txt");
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
     opt_index[2]=0;
     opt_index[5]=1;
     opt_index[9]=2;
     opt_index[12]=3;
}

void makeinverseOPTIndexing()
{

      inverse_opt[0]=2;
      inverse_opt[1]=5;
      inverse_opt[2]=9;
      inverse_opt[3]=12;
}
void makeSimIndexing()
{

     sim_index[9]=0;
     sim_index[28]=1;
     sim_index[35]=2;
     sim_index[36]=3;
     sim_index[37]=4;
     sim_index[38]=5;
     sim_index[42]=6;
    
}

void makeinverseSimIndexing()
{

      inverse_sim[0]=9;
      inverse_sim[1]=28;
      inverse_sim[2]=35;
      inverse_sim[3]=36;
      inverse_sim[4]=37;
      inverse_sim[5]=38;
      inverse_sim[6]=42;
}
void makeHYSSRIndexing()
{
     hyssr_index[2]=0;
     hyssr_index[5]=1;
     hyssr_index[7]=2;
     hyssr_index[9]=3;
     hyssr_index[12]=4;
     hyssr_index[28]=5;
     hyssr_index[35]=6;
     hyssr_index[36]=7;
     hyssr_index[37]=8;
     hyssr_index[38]=9;
     hyssr_index[39]=10;
     hyssr_index[40]=11;
     hyssr_index[41]=12;
     hyssr_index[42]=13;
     hyssr_index[29]=14;
     hyssr_index[30]=15;
     hyssr_index[31]=16;
     hyssr_index[32]=17;
     hyssr_index[33]=18;    
     hyssr_index[46]=19;
     hyssr_index[4]=20;
     hyssr_index[17]=21;
     hyssr_index[6]=22;
     hyssr_index[18]=23;
     hyssr_index[13]=24;
     hyssr_index[19]=25;
     hyssr_index[20]=26;
     hyssr_index[21]=27;
     hyssr_index[16]=28;
     hyssr_index[22]=29;

}

void makeinverseHYSSRIndexing()
{

      inverse_hyssr[0]=2;
      inverse_hyssr[1]=5;
      inverse_hyssr[2]=7;
      inverse_hyssr[3]=9;
      inverse_hyssr[4]=12;
      inverse_hyssr[5]=28;
      inverse_hyssr[6]=35;
      inverse_hyssr[7]=36;
      inverse_hyssr[8]=37;
      inverse_hyssr[9]=38;
      inverse_hyssr[10]=39;
      inverse_hyssr[11]=40;
      inverse_hyssr[12]=41;
      inverse_hyssr[13]=42;
      inverse_hyssr[14]=29;
      inverse_hyssr[15]=30;
      inverse_hyssr[16]=31;
      inverse_hyssr[17]=32;
      inverse_hyssr[18]=33;
      inverse_hyssr[19]=46;
      inverse_hyssr[20]=4;
      inverse_hyssr[21]=17;
      inverse_hyssr[22]=6;
      inverse_hyssr[23]=18;
      inverse_hyssr[24]=13;
      inverse_hyssr[25]=19;
      inverse_hyssr[26]=20;
      inverse_hyssr[27]=21;
      inverse_hyssr[28]=16;
      inverse_hyssr[29]=22;
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

void modelALB(int t)
        {

                vector<vector<double> > ALB_fc=utils::readMatrixFromDataFile("../data/ALB/ALB_fc.txt");
                vector<vector<double> > ALB_ARC=utils::readMatrixFromDataFile("../data/ALB/ALB_ARC.txt");
                vector<vector<double> > ALB_CRC=utils::readMatrixFromDataFile("../data/ALB/ALB_CRC.txt");
                
                int julian = (int)julians[t-1];
                double ALB_upper = ALB_fc[julian-1][0];
                double ALB_lower = ALB_fc[julian-1][1];
                
                vector<double> ALB_sfore(sim_days+1,0);
                vector<double> ALB_ORC(365,0);
                double s_ti=storage[7][t-1];
                ALB_sfore[t]=min(((inflow[7][t]-4000)+ discharge[19][t] + discharge[45][t])*cfs_kaf_conversion_factor+ s_ti,s_max[7]);
    
                double ALB_VRC= ALB_sfore[t];
                if(julian>=213 && julian <=365)
                {
                    ALB_ORC[julian-1]= max(ALB_CRC[julian-1][1],ALB_ARC[julian-1][1]);
                    ALB_ORC[julian-1]= min(ALB_ORC[julian-1],ALB_upper);
                    ALB_ORC[julian-1]= max(ALB_ORC[julian-1],ALB_lower);
                }
                else if(julian >=1 && julian <= 212)
                {
                    ALB_ORC[julian-1]= min(ALB_VRC,max(ALB_CRC[julian-1][1],ALB_ARC[julian-1][1]));
                    ALB_ORC[julian-1]= min(ALB_ORC[julian-1],ALB_upper);
                    ALB_ORC[julian-1]= max(ALB_ORC[julian-1],ALB_lower);
                }
    
                double d_ti;
                d_ti = max(4000*cfs_kaf_conversion_factor, s_ti +(inflow[7][t]+ discharge[19][t] + discharge[45][t])*cfs_kaf_conversion_factor - ALB_ORC[julian-1]);
                d_ti=min(d_ti,max_discharge_alb[months[t-1]-1]*cfs_kaf_conversion_factor);

               double overflow = max(0.0,d_ti - power_max[7]);
               double stor_avail = ALB_upper - s_ti;
    
                if(overflow>0 and stor_avail<=0)
                {
                    if(s_ti + (inflow[7][t]+ discharge[19][t] + discharge[45][t])*cfs_kaf_conversion_factor - d_ti > s_max[7])
                        d_ti = s_ti + (inflow[7][t]+ discharge[19][t] + discharge[45][t])*cfs_kaf_conversion_factor - s_max[7];
                    powerflow[7][t] = power_max[7];
                    spill[7][t] = d_ti - power_max[7];
                }
                else if(overflow>0 and stor_avail>0)
                {
                    d_ti = max(d_ti - stor_avail,power_max[7]);
                    if(s_ti + (inflow[7][t]+ discharge[19][t] + discharge[45][t])*cfs_kaf_conversion_factor - d_ti > s_max[7])
                        d_ti = s_ti + (inflow[7][t]+ discharge[19][t] + discharge[45][t])*cfs_kaf_conversion_factor - s_max[7];
                    powerflow[7][t] = power_max[7];
                    spill[7][t] = max(0.0,d_ti - power_max[7]);
                }
                else
                {
                    if(s_ti + (inflow[7][t]+ discharge[19][t] + discharge[45][t])*cfs_kaf_conversion_factor - d_ti > s_max[7])
                        d_ti = s_ti + (inflow[7][t]+ discharge[19][t] + discharge[45][t])*cfs_kaf_conversion_factor - s_max[7];
                    powerflow[7][t] = d_ti;
                }

            storage[7][t]=s_ti +(inflow[7][t]+ discharge[19][t] + discharge[45][t])*cfs_kaf_conversion_factor-d_ti;
            discharge[7][t]=d_ti*(1/cfs_kaf_conversion_factor);
            powerflow[7][t] *= (1/cfs_kaf_conversion_factor);
            spill[7][t]*=(1/cfs_kaf_conversion_factor);
            calcGeneration(7,t);

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

                        if (res_no==9)
                                s_ti += spokane_flow[t-1]*cfs_kaf_conversion_factor;

                        total_inflow[res_no][t]=(s_ti-storage[res_no][t-1])*(1/cfs_kaf_conversion_factor);
                    
                        
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
                                cout<<"discharge[res_no][t] "<<res_no<<endl;
                                cout<<"spill[res_no][t] "<<res_no<<endl;



                       }
                        */


                        calcGeneration(res_no,t);
                        //generation_sim[t-1][sim_index[res_no]]=generation[res_no][t];

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
                                /*
                                if (t<5)
                                {

                                cout<<"res_no "<<res_no<<endl;
                                cout<<"con_num "<<con_num<<endl;
                                cout<<"con_time "<<s_max[res_no]<<endl;
                                cout<<"power_max[res_no] "<<power_max[res_no]<<endl;
                                cout<<"discharge[res_no][t] "<<res_no<<endl;
                                cout<<"spill[res_no][t] "<<res_no<<endl;



                               }*/

                        }
                        discharge[res_no][t]=max(min_discharge,d_ti)*(1/cfs_kaf_conversion_factor);
                        storage[res_no][t]=(s_ti-d_ti);
                        powerflow[res_no][t] = min(d_ti,power_max[res_no])*(1/cfs_kaf_conversion_factor);
                        spill[res_no][t] = max(0.0,d_ti - power_max[res_no])*(1/cfs_kaf_conversion_factor);
                        calcGeneration(res_no,t);
                        //generation_sim[t-1][sim_index[res_no]]=generation[res_no][t];

                        /*
                        if (t<5)
                       {

                                cout<<"res_no "<<res_no<<endl;
                                cout<<"s_ti "<<s_ti<<endl;
                                cout<<"s_max[res_no] "<<s_max[res_no]<<endl;
                                cout<<"power_max[res_no] "<<power_max[res_no]<<endl;
                                cout<<"discharge[res_no][t] "<<res_no<<endl;
                                cout<<"spill[res_no][t] "<<res_no<<endl;



                       }*/

        }



        void calcGeneration(int res_no, int t)
        {

                        if(hydraulic_head[res_no]==-4)
                {
                          double res_stor =storage[res_no][t]*1000;
                          double smax = s_max[res_no]*1000;
                          double res_forebay = stageHYSSR(res_no,res_stor);
                          double res_tailwater = tailwaterHYSSR(res_no,discharge[res_no][t]);
                          //cout<<"res_tailwater_is"<<res_tailwater<<endl;
                          //Convert head to meters.
                          if(res_no==41)
                          res_forebay = 13.891*log(res_stor) - 20.692;
                          if(shoal_reservoirs.find(res_no)!=shoal_reservoirs.end())
                          {
                          stage_hyssr[hyssr_index[res_no]][t]=res_forebay;
                          tailwater_hyssr[hyssr_index[res_no]][t]=res_tailwater;
                          heads[hyssr_index[res_no]][t] = max(0.0,(res_forebay - res_tailwater)*head_metre_conversion);
                          generation[res_no][t] = generationHYSSR(heads[hyssr_index[res_no]][t], powerflow[res_no][t]*cfs_kaf_conversion_factor,res_no);
                          }
                  
                }
                  else
                {
                          generation[res_no][t] = generationHYSSR(hydraulic_head[res_no], powerflow[res_no][t]*cfs_kaf_conversion_factor,res_no);
                }
                if( hyssr_reservoirs.find(res_no)!=hyssr_reservoirs.end())
                       generation_fcrps[t-1][hyssr_index[res_no]]=generation[res_no][t];
                if(shoal_reservoirs.find(res_no)!=shoal_reservoirs.end())
                      generation_shoal[t-1][hyssr_index[res_no]]=generation[res_no][t];


        }

void fillSpokaneFlow()
        {

                        for(int t=0;t<sim_days+1;t++)
                        {
                                spokane_flow.push_back(inflow[8][t]+inflow[23][t]+inflow[24][t]+inflow[25][t]+inflow[26][t]+inflow[27][t]);
                        }

        }


void simulate(unsigned int rank_no,unsigned int master_no)
{

                        vector<double> input; 
                        vector<double> uu;
                        //cout<<"inside simulate"<<rank_no<<" "<<master_no<<endl;
        

                
                       for(int t=1;t<=sim_days;t++)
                       { 
                            
                            int julian = (int)julians[t-1];

                            input.push_back( storage[2][t-1] );
                            input.push_back( storage[5][t-1] );
                            input.push_back( storage[9][t-1] );
                            input.push_back( storage[12][t-1] );
                            input.push_back( inflow[2][t-1] );
                            input.push_back( inflow[5][t-1] );
                            input.push_back( inflow[9][t-1] );
                            input.push_back( inflow[12][t-1] );
                            input.push_back( sin( (2*pi*julian/365)-phi));
                            
                        

                            rbfs_inputs[0][t-1]=storage[2][t-1];
                            rbfs_inputs[1][t-1]=storage[5][t-1];
                            rbfs_inputs[2][t-1]=storage[9][t-1];
                            rbfs_inputs[3][t-1]=storage[12][t-1];
                            rbfs_inputs[4][t-1]=inflow[2][t-1];
                            rbfs_inputs[5][t-1]=inflow[5][t-1];
                            rbfs_inputs[6][t-1]=inflow[9][t-1];
                            rbfs_inputs[7][t-1]=inflow[12][t-1];
                            rbfs_inputs[8][t-1]=sin( (2*pi*julian/365)-phi);
                        





                            uu = mPolicy[0]->get_NormOutput(input);
                            //cout<<"uu.size() "<<uu.size()<<" "<<rank_no<<" "<<master_no<<endl; 
                            rbfs_releases[0][t-1]=uu[0];
                            rbfs_releases[1][t-1]=uu[1];
                            rbfs_releases[2][t-1]=uu[2];
                            rbfs_releases[3][t-1]=uu[3];

                            for (auto res : reservoirs_network)
                            {
                            if(opt_reservoirs.find(res)!=opt_reservoirs.end())
                                model_storage_res(res,t,uu);
                            else
                                {
                                    if(res==7)
                                        modelALB(t);
                                    else
                                        model_run_of_river(res,t);
                                }
                            }
                             //cout<<"input.size() "<<input.size()<<" "<<rank_no<<" "<<master_no<<endl;
                             input.clear();
                             uu.clear();


                         }


}



void initialize()
{
    fillSpokaneFlow();
    total_inflow[9][0]=modelled_discharge[1][0]+modelled_discharge[17][0]+modelled_discharge[22][0]+spokane_flow[0]+inflow[9][1];
    total_inflow[35][0]=modelled_discharge[12][0]+modelled_discharge[43][0]+inflow[35][1]
                        + inflow[49][1]+ inflow[50][1]+ inflow[51][1]+ inflow[52][1]+ inflow[53][1];
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
        utils::transpose(inflow);
        utils::transpose(maxhydro_matrix);
        utils::transpose(maxstorage_matrix);
        utils::transpose(power_constants);
        utils::transpose(discharge);
        utils::transpose(storage);
        utils::transpose(modelled_demands);
        utils::transpose(modelled_cahydro);
        utils::transpose(modelled_wind);
        utils::transpose(modelled_solar);
        utils::transpose(extra_generation);
        power_max=maxhydro_matrix[1];
        s_max=maxstorage_matrix[1];
        hydraulic_head=power_constants[1];
        powerhouse_efficiency=power_constants[2];
        initializeGCL();
        
}

void writeSimulationOutput(unsigned int rank_no,unsigned int master_no)
{

        cout<<"inside writing simulate"<<endl;
        utils::writeMatrix(generation,"../output_opt/policy"+to_string(rank_no)+"/generation_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(storage,"../output_opt/policy"+to_string(rank_no)+"/storage_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(discharge,"../output_opt/policy"+to_string(rank_no)+"/discharge_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(powerflow,"../output_opt/policy"+to_string(rank_no)+"/powerflow_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(spill,"../output_opt/policy"+to_string(rank_no)+"/spill_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(heads,"../output_opt/policy"+to_string(rank_no)+"/head_hyssr"+to_string(master_no)+".txt", sim_days+1, M_H);
        utils::writeMatrix(stage_hyssr,"../output_opt/policy"+to_string(rank_no)+"/stage_hyssr"+to_string(master_no)+".txt", sim_days+1, M_H);
        utils::writeMatrix(tailwater_hyssr,"../output_opt/policy"+to_string(rank_no)+"/tailwater_hyssr"+to_string(master_no)+".txt", sim_days+1, M_H);
        utils::writeMatrix(rbfs_releases,"../output_opt/policy"+to_string(rank_no)+"/rbfs_releases"+to_string(master_no)+".txt", sim_days, 4);
        utils::writeMatrix(rbfs_inputs,"../output_opt/policy"+to_string(rank_no)+"/rbfs_inputs"+to_string(master_no)+".txt", sim_days, 10);
        utils::writeMatrix(modelled_discharge,"../output_opt/policy"+to_string(rank_no)+"/modelled_discharge"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(modelled_storage,"../output_opt/policy"+to_string(rank_no)+"/modelled_storage"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(modelled_generation,"../output_opt/policy"+to_string(rank_no)+"/modelled_generation"+to_string(master_no)+".txt", sim_days+1, M+1);
        //cout<<"write simulatio output "<<to_string(rank_no)+" "<<to_string(master_no)<<endl;
        utils::transpose(generation_fcrps);
        utils::writeMatrix(generation_fcrps,"../output_opt/policy"+to_string(rank_no)+"/generation_fcrps"+to_string(master_no)+".txt", sim_days, M_HR);
        utils::transpose(generation_fcrps);
        utils::transpose(generation_shoal);
        utils::writeMatrix(generation_shoal,"../output_opt/policy"+to_string(rank_no)+"/generation_shoal"+to_string(master_no)+".txt", sim_days, M_H);
        utils::transpose(generation_shoal);

}

vector<double> minEnvironmentConstraints()
{
        vector<double> envConstraints(366,0);
        for(int t=0;t<sim_days;t++)
        {
            int julian=(int)julians[t];
            double sum=0;
            for(int i=0;i<spills_reservoirs.size();i++)
            {
                    double spills_kcfs=spill[inverse_sim[i]][t+1]/1000;
                    if(spills_kcfs<min_spills[months[t]-1][i])
                      sum+=pow(spills_kcfs-min_spills[months[t]-1][i],2);
                    else if(spills_kcfs>max_spills[months[t]-1][i])
                      sum+=pow(spills_kcfs-max_spills[months[t]-1][i],2);
                    else
                      sum+=0;
            }
            envConstraints[julian]+=sum;
            
        }
        std::transform(envConstraints.begin(), envConstraints.end(), envConstraints.begin(),[=](double i) { return i/sim_years; });
        envConstraints.erase (envConstraints.begin());
        return envConstraints;


}

vector<double> maxRenewablesGeneration()
{

        vector<double> renewables(366,0);
        for(int t=0;t<sim_days;t++)
        {
        
            int julian=(int)julians[t];
            renewables[julian]+=utils::columnSumRow(generation_fcrps,t);
        }
        std::transform(renewables.begin(), renewables.end(), renewables.begin(),[=](double i) { return i/sim_years; });
        renewables.erase (renewables.begin());
        return renewables;

}

vector<double> floodLevel(unsigned int rank_no,unsigned int master_no,vector<double> discharge_bon)
{

         vector<double> flood_level(366,0);
         //vector<double> discharge_bon=discharge[42];
         vector<double> waterlm_coeff=utils::loadVectorFromDataFile("../data/flood_level/water_lunar.txt");
         vector<double> arima=utils::loadVectorFromDataFile("../data/flood_level/arima.txt");
         vector<double> Y_t(sim_days,0);
         vector<double> e_t(sim_days,0);
         vector<double> w_t(sim_days,0);
         double eps_minus2=0.0;
         double eps_minus1=0.0;
         std::default_random_engine generator;
         std::normal_distribution<double> distribution(0.0,0.09390710208);
         for( int t=0;t<sim_days;t++)
         {
           int julian = (int)julians[t];
           w_t[t]=distribution(generator);
           if(t==0)
                 e_t[0]=arima[0]*eps_minus1+arima[1]*eps_minus2+w_t[t];
           else if(t==1)
                 e_t[t]=arima[0]*e_t[t-1]+arima[1]*eps_minus1+w_t[t];
           else
                 e_t[t]=arima[0]*e_t[t-1]+arima[1]*e_t[t-2]+w_t[t];
    
           Y_t[t] =waterlm_coeff[0] + waterlm_coeff[1]*log(discharge_bon[t+1]/1000) + waterlm_coeff[2]*cos(2*pi*(t+1)/365)+waterlm_coeff[3]*sin(2*pi*(t+1)/365)+waterlm_coeff[4]*sin(2*pi*(t+1)/14.75)+waterlm_coeff[5]*cos(2*pi*(t+1)/14.75)+waterlm_coeff[6]*(7665)+e_t[t];
           //cout<<"errors"<<endl;
           //if(t<100)
           //cout<<w_t[t]<<endl;
           /*
           if(t>3550 and t<3600)
           {
                   cout<<"intercept "<<waterlm_coeff[0]<<endl;
                   cout<<"discharges "<<waterlm_coeff[1]*log(discharge_bon[t]/1000)<<endl;
                   cout<<"365 cos "<<waterlm_coeff[2]*cos(2*pi*t/365)<<endl;
                   cout<<"365 sin "<<waterlm_coeff[3]*sin(2*pi*t/365)<<endl;
                   cout<<"14.7 sin "<<waterlm_coeff[4]*sin(2*pi*t/14.75)<<endl;
                   cout<<"14.7 cos "<<waterlm_coeff[5]*cos(2*pi*t/14.75)<<endl;
                   cout<<"time "<<waterlm_coeff[6]*(7665)<<endl;
                   cout<<"error "<<e_t[t]<<endl;
                   cout<<"predicted "<<Y_t[t]<<endl;
           
           
           
           }
           */
           flood_level[julian]=max(flood_level[julian],exp(Y_t[t]));
          }
          
          std::transform(Y_t.begin(), Y_t.end(), Y_t.begin(),[=](double i) { return exp(i); });
          
          if(write_output)
          {
          vector<vector<double> >Yt_matrix(1,vector<double> (sim_days,0));
          Yt_matrix[0]=Y_t;
          utils::writeMatrix(Yt_matrix,"../output_opt/policy"+to_string(rank_no)+"/flood_level_sim"+to_string(master_no)+".txt", sim_days, 1);
          }
          //std::transform(flood_level.begin(), flood_level.end(), flood_level.begin(),[=](double i) { return i/sim_years; });
          flood_level.erase (flood_level.begin());
          return flood_level;

}

vector<double> getCaisoPrices(unsigned int rank_no,unsigned int master_no)
{

          int j;
          
          vector<double> caiso_prices_coeff=utils::loadVectorFromDataFile("../data/cons_opt/models/ca_linear_coef.txt");
          vector<double> arima_caiso_coeff=utils::loadVectorFromDataFile("../data/cons_opt/models/ca_arima_coef.txt");
          
          vector<vector<double> > residuals_matrix=utils::readMatrixFromDataFile("../data/cons_opt/models/CA_cons_residuals_cpp.txt");
          vector<double> quantiles_list=utils::loadVectorFromDataFile("../data/cons_opt/models/ca_quantiles_list.txt");
          vector<double> Y_t(sim_days,0);
          vector<double> e_t(sim_days+2,0);
          vector<double> w_t(sim_days+2,0);
          vector<int> bucket_index(sim_days,-1);
          w_t[0]=0;//taken the last two values from ARIMA(1,0,2) model
          w_t[1]=0;
          e_t[0]=0;//taken the last value from the predictor's model
          e_t[1]=0;
          double fitted_value;
          
          for( int t=0;t<sim_days;t++)
         {
           //w_t[t+2]=distribution(generator);
           e_t[t+2]=arima_caiso_coeff[0]*e_t[(t+2)-1]+arima_caiso_coeff[1]*e_t[(t+2)-2]+arima_caiso_coeff[2]*w_t[(t+2)-1]+arima_caiso_coeff[3]*w_t[(t+2)-2];
           Y_t[t] =  caiso_prices_coeff[0] 
                   + caiso_prices_coeff[1]*modelled_demands[0][t]
                   + caiso_prices_coeff[2]*modelled_demands[4][t]
                   + caiso_prices_coeff[3]*modelled_demands[6][t]
                   + caiso_prices_coeff[4]*modelled_demands[8][t]
                   + caiso_prices_coeff[5]*modelled_demands[10][t]
                   + caiso_prices_coeff[6]*modelled_demands[11][t]
                   + caiso_prices_coeff[7]*modelled_demands[13][t]
                   + caiso_prices_coeff[8]*modelled_demands[14][t]
                   + caiso_prices_coeff[9]*modelled_solar[0][t]
                   + caiso_prices_coeff[10]*modelled_wind[0][t]
                   + caiso_prices_coeff[11]*modelled_cahydro[0][t]
                   + caiso_prices_coeff[12]*(utils::columnSumRow(generation_fcrps,t)+extra_generation[0][t])
                   + caiso_prices_coeff[13]*sin( 2*pi*julians[t]/7)
                   + caiso_prices_coeff[14]*cos( 2*pi*julians[t]/7);

            fitted_value= Y_t[t]+e_t[t+2];
            
                    for (j=0;j<quantiles_list.size();j++)
                      {
                        if(fitted_value<quantiles_list[j]){
                          break;
                        }
                        
                      }
            //cout<<"fitted_value "<<fitted_value<<" rank_no "<<rank_no<<" master_no "<<master_no<<endl;
            
            if(j>=20)
                j=19;
            bucket_index[t]=j; 
                   
            w_t[t+2]=residuals_matrix[j][rand() % residuals_matrix[j].size()];
            e_t[t+2]=e_t[t+2]+w_t[t+2];      
            Y_t[t]=Y_t[t]+e_t[t+2];
           
           //Y_t[t]=max(10.0,Y_t[t]);
        }
        if(write_output)
        {
            vector<vector<double> >Yt_matrix(1,vector<double> (sim_days,0));
            Yt_matrix[0]=Y_t;
            //utils::writeMatrix(Yt_matrix,"../output/Caiso_prices0.txt", sim_days, 1);
            utils::writeMatrix(Yt_matrix,"../output_opt/policy"+to_string(rank_no)+"/Caiso_prices"+to_string(master_no)+".txt", sim_days, 1);
        }
        return Y_t;


}

vector<double> getMidCPrices(unsigned int rank_no,unsigned int master_no)
{

          int j;
          vector<double> caiso_prices_coeff=utils::loadVectorFromDataFile("../data/cons_opt/models/midc_linear_coef.txt");
          vector<double> arima_caiso_coeff=utils::loadVectorFromDataFile("../data/cons_opt/models/midc_arima_coef.txt");
          
          vector<vector<double> > residuals_matrix=utils::readMatrixFromDataFile("../data/cons_opt/models/MidC_cons_residuals_cpp.txt");
          vector<double> quantiles_list=utils::loadVectorFromDataFile("../data/cons_opt/models/midc_quantiles_list.txt");
          vector<double> Y_t(sim_days,0);
          vector<double> e_t(sim_days+2,0);
          vector<double> w_t(sim_days+2,0);
          vector<int> bucket_index(sim_days,-1);
          w_t[0]=0;//taken the last two values from ARIMA(1,0,2) model
          w_t[1]=0;
          e_t[0]=0;//taken the last value from the predictor's model
          e_t[1]=0;
          double fitted_value;
          
          for( int t=0;t<sim_days;t++)
         {
           //w_t[t+2]=distribution(generator);
           e_t[t+2]=arima_caiso_coeff[0]*e_t[(t+2)-1]+arima_caiso_coeff[1]*e_t[(t+2)-2]+arima_caiso_coeff[2]*w_t[(t+2)-1]+arima_caiso_coeff[3]*w_t[(t+2)-2];
           Y_t[t] =  caiso_prices_coeff[0] 
                   + caiso_prices_coeff[1]*modelled_demands[0][t]
                   + caiso_prices_coeff[2]*modelled_demands[4][t]
                   + caiso_prices_coeff[3]*modelled_demands[6][t]
                   + caiso_prices_coeff[4]*modelled_demands[7][t]
                   + caiso_prices_coeff[5]*modelled_demands[8][t]
                   + caiso_prices_coeff[6]*modelled_demands[11][t]
                   + caiso_prices_coeff[7]*modelled_demands[13][t]
                   + caiso_prices_coeff[8]*modelled_demands[14][t]
                   + caiso_prices_coeff[9]*modelled_solar[0][t]
                   + caiso_prices_coeff[10]*modelled_wind[0][t]
                   + caiso_prices_coeff[11]*modelled_cahydro[0][t]
                   + caiso_prices_coeff[12]*(utils::columnSumRow(generation_fcrps,t)+extra_generation[0][t])
                   + caiso_prices_coeff[13]*sin( 2*pi*julians[t]/7)
                   + caiso_prices_coeff[14]*cos( 2*pi*julians[t]/7);

            fitted_value= Y_t[t]+e_t[t+2];
                    for (j=0;j<quantiles_list.size();j++)
                      {
                        if(fitted_value<quantiles_list[j]){
                          break;
                        }
                        
                      }
            if(j>=20)
                j=19;
            bucket_index[t]=j;        
            w_t[t+2]=residuals_matrix[j][rand() % residuals_matrix[j].size()];
            e_t[t+2]=e_t[t+2]+w_t[t+2];      
            Y_t[t]=Y_t[t]+e_t[t+2];
           
           Y_t[t]=max(10.0,exp(Y_t[t]));
        }
        if(write_output)
        {
            vector<vector<double> >Yt_matrix(1,vector<double> (sim_days,0));
            Yt_matrix[0]=Y_t;
            //utils::writeMatrix(Yt_matrix,"../output/MidC_prices0.txt", sim_days, 1);
            utils::writeMatrix(Yt_matrix,"../output_opt/policy"+to_string(rank_no)+"/MidC_prices"+to_string(master_no)+".txt", sim_days, 1);
        }
        
        return Y_t;
}


vector<double> getRevenueBPA(unsigned int rank_no,unsigned int master_no)
{
   
   
   vector<vector<double> > MidC(1,vector<double> (sim_days,0));
   vector<vector<double> > CAISO(1,vector<double> (sim_days,0));
   MidC[0]=getMidCPrices(rank_no,master_no);
   CAISO[0]=getCaisoPrices(rank_no,master_no);
   vector<double> midc_prices = {MidC[0].begin() + 122, MidC[0].end() - 243}; 
   vector<double> caiso_prices = {CAISO[0].begin() + 122, CAISO[0].end() - 243}; 

   
   vector<double> load_vector=utils::loadVectorFromDataFile("../data/Revenue Model/load_2018.txt");
   double custom_redux=0;
   double PF_load_y=load_vector[13]-custom_redux*load_vector[13];
   double IP_load_y=load_vector[3]-custom_redux*load_vector[3];
   double ET_load_y=load_vector[14];
   
   
   
   vector<vector<double> > extra_bpa_dams=utils::loadMatrix("../data/optimize_four/extra_bpa.txt",bpa_rev_sim_days,1,365,start_col);
   utils::transpose(extra_bpa_dams);
   //vector<double> extra_gen=utils::columnSumMatrix(extra_bpa_dams);
   vector<double> bpa_gen=utils::columnSumMatrix(generation_shoal);
   vector<double> bpa_gen_sub = {bpa_gen.begin() + 122, bpa_gen.end() - 243}; 
   vector<double> BPA_hydro(bpa_gen_sub.size(),0);
   std::transform(bpa_gen_sub.begin(),bpa_gen_sub.end(), extra_bpa_dams[0].begin(),BPA_hydro.begin(),std::plus<double>());
   std::transform(BPA_hydro.begin(), BPA_hydro.end(), BPA_hydro.begin(),[=](double i) { return i/24; });
   
   
   //cout<<" in revenue model 1"<<endl;
   vector<double> net_resources=utils::loadVectorFromDataFile("../data/Revenue Model/BPA_net_resources_2018.txt");
   double Nuc_y=net_resources[7];
   double Wind_y=net_resources[8];
   double Purch_y=net_resources[10];
   
   vector<double> PF_rates=utils::loadVectorFromDataFile("../data/Revenue Model/PF_rates_2018.txt");
   vector<double> IP_rates=utils::loadVectorFromDataFile("../data/Revenue Model/IP_rates_2018.txt");
   /*
   vector<vector<double> > BPAT_load=utils::loadMatrix("../data/cons_opt/Datafiles/daily_load_bpa.txt",sim_days,1,start_row,start_col);
   BPAT_load= utils::matrixReshape( BPAT_load,sim_days,1);
   utils::transpose(BPAT_load);
   */
   vector<double> load_bpa=utils::loadVectorFromDataFile("../data/optimize_four/daily_load_bpa.txt");
   vector<double> BPAT_load={load_bpa.begin() + 365, load_bpa.begin() + 365+ bpa_rev_sim_days};
   cout<<"BPAT_load.size() "<<BPAT_load.size()<<endl;
   //utils::transpose(modelled_wind);
   vector<double> wind_bpa=utils::loadVectorFromDataFile("../data/optimize_four/wind_surrogate_bpa.txt");
   vector<double> BPAT_wind={wind_bpa.begin() + 365, wind_bpa.begin() +365 +bpa_rev_sim_days}; 
   cout<<"BPAT_wind.size() "<<BPAT_wind.size()<<endl;
   //std::transform(modelled_wind[0].begin(), modelled_wind[0].end(), BPAT_wind.begin(),[=](double i) { return i *(0.766/1.766); });
   

   double costs_y=utils::loadVectorFromDataFile("../data/Revenue Model/BPA_yearly_costs.txt")[8];
   
   vector<double> load_ratio(BPAT_load.size(),0);
   double BPAT_load_mean=meanVector(BPAT_load);
   //cout<<BPAT_load_mean<<endl;
   std::transform(BPAT_load.begin(), BPAT_load.end(), load_ratio.begin(),[=](double i) { return i /BPAT_load_mean; });
   
   vector<double> wind_ratio(BPAT_wind.size(),0);
   double BPAT_wind_mean=meanVector(BPAT_wind);
   std::transform(BPAT_wind.begin(), BPAT_wind.end(), wind_ratio.begin(),[=](double i) { return i /BPAT_wind_mean; });
   
   vector<vector<double> >PF_load(1,vector<double> (BPAT_wind.size(),0));
   vector<vector<double> >IP_load(1,vector<double> (BPAT_wind.size(),0));
   std::transform(load_ratio.begin(), load_ratio.end(), PF_load[0].begin(),[=](double i) { return i *PF_load_y; });
   std::transform(load_ratio.begin(), load_ratio.end(), IP_load[0].begin(),[=](double i) { return i *IP_load_y; });
   //cout<<PF_load.size()<<endl;
   //cout<<PF_load[0].size()<<endl;
   //cout<<"load_ratio is"<<load_ratio.size()<<endl;
   //utils::transpose(PF_load);
   PF_load=utils::matrixReshape( PF_load,365,bpa_rev_sim_years);
   IP_load=utils::matrixReshape( IP_load,365,bpa_rev_sim_years);
   //cout<<PF_load.size()<<endl;
   //cout<<PF_load[0].size()<<endl;
   utils::transpose(PF_load);
   utils::transpose(IP_load);
   
   
   double PF_load_avg=meanVector(rowSumsMatrix(PF_load));
   double IP_load_avg=meanVector(rowSumsMatrix(IP_load));
   //cout<<PF_load_avg<<endl;
   //cout<<IP_load_avg<<endl;
   
   
   vector<double> ET_load(load_ratio.size(),0);
   vector<double> Purch(load_ratio.size(),0);
   vector<double> Wind(load_ratio.size(),0);
   vector<double> Nuc(load_ratio.size(),1);
   std::transform(load_ratio.begin(), load_ratio.end(), ET_load.begin(),[=](double i) { return i *ET_load_y; });
   std::transform(load_ratio.begin(), load_ratio.end(), Purch.begin(),[=](double i) { return i *Purch_y; });
   std::transform(wind_ratio.begin(), wind_ratio.end(), Wind.begin(),[=](double i) { return i *Wind_y; });
   std::transform(Nuc.begin(), Nuc.end(), Nuc.begin(),[=](double i) { return i *Nuc_y; });

    /*
    vector<vector<double> > MidC=utils::loadMatrix("../output/MidC_prices.txt",sim_days,1,0,start_col);
    utils::transpose(MidC);
    
    vector<vector<double> > CAISO=utils::loadMatrix("../output/Caiso_prices.txt",sim_days,1,0,start_col);
    utils::transpose(CAISO);
    */
    
    double ExR=0.71;   
    double TA=1000;
    vector<double> temp(load_ratio.size(),0);
    vector<double> temp1(load_ratio.size(),0);
    vector<double> temp2(load_ratio.size(),0);
    vector<double> trans_losses(load_ratio.size(),0);
    vector<double> BPA_res(load_ratio.size(),0);
    //trans_losses=3*(Wind + BPA_hydro + Nuc)/100;
    std::transform(Wind.begin(), Wind.end(), Nuc.begin(),temp.begin(),std::plus<double>());
    std::transform(BPA_hydro.begin(),BPA_hydro.end(), temp.begin(),temp1.begin(),std::plus<double>());
    std::transform(temp1.begin(), temp1.end(), trans_losses.begin(),[=](double i) { return i *0.03; });
    std::transform(Purch.begin(),Purch.end(), temp1.begin(),temp2.begin(),std::plus<double>());
    std::transform(temp2.begin(),temp2.end(), trans_losses.begin(),BPA_res.begin(),std::minus<double>());
    //BPA_res=pd.DataFrame(data=(Wind + BPA_hydro + Purch + Nuc)-trans_losses)
    
    
    utils::transpose(PF_load);
    utils::transpose(IP_load);
    PF_load=utils::matrixReshape( PF_load,365*bpa_rev_sim_years,1);
    IP_load=utils::matrixReshape( IP_load,365*bpa_rev_sim_years,1);
    utils::transpose(PF_load);
    utils::transpose(IP_load);
    
    vector<double> SD(load_ratio.size(),0);
    std::transform(PF_load[0].begin(), PF_load[0].end(), IP_load[0].begin(),temp.begin(),std::plus<double>());
    std::transform(ET_load.begin(), ET_load.end(), temp.begin(),temp1.begin(),std::plus<double>());
    std::transform(BPA_res.begin(), BPA_res.end(), temp1.begin(),SD.begin(),std::minus<double>());
    //for(int i=0;i<1400;i++)
    //cout<<SD[i]<<endl;
    
    vector<vector<double> > months_matrix=utils::loadMatrix("../data/Revenue Model/months.txt",bpa_rev_sim_days,1,365,start_col);
    utils::transpose(months_matrix);
    vector<double> months_rev=months_matrix[0];
    makePF_rates_monthsIndexing();
    double CRAC=0;
      int num_years=bpa_rev_sim_years+1;
      double RatePF,RateIP,Ex,losses,Net_res1; 
      double start_res=158.7*pow(10,6);
      double starting_BA = 2.421*pow(10,9);
      double Treas_fac1=320*pow(10,6);
      double Treas_fac2=430*pow(10,6);
      double trans_BA= 9.782*pow(10,6)*0.4;
      double Used_TF= 0; 
      vector<double> PF_rev(load_ratio.size(),0);
      vector<double> IP_rev(load_ratio.size(),0);
      vector<double> P(load_ratio.size(),0);
      vector<double> SS(load_ratio.size(),0);
      vector<double> BPA_rev_d(load_ratio.size(),0);
      vector<double> crac(load_ratio.size(),0);
      vector<double> BPA_Net_rev_y(num_years,0);
      vector<double> Reserves(num_years+1,0);
      vector<double> Remaining_BA(num_years+1,0);
      vector<double> TF1(num_years+1,0);
      vector<double> TF2(num_years+1,0);
      vector<double> CRAC_y(num_years+1,0);
      vector<double> CRAC_rev(num_years,0);
      vector<double> TTP(num_years,-1);
      Reserves[1]=start_res;
      Remaining_BA[1]=starting_BA;
      
      double p=10; 
      double p2=32;
      double d=1;
      double e=1;
      
      sdNegativeSumOpt=0;
      for(int i=0;i<SD.size();i++)
      {
              
           RatePF=PF_rates[PF_rates_months[months_rev[i]]];
           RatePF+=CRAC;
           PF_rev[i]=PF_load[0][i]*RatePF*24;
           RateIP=IP_rates[PF_rates_months[months_rev[i]]];
           RateIP+=CRAC;
           IP_rev[i]=IP_load[0][i]*RateIP*24;
           //cout<<"CRAC is"<<CRAC<<endl;
           crac[i]=CRAC;
           
           if(SD[i]<0)
           {
                  if(caiso_prices[i]>midc_prices[i])
                      P[i]=SD[i]*midc_prices[i]*24;
                  else
                      P[i]=SD[i]*caiso_prices[i]*24;
                  SS[i]=0;
                  sdNegativeSumOpt+=SD[i];
           }
           else
           {
                  P[i]=0;
                  if(caiso_prices[i]>midc_prices[i])
                      Ex=min(SD[i],TA);
                  else
                      Ex=0;
                  SS[i]=ExR*(Ex*caiso_prices[i]*24)+(SD[i]-Ex)*midc_prices[i]*24;
           
           }
           BPA_rev_d[i]= PF_rev[i] + IP_rev[i] + SS[i] + P[i];
           
            if(((i+1)%365)==0)
            {
                  int year=(i+1)/365;
                  //cout<<"year is"<<year<<endl;
                  bool bol=(year%2 == 0);
                  double PF_load_i = accumulate(PF_load[0].begin()+(year-1)*365,PF_load[0].begin()+(year)*365,0.0);     
                  double IP_load_i = accumulate(IP_load[0].begin()+(year-1)*365,IP_load[0].begin()+(year)*365,0.0);
                  double tot_load_i = PF_load_i + IP_load_i;
                  BPA_Net_rev_y[year]=accumulate(BPA_rev_d.begin()+(i-364),BPA_rev_d.begin()+i+1,0.0)- costs_y;
                  if(BPA_Net_rev_y[year]<0.0)
                  {
                             losses=-BPA_Net_rev_y[year];
                             Net_res1= Reserves[year] - losses;
                             Reserves[year+1] = max(Net_res1, 0.0);
                             
                             if(Net_res1 < 0.0)
                             {
                                   losses=-Net_res1;
                                   if((Remaining_BA[year] - Used_TF) > 750*pow(10,6))
                                    {
                                    TF1[year+1]=min(losses , Treas_fac1-TF1[year]*bol);
                                    Used_TF+=TF1[year+1];
                                          if ((Treas_fac1-TF1[year]*bol - losses)<0)
                                          {
                                                 losses= - (Treas_fac1-TF1[year]*bol - losses);
                                                 //cout<<"first crac"<<endl;
                                                 //cout<<"losses is"<<losses<<endl;
                                                 //cout<<"tot load"<<tot_load_i<<endl;
                                                 CRAC+=calculate_CRAC(losses, tot_load_i);
                                                 //set max crac as +5$/MWh per ensemble
                                                 CRAC=min(CRAC, 5.0);
                                                 TF2[year+1]=min(losses , Treas_fac2-TF2[year]*bol);
                                                 Used_TF+=TF2[year+1];
                                                 if((Treas_fac2-TF2[year]*bol - losses) <= 0)
                                                    TTP[year]=-4;
                                          }
                                    }
                                    else
                                   {
                                         CRAC+=calculate_CRAC(losses, tot_load_i);
                                         CRAC=min(CRAC, 5.0);
                                         TTP[year]=losses;
                                   }
                            }
                            
                  }
                  else
                  {
                            
                            Reserves[year+1]= min(Reserves[year] + 0.01*p*BPA_Net_rev_y[year],608691000-Treas_fac1);
                            //cout<<"Debt optimization, added:"<<endl;
                            //cout<<0.01*p2*BPA_Net_rev_y[year]<<endl;
                  
                  }
                            
                  CRAC_y[year+1]=CRAC;
                  Remaining_BA[year+1]= max(0.0, max(Remaining_BA[year] - 484*pow(10,6) + trans_BA, Remaining_BA[year] - 484*pow(10,6) + trans_BA + 0.01*p2*BPA_Net_rev_y[year]));
                  if(year%20==0)
                  {
                          Reserves[year+1]=start_res;
                          CRAC=0;
                          CRAC_y[year+1]=0;
                          Used_TF=0;
                          Remaining_BA[year+1] = starting_BA;
                  }
                  
                   
            }      
              //cout<<"i is "<<i<<endl;   
           
      }
        if(write_output)
        {
            vector<vector<double> >Yt_matrix(1,vector<double> (num_years,0));
            Yt_matrix[0]=BPA_Net_rev_y;
            utils::writeMatrix(Yt_matrix,"../output_opt/policy"+to_string(rank_no)+"/BPA_Net_Revenue"+to_string(master_no)+".txt", num_years, 1);
            
            
            vector<vector<double> >Sd_matrix(1,vector<double> (load_ratio.size(),0));
            Sd_matrix[0]=SD;
            utils::writeMatrix(Sd_matrix,"../output_opt/policy"+to_string(rank_no)+"/SD"+to_string(master_no)+".txt", num_years, 1);
        }     
        
        BPA_Net_rev_y.erase (BPA_Net_rev_y.begin());
        return BPA_Net_rev_y;
}

vector<double> getObjectives(unsigned int rank_no,unsigned int master_no)
{
        //cout<<"before simulate"<<rank_no<<" "<<master_no<<endl;
        simulate(rank_no,master_no);
        //cout<<"after getobjectives"<<rank_no<<" "<<master_no<<endl;
        if(write_output)
        writeSimulationOutput(rank_no,master_no);
        vector<double> J;
        J.push_back(meanVector(minEnvironmentConstraints()));
        J.push_back(-1*meanVector(maxRenewablesGeneration()));
        vector<double> floodlevel=floodLevel(rank_no,master_no,discharge[42]);
        //vector<double> floodlevel_historical=floodLevel(rank_no,master_no,modelled_discharge[42]);
        //cout<<"historical flood "<<*max_element(floodlevel_historical.begin(), floodlevel_historical.end())<<endl;
        J.push_back(*max_element(floodlevel.begin(), floodlevel.end()));
        J.push_back((meanVector(getRevenueBPA(rank_no,master_no))/1000000)*-1);

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
                makeHYSSRIndexing();
                makeinverseHYSSRIndexing();
                initialize();
                makeFlowRoute();
                //cout<<"INITIALIZE RBFS "<<to_string(rank_no)<<" "<<to_string(master_no)<<endl;
                }
                
                vector<vector<double> > min_max_inflow_rbfs_range=utils::readMatrixFromDataFile("../data/optimize_four/inflows_rbfs_ranges.txt");

                
                //cout<<to_string(rank_no)+"inflow range size "<<to_string(master_no)<< "min_max_inflow_rbfs_range[0]"<< min_max_inflow_rbfs_range[5][1]<<endl;
                vector<vector<double> > min_max_releases_range=utils::readMatrixFromDataFile("../data/optimize_four/discharge_rbfs_ranges.txt");
                
                for(int i=0;i<1;i++)
                {   
                p_param[i].tPolicy=1;
                p_param[i].policyInput=9;
                p_param[i].policyOutput=4;
                p_param[i].policyStr=14;
                param_function* mp1 = new rbf(p_param[i].policyInput,p_param[i].policyOutput,p_param[i].policyStr);
                //param_function* mp1 = new ncRBF(p_param[i].policyInput,p_param[i].policyOutput,p_param[i].policyStr);
                mPolicy[i] = mp1;
                }
                for(int i=0;i<1;i++)
                {   
                
                        
                        p_param[i].mIn.push_back(0);
                        p_param[i].mIn.push_back(0);
                        p_param[i].mIn.push_back(0);
                        p_param[i].mIn.push_back(0);
                        p_param[i].mIn.push_back(min_max_inflow_rbfs_range[0][0]);
                        p_param[i].mIn.push_back(min_max_inflow_rbfs_range[1][0]);
                        p_param[i].mIn.push_back(min_max_inflow_rbfs_range[2][0]);
                        p_param[i].mIn.push_back(min_max_inflow_rbfs_range[3][0]);
                        p_param[i].mIn.push_back(-1);
                        
                        
                        p_param[i].MIn.push_back(s_max[2]);
                        p_param[i].MIn.push_back(s_max[5]);
                        p_param[i].MIn.push_back(s_max[9]);
                        p_param[i].MIn.push_back(s_max[12]);
                        p_param[i].MIn.push_back(min_max_inflow_rbfs_range[0][1]);
                        p_param[i].MIn.push_back(min_max_inflow_rbfs_range[1][1]);
                        p_param[i].MIn.push_back(min_max_inflow_rbfs_range[2][1]);
                        p_param[i].MIn.push_back(min_max_inflow_rbfs_range[3][1]);
                        p_param[i].MIn.push_back(1);
                        

                        p_param[i].mOut.push_back(0);
                        p_param[i].mOut.push_back(0);
                        p_param[i].mOut.push_back(0);
                        p_param[i].mOut.push_back(0);
                        p_param[i].MOut.push_back(min_max_releases_range[0][1]);
                        p_param[i].MOut.push_back(min_max_releases_range[1][1]);
                        p_param[i].MOut.push_back(min_max_releases_range[2][1]);
                        p_param[i].MOut.push_back(min_max_releases_range[3][1]);
                        
                        mPolicy[i]->setMinInput(p_param[i].mIn); 
                        mPolicy[i]->setMinOutput(p_param[i].mOut);
                        mPolicy[i]->setMaxInput(p_param[i].MIn); 
                        mPolicy[i]->setMaxOutput(p_param[i].MOut);
                        
                }
                //cout<<"after initialize"<<endl;
        }


void evaluate(double* var, double* obj, double* constrs, unsigned int inp_rank_no,unsigned int inp_master_no){

        unsigned int rank_no= inp_rank_no,master_no=inp_master_no;
        int decsvars=281;
        //cout<<"evaluate first line"<<rank_no<<" "<<master_no<<endl;
        initializeRBFs(rank_no,master_no);
        //cout<<"evaluate initializeRBFs"<<rank_no<<" "<<master_no<<endl;
        //cout<<"inside evaluate"<<endl;
        unsigned int p2 = 0;
        int getArrayLength;
        for(int k=0;k<1;k++)
                {
                
                double fullVar[p_param[k].policyStr*(2*p_param[k].policyInput+p_param[k].policyOutput)];
                //double fullVar[p_param[k].policyStr*(2*p_param[k].policyInput+p_param[k].policyOutput)];
                unsigned int p1 = 0;
                
                //ncRBF constant
                /*
                if (var[304]<0.5)
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
                }*/
                
                /*
                for(unsigned int l = 0; l < p_param[k].policyOutput; l++)
                        {
                                fullVar[p1++] = 0.0;
                                p2++;
                        }
                */
                


                for(unsigned int i = 0; i < p_param[k].policyStr; i++){
                
                        for(unsigned int j = 0; j < p_param[k].policyInput-1; j++)
                        {
                                fullVar[p1++] = var[p2++];
                                fullVar[p1++] = var[p2++];
                        }
                        for(unsigned int j = 0; j < 1 ; j++)
                        {
                                // set centers to 0 and radii to 1 for sin and cos inputs
                                fullVar[p1++] = 0.0;
                                fullVar[p1++] = 1.0;
                        }
                        for(unsigned int l = 0; l < p_param[k].policyOutput; l++)
                        {
                                fullVar[p1++] = var[p2++];
                        }
                        
                }
                phi=p2++;

                mPolicy[k]->setParameters(fullVar);
                
                //cout<<"mPolicy[k]->param "<< mPolicy[k]->param << endl;

                
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

                }
                
                //cout<<"before getobjectives"<<rank_no<<" "<<master_no<<endl;
                vector<double> J ;
                J = getObjectives(rank_no,master_no);
                //cout<<"after getobjectives"<<rank_no<<" "<<master_no<<endl;
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

                //mPolicy[0]->getParameters();
                
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