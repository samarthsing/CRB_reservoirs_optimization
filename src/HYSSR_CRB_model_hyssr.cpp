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
//const int M_H=12;
const int M_HR=23;
int sim_years=99;
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
int extra_start=0;
double sdNegativeSumOpt;

vector<vector<double> > modelled_discharge;
//vector<vector<double> > modelled_powerflow;
vector<vector<double> > modelled_storage;
vector<vector<double> > modelled_generation;
vector<vector<double> > inflow;
vector<vector<double> > synth_weather;
vector<vector<double> > discharge;
vector<vector<double> > storage;
//Data structure to model the prices
vector<vector<double> > modelled_demands;
vector<vector<double> > extra_generation;
vector<vector<double> > ca_surrogate;
vector<vector<double> > pnw_surrogate;
vector<vector<double> > ca_surrogate_df;
vector<vector<double> > pnw_surrogate_df;

void readDataFiles(unsigned int rank_no,unsigned int master_no)
{
int folder_no=master_no;
//cout<<"reading files begin"<<endl;
modelled_discharge=utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/discharge.txt",sim_days+1,M,365-1+extra_start,start_col);//water year
//modelled_powerflow=utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/powerflow.txt",sim_days+1,M,365-1+extra_start,start_col);//water year
modelled_storage=utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/storage.txt",sim_days+1,M,365-1+extra_start,start_col);//water year
modelled_generation=utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/generation.txt",sim_days+1,M+1,365-1+extra_start,start_col);//water year
inflow=utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/inflows.txt",sim_days+1,M_comp,start_row+365-1+extra_start,start_col);//jan 1
synth_weather=utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/synth_weather.txt",sim_days+1,34,start_row+365-1+extra_start,start_col);//jan 1

//Data structure to model the prices
modelled_demands    = utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/daily_load_data.txt",sim_days,6,start_row+extra_start,start_col);//jan 1
pnw_surrogate_df       = utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/pnw_surrogate_data_reduced_fossils.txt",sim_days,10,start_row+extra_start,start_col);//jan 1
ca_surrogate_df        = utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/ca_surrogate_data.txt",sim_days,10,start_row+extra_start,start_col);//jan 1
extra_generation    = utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/extra_gen.txt",sim_days,1,start_row+extra_start,start_col);//jan 1

discharge=modelled_discharge;
storage=modelled_storage;
pnw_surrogate=pnw_surrogate_df;
ca_surrogate=ca_surrogate_df;
//cout<<ca_surrogate.size()<<endl;
//cout<<ca_surrogate[0].size()<<endl;

//cout<<"reading files end"<<endl;
}


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
//vector<vector<double> > discharge=modelled_discharge;
//vector<vector<double> > storage(M,vector<double> (sim_days+1,0));

/*
//Data structure to model the prices
//vector<vector<double> > modelled_demands    = utils::loadMatrix("../data/optimize_four/demands_surrogate.txt",sim_days,15,start_row,start_col);//jan 1
vector<vector<double> > modelled_demands    = utils::loadMatrix("../data/optimize_four_1000/"+to_string(master_no)+"/daily_load_data.txt",sim_days,6,start_row+extra_start,start_col);//jan 1
vector<vector<double> > modelled_solar      = utils::loadMatrix("../data/optimize_four_1000/"+to_string(master_no)+"/solar_surrogate.txt",sim_days,1,start_row+extra_start,start_col);//jan 1
vector<vector<double> > modelled_wind       = utils::loadMatrix("../data/optimize_four_1000/"+to_string(master_no)+"/wind_surrogate.txt",sim_days,1,start_row+extra_start,start_col);//jan 1
vector<vector<double> > modelled_cahydro    = utils::loadMatrix("../data/optimize_four_1000/"+to_string(master_no)+"/CA_hydro_surrogate.txt",sim_days,1,start_row+extra_start,start_col);//jan 1
vector<vector<double> > extra_generation    = utils::loadMatrix("../data/optimize_four_1000/"+to_string(master_no)+"/extra_gen.txt",sim_days,1,start_row+extra_start,start_col);//jan 1
*/
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


//vector<vector<double> > storage=modelled_storage;
vector<vector<double> > total_inflow=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > generation=vector<vector<double> >(M+1,vector<double> (sim_days+1,0));
vector<vector<double> > spill=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > powerflow=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > tailwater=vector<vector<double> >(M,vector<double> (sim_days+1,0));
vector<vector<double> > stage_hyssr=vector<vector<double> >(M_H,vector<double> (sim_days+1,0));
vector<vector<double> > tailwater_hyssr=vector<vector<double> >(M_H,vector<double> (sim_days+1,0));
vector<vector<double> > heads=vector<vector<double> >(M_H,vector<double> (sim_days+1,0));


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
//set<int> shoal_reservoirs={2,7,9,12,35,36,37,38,39,40,41,42};
//set<int> hyssr_reservoirs={2,5,7,9,12,28,29,30,31,32,33,35,36,37,38,39,40,41,42,46,4,17,6,18,13,19,20,21,16,22};
set<int> hyssr_reservoirs={2,5,7,9,12,28,29,30,31,32,33,35,36,37,38,39,40,41,42,13,19,20,21};
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


vector<vector<double> > min_spills=utils::readMatrixFromDataFile("../data/HYSSR/min_constraints/min_spills_opt.txt");
vector<vector<double> > max_spills=utils::readMatrixFromDataFile("../data/HYSSR/min_constraints/max_spills_opt.txt");
vector<vector<double> > min_powerh=utils::readMatrixFromDataFile("../data/HYSSR/min_constraints/min_ph_opt.txt");
vector<double> max_powerh=utils::loadVectorFromDataFile("../data/HYSSR/min_constraints/max_ph_opt.txt");
vector<double> max_discharge_alb=utils::loadVectorFromDataFile("../data/ALB/alb_max_discharge.txt");
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
     opt_index[2]=0;
     opt_index[5]=1;
     opt_index[9]=2;
     opt_index[12]=3;
     //opt_index[9]=0;
}

void makeinverseOPTIndexing()
{

      inverse_opt[0]=2;
      inverse_opt[1]=5;
      inverse_opt[2]=9;
      inverse_opt[3]=12;
      //inverse_opt[0]=9;
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
     //hyssr_index[46]=19;
     //hyssr_index[4]=20;
     //hyssr_index[17]=21;
     //hyssr_index[6]=22;
     //hyssr_index[18]=23;
     hyssr_index[13]=19;
     hyssr_index[19]=20;
     hyssr_index[20]=21;
     hyssr_index[21]=22;
     //hyssr_index[16]=28;
     //hyssr_index[22]=29;

}
/*
void makeHYSSRIndexing()
{
     hyssr_index[2]=0;
     //hyssr_index[5]=1;
     hyssr_index[7]=1;
     hyssr_index[9]=2;
     hyssr_index[12]=3;
     //hyssr_index[28]=5;
     hyssr_index[35]=4;
     hyssr_index[36]=5;
     hyssr_index[37]=6;
     hyssr_index[38]=7;
     hyssr_index[39]=8;
     hyssr_index[40]=hyssr_index[39]+1;
     hyssr_index[41]=hyssr_index[40]+1;
     hyssr_index[42]=hyssr_index[41]+1;
     hyssr_index[29]=hyssr_index[42]+1;
     hyssr_index[30]=hyssr_index[29]+1;
     hyssr_index[31]=hyssr_index[30]+1;
     hyssr_index[32]=hyssr_index[31]+1;
     hyssr_index[33]=hyssr_index[32]+1;    
     hyssr_index[46]=hyssr_index[33]+1;
     hyssr_index[4]=hyssr_index[46]+1;
     hyssr_index[17]=hyssr_index[4]+1;
     hyssr_index[6]=hyssr_index[17]+1;
     hyssr_index[18]=hyssr_index[6]+1;
     hyssr_index[13]=hyssr_index[18]+1;
     hyssr_index[19]=hyssr_index[13]+1;
     hyssr_index[20]=hyssr_index[19]+1;
     hyssr_index[21]=hyssr_index[20]+1;
     hyssr_index[16]=hyssr_index[21]+1;
     hyssr_index[22]=hyssr_index[16]+1;

}*/

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
      //inverse_hyssr[19]=46;
      //inverse_hyssr[20]=4;
      //inverse_hyssr[21]=17;
      //inverse_hyssr[22]=6;
      //inverse_hyssr[23]=18;
      inverse_hyssr[19]=13;
      inverse_hyssr[20]=19;
      inverse_hyssr[21]=20;
      inverse_hyssr[22]=21;
      //inverse_hyssr[28]=16;
      //inverse_hyssr[29]=22;
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


double getGCL_fctarget(int julian,double TDA_for)
{
         double GCL_fctarget;
         if(TDA_for <= 57000)
                 GCL_fctarget = GCL_fc[julian-1][1];
         else if(TDA_for > 57000 && TDA_for <= 60000)
                 GCL_fctarget = GCL_fc[julian-1][1] + ((TDA_for - 57000)/(60000-57000))*(GCL_fc[julian-1][2] - GCL_fc[julian-1][1]);
         else if(TDA_for > 60000 && TDA_for <= 63250)
                 GCL_fctarget = GCL_fc[julian-1][2] + ((TDA_for - 60000)/(63250-60000))*(GCL_fc[julian-1][3] - GCL_fc[julian-1][2]);
         else if(TDA_for > 63250 and TDA_for <= 65000)
                 GCL_fctarget = GCL_fc[julian-1][3] + ((TDA_for - 63250)/(65000-63250))*(GCL_fc[julian-1][4] - GCL_fc[julian-1][3]);
         else if(TDA_for > 65000 and TDA_for <= 67660)
                 GCL_fctarget = GCL_fc[julian-1][4] + ((TDA_for - 65000)/(67660-65000))*(GCL_fc[julian-1][5] - GCL_fc[julian-1][4]);
         else if(TDA_for > 67660 and TDA_for <= 71000)
                 GCL_fctarget = GCL_fc[julian-1][5] + ((TDA_for - 67660)/(71000-67660))*(GCL_fc[julian-1][6] - GCL_fc[julian-1][5]);
         else if(TDA_for > 71000 and TDA_for <= 75000)
                 GCL_fctarget = GCL_fc[julian-1][6] + ((TDA_for - 71000)/(75000-71000))*(GCL_fc[julian-1][7] - GCL_fc[julian-1][6]);
         else if(TDA_for > 75000 and TDA_for <= 87500)
                 GCL_fctarget = GCL_fc[julian-1][7] + ((TDA_for - 75000)/(87500-75000))*(GCL_fc[julian-1][8] - GCL_fc[julian-1][7]);
         else if(TDA_for > 87500 and TDA_for <= 100000)
                 GCL_fctarget = GCL_fc[julian-1][8] + ((TDA_for - 87500)/(100000-87500))*(GCL_fc[julian-1][9] - GCL_fc[julian-1][8]);
         else if(TDA_for > 100000)
                 GCL_fctarget = GCL_fc[julian-1][9];
   
         return  GCL_fctarget;             
}


double getGCL_VRCLL(int julian)
{
              double GCL_VRCLL;
              double vrcll_conversion_factor=2.4466;
              if(julian >=1 && julian <=31)
                     GCL_VRCLL=1778.9*vrcll_conversion_factor;
              else if(julian >=32 && julian<=59)
                     GCL_VRCLL=1054.5*vrcll_conversion_factor;
              else if(julian >=60 && julian <=90)
                     GCL_VRCLL=418.7*vrcll_conversion_factor;
              else if(julian >=91 && julian <=120)
                     GCL_VRCLL=418.7*vrcll_conversion_factor;
              else if(julian >=121 && julian <=151)
                     GCL_VRCLL=843.7*vrcll_conversion_factor;
              else if(julian >=152 && julian <=181)
                     GCL_VRCLL=2411.3*vrcll_conversion_factor;
              else if(julian >=182 && julian <=212)
                     GCL_VRCLL=2614.3*vrcll_conversion_factor;
                     
              return GCL_VRCLL;
}

void modelGCL(int i,int t,vector<double> GCL_sfore,vector<double> Storage_left,vector<double> GCL_VRC,vector<double> GCL_ORC,vector<double> spokane_flow)
{

          int julian = (int)julians[t-1];
          int year = (int)years[t-1];
          int floodyear;
          if(julian<182)
               floodyear = year;
          else
               floodyear = int(min(year+1,no_years));
            
          double ICF = ICFs[floodyear-1];
          double fillstart = max(ICF - 20,0.0);
          bool evac;                                     
          if(fillstart>0)
               evac = julian<fillstart || julian>334;
          else
               evac = julian<182 || julian>334;
                
                
          bool spring = julian<182;
          double TDA_for;
          if(spring == true)
               TDA_for = accumulate(TDA_unreg.begin()+(year-2)*365+212,TDA_unreg.begin()+(year-2)*365+365,decltype(TDA_unreg)::value_type(0));
          else
               TDA_for = accumulate(TDA_unreg.begin()+(year-1)*365+212,TDA_unreg.begin()+(year-1)*365+365,decltype(TDA_unreg)::value_type(0));
          TDA_for = TDA_for*cfs_kaf_conversion_factor;
          
          double min_discharge=30000*cfs_kaf_conversion_factor;
          double GCL_fctarget=getGCL_fctarget(julian,TDA_for);
          double VPDR_GCL=min_discharge;
          GCL_sfore[t]=min((inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor+ storage[i][t-1] - VPDR_GCL,s_max[i]);
          GCL_VRC[julian-1]= GCL_sfore[t];
          
          
          double GCL_VRCLL=getGCL_VRCLL(julian);

          if(julian>=213 && julian <=365)
          {
                  GCL_ORC[julian-1]= max(GCL_CRC[julian-1][2]*RC_conversion_factor,GCL_ARC[julian-1][2]*RC_conversion_factor);
                  GCL_ORC[julian-1]= min(GCL_ORC[julian-1],GCL_fctarget);
          }
          else if(julian >=1 && julian <= 212)
          {
                  GCL_ORC[julian-1]= min(GCL_VRC[julian-1],max(GCL_CRC[julian-1][2]*RC_conversion_factor,GCL_ARC[julian-1][2]*RC_conversion_factor));
                  GCL_ORC[julian-1]= min(GCL_ORC[julian-1],GCL_fctarget);
                  GCL_ORC[julian-1]= max(GCL_ORC[julian-1],GCL_VRCLL);
          }
          else
                 GCL_ORC[julian-1]=GCL_fctarget;
          
          Storage_left[t]=GCL_fctarget-GCL_ORC[julian-1];
          double s_ti=storage[i][t-1];
          double d_ti;
          if(evac == true)
          {
              d_ti = max(min_discharge, s_ti + (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor -  GCL_ORC[julian-1]);
              d_ti=min(d_ti,max_discharge_gcl[months[t-1]-1]*cfs_kaf_conversion_factor); 
              double overflow = max(0.0,d_ti - power_max[9]);
              double stor_avail = GCL_fctarget - s_ti;
              if(overflow>0 && stor_avail<=0)
              {
                  
                 if(s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor - d_ti > s_max[9])
                        d_ti = s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor  - s_max[9];

                  powerflow[i][t] = power_max[9];
                  spill[i][t] = d_ti - power_max[9];
              }
              else if(overflow>0 && stor_avail>0)
              {
                  d_ti = max(d_ti - stor_avail,power_max[9]);
                  if(s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor - d_ti > s_max[9])
                        d_ti = s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor  - s_max[9];
                  powerflow[i][t] = power_max[9];
                  spill[i][t] = max(0.0,d_ti - power_max[9]);
              }
              else
                  powerflow[i][t] = d_ti;
              
              s_ti = s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor - d_ti;
          }
           else
          {
              d_ti = max(min_discharge, s_ti + (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor -  GCL_ORC[julian-1]);
              d_ti=min(d_ti,max_discharge_gcl[months[t-1]-1]*cfs_kaf_conversion_factor); 
              double overflow = max(0.0,d_ti - power_max[9]);
              double stor_avail = s_max[i] - s_ti;
              if(overflow>0 && stor_avail<=0)
              {
                  if(s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor - d_ti > s_max[9])
                        d_ti = s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor  - s_max[9];
                  powerflow[i][t] = power_max[9];
                  spill[i][t] = d_ti - power_max[9];
              }
              else if(overflow>0 && stor_avail>0)
              {
                    d_ti = max(d_ti - stor_avail,power_max[9]);
                    if(s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor - d_ti > s_max[9])
                        d_ti = s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor  - s_max[9];
                    powerflow[i][t] = power_max[9];
                    spill[i][t] = max(0.0,d_ti - power_max[9]);
              }
              else
                    if(s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor - d_ti > s_max[9])
                        d_ti = s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor  - s_max[9];
                    powerflow[i][t] = d_ti;

                  
                  s_ti = s_ti+ (inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1])*cfs_kaf_conversion_factor - d_ti;

          }   
             total_inflow[i][t]=(inflow[i][t] + spokane_flow[t-1] + discharge[1][t-1] + discharge[17][t-1] + discharge[22][t-1]);
             discharge[i][t]=max(min_discharge,d_ti)*(1/cfs_kaf_conversion_factor);
             storage[i][t]=s_ti;
             powerflow[i][t]=powerflow[i][t]*(1/cfs_kaf_conversion_factor);
             spill[i][t]=spill[i][t]*(1/cfs_kaf_conversion_factor); 
             calcGeneration(9,t);

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
                        //discharge[res_no][t] = min( max_discharge , max( min_discharge , uu[opt_index[res_no]]*cfs_kaf_conversion_factor) );
                        discharge[res_no][t]=modelled_discharge[res_no][t]*cfs_kaf_conversion_factor;
                        powerflow[res_no][t] = min(discharge[res_no][t],power_max[res_no]);
                        spill[res_no][t]    = max(0.0,discharge[res_no][t]- powerflow[res_no][t]);
                        storage[res_no][t]   = max(0.0,(s_ti-discharge[res_no][t]));
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
                        discharge[res_no][t]=modelled_discharge[res_no][t];
                        storage[res_no][t]=max(0.0,(s_ti-d_ti));
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
                       vector<double> uu;
                       vector<double> input;
                       vector<double> GCL_sfore(sim_days+1,0);
                       vector<double> Storage_left(sim_days+1,0);
                       vector<double> GCL_VRC(365,0);
                       vector<double> GCL_ORC(365,0);
                       for(int t=1;t<=sim_days;t++)
                       { 
                            
                            int julian = (int)julians[t-1];

                            for (auto res : reservoirs_network)
                            {
                            if(opt_reservoirs.find(res)!=opt_reservoirs.end())
                                model_storage_res(res,t,uu);
                            else
                                {
                                    if(res==7)
                                        modelALB(t);
                                    else if(res==9)
                                        modelGCL(res,t,GCL_sfore,Storage_left,GCL_VRC,GCL_ORC,spokane_flow);
                                    else
                                        model_run_of_river(res,t);
                                }
                            }

                         }


}



void initialize()
{
    fillSpokaneFlow();
    total_inflow[9][0]=modelled_discharge[1][0]+modelled_discharge[17][0]+modelled_discharge[22][0]+spokane_flow[0]+inflow[9][1];
    total_inflow[35][0]=modelled_discharge[12][0]+modelled_discharge[43][0]+inflow[35][1]
                        + inflow[49][1]+ inflow[50][1]+ inflow[51][1]+ inflow[52][1]+ inflow[53][1];
}

void initializeGCL(unsigned int rank_no,unsigned int master_no)
{

        int folder_no=master_no;
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
        
        vector<vector<double> > ICF_file=utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/icf.txt",sim_days/365,1,1,0);
        utils::transpose(ICF_file);
        ICFs=ICF_file[0];

}

void preprocess(unsigned int rank_no,unsigned int master_no)
{

        
        utils::transpose(modelled_discharge);
        //utils::transpose(modelled_powerflow);
        utils::transpose(modelled_storage);
        utils::transpose(modelled_generation);
        utils::transpose(inflow);
        utils::transpose(synth_weather);
        utils::transpose(maxhydro_matrix);
        utils::transpose(maxstorage_matrix);
        utils::transpose(power_constants);
        utils::transpose(discharge);
        utils::transpose(storage);
        utils::transpose(modelled_demands);
        utils::transpose(extra_generation);
        utils::transpose(ca_surrogate_df);
        utils::transpose(pnw_surrogate_df);
        utils::transpose(ca_surrogate);
        utils::transpose(pnw_surrogate);
        power_max=maxhydro_matrix[1];
        s_max=maxstorage_matrix[1];
        hydraulic_head=power_constants[1];
        powerhouse_efficiency=power_constants[2];
        initializeGCL(rank_no, master_no);
        
}

void writeSimulationOutput(unsigned int rank_no,unsigned int master_no)
{

        //cout<<"inside writing simulate"<<endl;
        utils::writeMatrix(generation,"../output_opt_historical/policy"+to_string(rank_no)+"/generation_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(storage,"../output_opt_historical/policy"+to_string(rank_no)+"/storage_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(discharge,"../output_opt_historical/policy"+to_string(rank_no)+"/discharge_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(powerflow,"../output_opt_historical/policy"+to_string(rank_no)+"/powerflow_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        //utils::writeMatrix(modelled_powerflow,"../output_opt_historical/policy"+to_string(rank_no)+"/modelled_powerflow"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(spill,"../output_opt_historical/policy"+to_string(rank_no)+"/spill_hyssr"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(heads,"../output_opt_historical/policy"+to_string(rank_no)+"/head_hyssr"+to_string(master_no)+".txt", sim_days+1, M_H);
        utils::writeMatrix(stage_hyssr,"../output_opt_historical/policy"+to_string(rank_no)+"/stage_hyssr"+to_string(master_no)+".txt", sim_days+1, M_H);
        utils::writeMatrix(tailwater_hyssr,"../output_opt_historical/policy"+to_string(rank_no)+"/tailwater_hyssr"+to_string(master_no)+".txt", sim_days+1, M_H);
        utils::writeMatrix(modelled_discharge,"../output_opt_historical/policy"+to_string(rank_no)+"/modelled_discharge"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(modelled_storage,"../output_opt_historical/policy"+to_string(rank_no)+"/modelled_storage"+to_string(master_no)+".txt", sim_days+1, M);
        utils::writeMatrix(modelled_generation,"../output_opt_historical/policy"+to_string(rank_no)+"/modelled_generation"+to_string(master_no)+".txt", sim_days+1, M+1);
        utils::writeMatrix(inflow,"../output_opt_historical/policy"+to_string(rank_no)+"/inflow"+to_string(master_no)+".txt", sim_days+1, M_comp);
        utils::writeMatrix(total_inflow,"../output_opt_historical/policy"+to_string(rank_no)+"/total_inflow"+to_string(master_no)+".txt", sim_days+1, M);
        //cout<<"write simulatio output "<<to_string(rank_no)+" "<<to_string(master_no)<<endl;
        utils::transpose(generation_fcrps);
        utils::writeMatrix(generation_fcrps,"../output_opt_historical/policy"+to_string(rank_no)+"/generation_fcrps"+to_string(master_no)+".txt", sim_days, M_HR);
        utils::transpose(generation_fcrps);
        utils::transpose(generation_shoal);
        utils::writeMatrix(generation_shoal,"../output_opt_historical/policy"+to_string(rank_no)+"/generation_shoal"+to_string(master_no)+".txt", sim_days, M_H);
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
          utils::writeMatrix(Yt_matrix,"../output_opt_historical/policy"+to_string(rank_no)+"/flood_level_sim"+to_string(master_no)+".txt", sim_days, 1);
          }
          //std::transform(flood_level.begin(), flood_level.end(), flood_level.begin(),[=](double i) { return i/sim_years; });
          flood_level.erase (flood_level.begin());
          return flood_level;

}

/*
vector<double> getCaisoPrices(unsigned int rank_no,unsigned int master_no)
{

          int j;
          
          vector<double> caiso_prices_coeff=utils::loadVectorFromDataFile("../data/cons_opt/models_surrogate_wload/ca_linear_coef_list.txt");
          vector<double> arima_caiso_coeff=utils::loadVectorFromDataFile("../data/cons_opt/models_surrogate_wload/ca_arima_coef_list.txt");
          
          vector<vector<double> > residuals_matrix=utils::readMatrixFromDataFile("../data/cons_opt/models_surrogate_wload/ca_cons_residuals_cpp.txt");
          vector<double> quantiles_list=utils::loadVectorFromDataFile("../data/cons_opt/models_surrogate_wload/ca_quantiles_list.txt");
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
                   + caiso_prices_coeff[2]*modelled_demands[3][t]
                   + caiso_prices_coeff[3]*modelled_demands[4][t]
                   + caiso_prices_coeff[4]*modelled_demands[5][t]
                   + caiso_prices_coeff[5]*modelled_solar[0][t]
                   + caiso_prices_coeff[6]*modelled_wind[0][t]
                   + caiso_prices_coeff[7]*modelled_cahydro[0][t]
                   + caiso_prices_coeff[8]*(utils::columnSumRow(generation_fcrps,t)+extra_generation[0][t])
                   + caiso_prices_coeff[9]*sin( 2*pi*julians[t]/7)
                   + caiso_prices_coeff[10]*cos( 2*pi*julians[t]/7);

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
            srand(t);
            w_t[t+2]=residuals_matrix[j][rand() % residuals_matrix[j].size()];
            e_t[t+2]=e_t[t+2]+w_t[t+2];      
            Y_t[t]=Y_t[t]+e_t[t+2];
           
            Y_t[t]=max(10.0,Y_t[t]);
        }
        if(write_output)
        {
            vector<vector<double> >Yt_matrix(1,vector<double> (sim_days,0));
            Yt_matrix[0]=Y_t;
            //utils::writeMatrix(Yt_matrix,"../output/Caiso_prices0.txt", sim_days, 1);
            utils::writeMatrix(Yt_matrix,"../output_opt_historical/policy"+to_string(rank_no)+"/Caiso_prices"+to_string(master_no)+".txt", sim_days, 1);
        }
        return Y_t;


}*/

vector<vector<double >> getMidCaisoPrices(unsigned int rank_no,unsigned int master_no)
{

        cout<<"inside getMidCaisoPrices"<<endl;
        cout<<pnw_surrogate.size()<<endl;
        cout<<pnw_surrogate[0].size()<<endl;

          int j;
          double min_fossils_constraint=5460;
          vector<double> prices_coeff=utils::loadVectorFromDataFile("../data/surrogates_with_curtailment/MidC/midc_linear_coef_cpp.txt");
          vector<double> arima_coeff=utils::loadVectorFromDataFile("../data/surrogates_with_curtailment/MidC/midc_arima_coef_cpp.txt");
          vector<vector<double> > residuals_matrix=utils::readMatrixFromDataFile("../data/surrogates_with_curtailment/MidC/midc_residuals_loadsim_cpp.txt");
          vector<double> quantiles_list=utils::loadVectorFromDataFile("../data/surrogates_with_curtailment/MidC/midc_quantiles_loadsim_cpp.txt");
          vector<double> Y_t(sim_days,0);
          vector<double> e_t(sim_days+2,0);
          vector<double> w_t(sim_days+2,0);
          vector<int> bucket_index(sim_days,-1);
          w_t[0]=0;//taken the last two values from ARIMA(2,0,2) model
          w_t[1]=0;
          e_t[0]=0;//taken the last value from the predictor's model
          e_t[1]=0;
          double fitted_value;
          

          vector<double> prices_coeff_ca=utils::loadVectorFromDataFile("../data/surrogates_with_curtailment/CAISO/caiso_linear_coef_cpp.txt");
          vector<double> arima_coeff_ca=utils::loadVectorFromDataFile("../data/surrogates_with_curtailment/CAISO/caiso_arima_coef_cpp.txt");
          vector<vector<double> > residuals_matrix_ca=utils::readMatrixFromDataFile("../data/surrogates_with_curtailment/CAISO/caiso_residuals_loadsim_cpp.txt");
          vector<double> quantiles_list_ca=utils::loadVectorFromDataFile("../data/surrogates_with_curtailment/CAISO/caiso_quantiles_loadsim_cpp.txt");
          vector<double> Y_t_ca(sim_days,0);
          vector<double> e_t_ca(sim_days+2,0);
          vector<double> w_t_ca(sim_days+2,0);
          vector<int> bucket_index_ca(sim_days,-1);
          w_t_ca[0]=0;//taken the last two values from ARIMA(1,0,2) model
          w_t_ca[1]=0;
          e_t_ca[0]=0;//taken the last value from the predictor's model
          e_t_ca[1]=0;
          double fitted_value_ca;



        for( int t=0;t<sim_days;t++)
         {

           //w_t[t+2]=distribution(generator);
           double changed_hydro=(utils::columnSumRow(generation_fcrps,t)+extra_generation[0][t]);
           double difference;
           double pnw_hydro=pnw_surrogate[5][t];
           double pnw_fossils=pnw_surrogate[7][t];
           double pnw_exports=pnw_surrogate[1][t];
           double ca_imports=ca_surrogate[6][t];
           double ca_fossils=ca_surrogate[7][t];

           if(t<=5)
           {

            cout<<"t "<<endl;
            cout<<"changed_hydro"<<changed_hydro<<endl;
            cout<<"pnw_hydro "<<pnw_hydro<<endl;
            cout<<"pnw_fossils "<<pnw_fossils<<endl;
            cout<<"pnw_exports "<<pnw_exports<<endl;
            cout<<"ca_imports "<<ca_imports<<endl;
            cout<<"ca_fossils "<<ca_fossils<<endl;
            cout<<"pnw_surrogate[0][t] "<<pnw_surrogate[0][t]<<endl;
           }






           if(changed_hydro<=pnw_hydro)
           {
                    difference=pnw_hydro-changed_hydro;
                    pnw_fossils=pnw_fossils+difference;
                    if(t<=5)
                       {
                        cout<<"inside first if"<<endl;
                        cout<<"t "<<endl;
                        cout<<"difference "<<difference<<endl;
                        cout<<"pnw_fossils "<<pnw_fossils<<endl;

                       }
           }
           else
           {
                    difference=changed_hydro-pnw_hydro;
                    pnw_fossils=pnw_fossils-difference;
                    if(t<=5)
                       {
                        cout<<"inside second if"<<endl;
                        cout<<"t "<<endl;
                        cout<<"difference "<<difference<<endl;
                        cout<<"pnw_fossils "<<pnw_fossils<<endl;

                       }
                    if(pnw_fossils<=min_fossils_constraint)
                    {

                        
                        double ca_exports_diff=min_fossils_constraint-pnw_fossils;
                        pnw_fossils=min_fossils_constraint;
                        ca_exports_diff=max(ca_exports_diff,100800.0);
                        pnw_exports=pnw_exports+ca_exports_diff;
                        ca_imports=ca_imports+ca_exports_diff;
                        ca_fossils=ca_fossils-ca_exports_diff;
                        if(t<=5)
                       {
                        cout<<"inside third if"<<endl;
                        cout<<"t "<<endl;
                        cout<<"difference "<<difference<<endl;
                        cout<<"pnw_fossils "<<pnw_fossils<<endl;
                        cout<<"ca_exports_diff "<<ca_exports_diff<<endl;
                        cout<<"pnw_exports "<<pnw_exports<<endl;
                        cout<<"ca_imports "<<ca_imports<<endl;

                       }
                        if(ca_fossils<=0)
                        {
                        double hydro_curt=-1*(ca_fossils);
                        changed_hydro=changed_hydro-hydro_curt;
                        pnw_exports=pnw_exports-hydro_curt;
                        ca_imports=ca_imports-hydro_curt;
                        ca_fossils=0;
                        }
                    }


            }

           e_t[t+2]=arima_coeff[0]*e_t[(t+2)-1]+arima_coeff[1]*e_t[(t+2)-2]+arima_coeff[2]*w_t[(t+2)-1]+arima_coeff[3]*w_t[(t+2)-2];
           Y_t[t] =  prices_coeff[0] 
                   + prices_coeff[1]*pnw_surrogate[0][t] 
                   + prices_coeff[2]*pnw_exports
                   + prices_coeff[3]*pnw_surrogate[4][t] 
                   + prices_coeff[4]*(changed_hydro+pnw_surrogate[9][t])
                   + prices_coeff[5]*pnw_surrogate[6][t] 
                   + prices_coeff[6]*sin( 2*pi*julians[t]/7)
                   + prices_coeff[7]*cos( 2*pi*julians[t]/7);
    

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
            srand(t);       
            w_t[t+2]=residuals_matrix[j][rand() % residuals_matrix[j].size()];
            e_t[t+2]=e_t[t+2]+w_t[t+2];      
            Y_t[t]=Y_t[t]+e_t[t+2];
            Y_t[t]=max(10.0,Y_t[t]);

            if(t<=5)
           {

            cout<<"t "<<endl;
            cout<<"changed_hydro"<<changed_hydro<<endl;
            cout<<"e_t[t+2] "<<e_t[t+2]<<endl;
            cout<<"w_t[t+2] "<<w_t[t+2]<<endl;
            cout<<"pnw_exports "<<pnw_exports<<endl;
            cout<<"Y_t[t] "<<Y_t[t]<<endl;

           }



            e_t_ca[t+2]=arima_coeff_ca[0]*e_t[(t+2)-1]+arima_coeff_ca[1]*e_t[(t+2)-2]+arima_coeff_ca[2]*w_t[(t+2)-1]+arima_coeff_ca[3]*w_t[(t+2)-2];
            Y_t_ca[t] =  prices_coeff_ca[0] 
                   + prices_coeff_ca[1]*ca_surrogate[0][t] 
                   + prices_coeff_ca[2]*ca_surrogate[1][t]
                   + prices_coeff_ca[3]*ca_surrogate[3][t] 
                   + prices_coeff_ca[4]*ca_surrogate[4][t]
                   + prices_coeff_ca[5]*(ca_surrogate[5][t]+ca_surrogate[9][t])
                   + prices_coeff_ca[6]*ca_imports 
                   + prices_coeff_ca[7]*sin( 2*pi*julians[t]/7)
                   + prices_coeff_ca[8]*cos( 2*pi*julians[t]/7);
    

            fitted_value_ca= Y_t_ca[t]+e_t_ca[t+2];
                    for (j=0;j<quantiles_list_ca.size();j++)
                      {
                        if(fitted_value_ca<quantiles_list_ca[j]){
                          break;
                        }
                        
                      }
            if(j>=20)
                j=19;
            bucket_index_ca[t]=j; 
            srand(t);       
            w_t_ca[t+2]=residuals_matrix_ca[j][rand() % residuals_matrix_ca[j].size()];
            e_t_ca[t+2]=e_t_ca[t+2]+w_t_ca[t+2];      
            Y_t_ca[t]=Y_t_ca[t]+e_t_ca[t+2];
            Y_t_ca[t]=max(10.0,Y_t_ca[t]);
            if(t<=5)
           {

            cout<<"t "<<endl;
            cout<<"changed_hydro"<<changed_hydro<<endl;
            cout<<"e_t_ca[t+2] "<<e_t_ca[t+2]<<endl;
            cout<<"w_t_ca[t+2] "<<w_t_ca[t+2]<<endl;
            cout<<"ca_imports "<<ca_imports<<endl;
            cout<<"Y_t_ca[t] "<<Y_t_ca[t]<<endl;

           }


            pnw_surrogate[1][t]=pnw_exports;
            pnw_surrogate[5][t]=changed_hydro;
            pnw_surrogate[7][t]=pnw_fossils;
            pnw_surrogate[8][t]=Y_t[t];

            ca_surrogate[6][t]=ca_imports;
            ca_surrogate[7][t]=ca_fossils;
            ca_surrogate[8][t]=Y_t_ca[t];     

        }

        vector<vector<double> >Yt_matrix(2,vector<double> (sim_days,0));
        Yt_matrix[0]=Y_t;
        Yt_matrix[1]=Y_t_ca;

        if(write_output)
        {

            utils::writeMatrix(Yt_matrix,"../output_opt_historical/policy"+to_string(rank_no)+"/prices"+to_string(master_no)+".txt", sim_days, 2);
            utils::writeMatrix(pnw_surrogate,"../output_opt_historical/policy"+to_string(rank_no)+"/pnw_surrogate"+to_string(master_no)+".txt", sim_days, 10);
            utils::writeMatrix(ca_surrogate,"../output_opt_historical/policy"+to_string(rank_no)+"/ca_surrogate"+to_string(master_no)+".txt", sim_days, 10);
        }
        
        return Yt_matrix;
}


vector<double> getRevenueBPA(unsigned int rank_no,unsigned int master_no)
{
   
   int folder_no=master_no;
   vector<vector<double> > MidC(1,vector<double> (sim_days,0));
   vector<vector<double> > CAISO(1,vector<double> (sim_days,0));
   //vector<double>MidC=utils::loadVectorFromDataFile("midc_prices_53.txt");
   vector<vector<double> > prices=getMidCaisoPrices(rank_no,master_no);
   MidC[0]=prices[0];
   CAISO[0]=prices[1];
   vector<double> midc_prices = {MidC[0].begin() + 122, MidC[0].end() - 243}; 
   //vector<double> midc_prices = {MidC.begin() + 122, MidC.end() - 243}; 
   vector<double> caiso_prices = {CAISO[0].begin() + 122, CAISO[0].end() - 243};
   //cout<<"meanVector(midc_prices) "<<meanVector(midc_prices)<<"rank_no "<<rank_no<<endl;
   //cout<<"meanVector(caiso_prices]) "<<meanVector(caiso_prices)<<"rank_no "<<rank_no<<endl;

   
   //vector<double> load_vector=utils::loadVectorFromDataFile("../data/Revenue Model/load_2018.txt");
   vector<double> load_vector=utils::loadVectorFromDataFile("../data/Revenue Model/load_avg.txt");
   double custom_redux=0;
   double PF_load_y=load_vector[13]-custom_redux*load_vector[13];
   double IP_load_y=load_vector[3]-custom_redux*load_vector[3];
   double ET_load_y=load_vector[14];
   
   
   
   vector<vector<double> > extra_bpa_dams=utils::loadMatrix("../data/optimize_four_1000/"+to_string(folder_no)+"/extra_bpa.txt",bpa_rev_sim_days,1,365+extra_start,start_col);
   utils::transpose(extra_bpa_dams);
   //cout<<"extra_bpa_dams[0].size() "<<extra_bpa_dams[0].size()<<"rank_no "<<rank_no<<endl;
   //cout<<"meanVector(extra_bpa_dams[0]) "<<meanVector(extra_bpa_dams[0])<<"rank_no "<<rank_no<<endl;
   //vector<double> extra_gen=utils::columnSumMatrix(extra_bpa_dams);
   vector<double> bpa_gen=utils::columnSumMatrix(generation_shoal);
   vector<double> bpa_gen_sub = {bpa_gen.begin() + 122, bpa_gen.end() - 243};
   //cout<<"bpa_gen_sub.size() "<<bpa_gen_sub.size()<<"rank_no "<<rank_no<<endl;
   //cout<<"meanVector(bpa_gen_sub) "<<meanVector(bpa_gen_sub)<<"rank_no "<<rank_no<<endl;
   vector<double> BPA_hydro(bpa_gen_sub.size(),0);
   std::transform(bpa_gen_sub.begin(),bpa_gen_sub.end(), extra_bpa_dams[0].begin(),BPA_hydro.begin(),std::plus<double>());
   std::transform(BPA_hydro.begin(), BPA_hydro.end(), BPA_hydro.begin(),[=](double i) { return i/24; });
   std::transform(BPA_hydro.begin(), BPA_hydro.end(), BPA_hydro.begin(),[](double value) { return std::min(value, 45000.0); });
   double maxhydro = *max_element(BPA_hydro.begin(), BPA_hydro.end());
   cout<<"maxhydro "<<maxhydro<<endl;
   //cout<<"BPA_hydro.size() "<<BPA_hydro.size()<<"rank_no "<<rank_no<<endl;
   //cout<<"meanVector(BPA_hydro) "<<meanVector(BPA_hydro)<<"rank_no "<<rank_no<<endl;
   
   //cout<<" in revenue model 1"<<endl;
   //vector<double> net_resources=utils::loadVectorFromDataFile("../data/Revenue Model/BPA_net_resources_2018.txt");
   vector<double> net_resources=utils::loadVectorFromDataFile("../data/Revenue Model/BPA_net_resources_avg.txt");
   double Nuc_y=net_resources[7];
   double Wind_y=net_resources[8];
   double Purch_y=net_resources[10];
   
   //vector<double> PF_rates=utils::loadVectorFromDataFile("../data/Revenue Model/PF_rates_2018.txt");
   //vector<double> IP_rates=utils::loadVectorFromDataFile("../data/Revenue Model/IP_rates_2018.txt");
   vector<double> PF_rates=utils::loadVectorFromDataFile("../data/Revenue Model/PF_rates_avg.txt");
   vector<double> IP_rates=utils::loadVectorFromDataFile("../data/Revenue Model/IP_rates_avg.txt");
   /*
   vector<vector<double> > BPAT_load=utils::loadMatrix("../data/cons_opt/Datafiles/daily_load_bpa.txt",sim_days,1,start_row,start_col);
   BPAT_load= utils::matrixReshape( BPAT_load,sim_days,1);
   utils::transpose(BPAT_load);
   */
   vector<double> load_bpa=modelled_demands[0];
   vector<double> BPAT_load={load_bpa.begin() + 365, load_bpa.begin() + 365+ bpa_rev_sim_days};
   //cout<<"BPAT_load.size() "<<BPAT_load.size()<<"rank_no "<<rank_no<<endl;
   //cout<<"meanVector(BPAT_load) "<<meanVector(BPAT_load)<<"rank_no "<<rank_no<<endl;
   //utils::transpose(modelled_wind);
   vector<double> wind_bpa=utils::loadVectorFromDataFile("../data/optimize_four_1000/"+to_string(folder_no)+"/wind_surrogate_bpa.txt");
   vector<double> BPAT_wind={wind_bpa.begin() + 365+extra_start, wind_bpa.begin() +365 +bpa_rev_sim_days+extra_start}; 
   //cout<<"BPAT_wind.size() "<<BPAT_wind.size()<<"rank_no "<<rank_no<<endl;
   //cout<<"meanVector(BPAT_wind) "<<meanVector(BPAT_wind)<<"rank_no "<<rank_no<<endl;
   //std::transform(modelled_wind[0].begin(), modelled_wind[0].end(), BPAT_wind.begin(),[=](double i) { return i *(0.766/1.766); });
   

   //double costs_y=utils::loadVectorFromDataFile("../data/Revenue Model/BPA_yearly_costs.txt")[8];
   double costs_y=2229980000.0;
   
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
                                   cout<<"Net_res1 is negative"<<endl;
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
      cout<<"inside BPA rev"<<endl;
        if(write_output)
        {
            vector<vector<double> >Yt_matrix(3,vector<double> (num_years,0));
            Yt_matrix[0]=BPA_Net_rev_y;
            Yt_matrix[1]=CRAC_rev;
            Yt_matrix[2]=TTP;
            utils::writeMatrix(Yt_matrix,"../output_opt_historical/policy"+to_string(rank_no)+"/BPA_Net_Revenue"+to_string(master_no)+".txt", num_years, 3);
            
            vector<vector<double> >Ytp1_matrix(5,vector<double> (num_years+1,0));
            Ytp1_matrix[0]=Reserves;
            Ytp1_matrix[1]=Remaining_BA;
            Ytp1_matrix[2]=TF1;
            Ytp1_matrix[3]=TF2;
            Ytp1_matrix[4]=CRAC_y;
            utils::writeMatrix(Ytp1_matrix,"../output_opt_historical/policy"+to_string(rank_no)+"/BPA_Net_Revenue_plus_one"+to_string(master_no)+".txt", num_years+1, 5);

           
            vector<vector<double> >Sd_matrix(8,vector<double> (load_ratio.size(),0));
            Sd_matrix[0]=SD;
            Sd_matrix[1]=PF_rev;
            Sd_matrix[2]=IP_rev;
            Sd_matrix[3]=P;
            Sd_matrix[4]=SS;
            Sd_matrix[5]=BPA_rev_d;
            Sd_matrix[6]=crac;
            Sd_matrix[7]=load_ratio;
            utils::writeMatrix(Sd_matrix,"../output_opt_historical/policy"+to_string(rank_no)+"/SD"+to_string(master_no)+".txt", load_ratio.size(), 8);
        }     
        
        BPA_Net_rev_y.erase (BPA_Net_rev_y.begin());
        return BPA_Net_rev_y;
}


double flood_level_frequency(vector<double> floodlevel)
{
        int count_flood=0;
        for(auto & flood : floodlevel)
        {
            if (flood>=17)
                count_flood+=1;
        }

    return (count_flood*1.0)/sim_years;
}

double minimize_fossils()
{
        vector<double> fossils(sim_days,0);
        //trans_losses=3*(Wind + BPA_hydro + Nuc)/100;
        std::transform(pnw_surrogate[7].begin(), pnw_surrogate[7].end(), ca_surrogate[7].begin(),fossils.begin(),std::plus<double>());

    return meanVector(fossils);
}

double waterTempConstraints(unsigned int rank_no,unsigned int master_no)
{

        vector<double> Y_t(sim_days);
        vector<double> watertmp_coeff=utils::loadVectorFromDataFile("../data/water_temp_model/water_temp_coef_list.txt");
        double waterTempViolations=0;
        for(int t=0;t<sim_days;t++)
        {
            Y_t[t] =  watertmp_coeff[0] 
                   + watertmp_coeff[1]*sin( 2*pi*julians[t]/365)
                   + watertmp_coeff[2]*cos( 2*pi*julians[t]/365)
                   + watertmp_coeff[3]*total_inflow[35][t]/1000
                   + watertmp_coeff[4]*((synth_weather[10][t]*9)/5+32)
                   + watertmp_coeff[5]*synth_weather[11][t];
            int julian=(int)julians[t];
            if(60<=julian<=273)
                {
                    if(Y_t[t]>=71.06)
                        waterTempViolations+=((Y_t[t]-71.06)*(Y_t[t]-71.06));
                }
            
        }
        if(write_output)
        {
            vector<vector<double> >Yt_matrix(1,vector<double> (sim_days,0));
            Yt_matrix[0]=Y_t;
            //utils::writeMatrix(Yt_matrix,"../output/MidC_prices0.txt", sim_days, 1);
            utils::writeMatrix(Yt_matrix,"../output_opt_historical/policy"+to_string(rank_no)+"/WaterTemp"+to_string(master_no)+".txt", sim_days, 1);
        }

        return waterTempViolations;

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
        //J.push_back(-1*meanVector(maxRenewablesGeneration()));
        vector<double> floodlevel=floodLevel(rank_no,master_no,discharge[42]);
        //vector<double> floodlevel_historical=floodLevel(rank_no,master_no,modelled_discharge[42]);
        //cout<<"historical flood "<<*max_element(floodlevel_historical.begin(), floodlevel_historical.end())<<endl;
        J.push_back(*max_element(floodlevel.begin(), floodlevel.end()));
        J.push_back((meanVector(getRevenueBPA(rank_no,master_no))/1000000)*-1);
        cout<<"sdNegativeSumOpt "<<sdNegativeSumOpt<<endl;
        J.push_back(minimize_fossils());
        //J.push_back(flood_level_frequency(floodlevel));
        //cout<<"flood frequency"<<flood_level_frequency(floodlevel)<<endl;
        //J.push_back(waterTempConstraints(rank_no,master_no));
        //cout<<"waterTempConstraints(rank_no,master_no)"<<waterTempConstraints(rank_no,master_no)<<endl;


        return J;


}



void evaluate( double* obj,unsigned int inp_rank_no,unsigned int inp_master_no){

                unsigned int rank_no= inp_rank_no,master_no=inp_master_no;
                
                readDataFiles(rank_no, master_no);
                preprocess(rank_no, master_no);
                makeSimIndexing();
                makeinverseSimIndexing();
                makeOPTIndexing();
                makeinverseOPTIndexing();
                makeHYSSRIndexing();
                makeinverseHYSSRIndexing();
                initialize();
                makeFlowRoute();
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
                    utils::writeMatrix(obj_matrix,"../output_opt_historical/policy"+to_string(rank_no)+"/objectives"+to_string(master_no)+".txt", Nobj, 1);

                    
                }



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

int main(int argc, char* argv[])
{
         
        MyClass obj1=MyClass();
        vector<double> objs(4,0);
        double* obj_array = &objs[0];
        //int master_rank=atoi(argv[1]);
        int master_rank=2;
        unsigned int master_no=master_rank;
        obj1.evaluate(obj_array,52,master_no);
         
         return 0;
}
