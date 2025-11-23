#include<iostream>
#include<cmath>
#include<cstdlib>
#include<vector>
#include<fstream>
#include<ctime>

using namespace std;
#include "unfragmented_heap.h"
#include "mtrand.h"

void read_params(ifstream& in, vector<double>& value, int& no)
{	string name;
	double d;
	int total=0;
	for (int i=0;i<20;++i)value.push_back(0);
	while (!in.eof())
	{	in>>name;
		if(name=="nS0")
		{	in>>d;
			value[0]=d;
			total++;
		}
		
		else if(name=="nM")
		{	in>>d;
			value[1]=d;
			total++;
		}

		else if(name=="nD")
		{	in>>d;
			value[2]=d;
			total++;
		}
		
		else if(name=="kpM")
		{	in>>d;
			value[3]=d;
			total++;
		}

		else if(name=="kmM")
		{	in>>d;
			value[4]=d;
			total++;
		}
		
		else if(name=="kcatM1")
		{	in>>d;
			value[5]=d;
			total++;
		}
	
		else if(name=="kcatM2")
		{	in>>d;
			value[6]=d;
			total++;
		}

		else if(name=="kpD")
		{	in>>d;
			value[7]=d;
			total++;
		}

		else if(name=="kmD")
		{	in>>d;
			value[8]=d;
			total++;
		}
	
		else if(name=="kcatD")	
		{	in>>d;
			value[9]=d;
			total++;
		}
		
		else if(name=="Q")
		{	in>>d;
			value[10]=d;
			total++;
		}
		
		else if(name=="d1")
		{	in>>d;
			value[11]=d;
			total++;
		}

		else if(name=="d2")
		{	in>>d;
			value[12]=d;
			total++;
		}
	}
		no=total-1;
}

int main (int argc, char ** argv)
{	
	if (argc != 5)
	{
		exit (1) ;
    }
	
	ifstream inFile;
	ofstream outFile, outFile1;
	double tmax = atof(argv[1]);
	int Nsim = atoi(argv[2]);
	inFile.open(argv[3]);
	outFile.open(argv[4]);
	int sd = (int)clock();
	MTRand mtrand;
	mtrand.seed(sd);
	
//Reading the parameters from variable file "var.txt"
	vector<double>parameter;
	int nparam;
	read_params(inFile,parameter,nparam);
	//cout<<"Number of parameters = "<<nparam<<endl;


	
//Rate constants
	double kpM	= (double)parameter[3];
	double kmM	= (double)parameter[4];
	double kcatM1	= (double)parameter[5];
	double kcatM2	= (double)parameter[6];	
	double kpD		= (double)parameter[7];
	double kmD		= (double)parameter[8];
	double kcatD	= (double)parameter[9];
	double Q	= (double)parameter[10];
	double d1	= (double)parameter[11];
	double d2	= (double)parameter[12];

 	int x = 0;
	
	double time, delta, r1, r2 ;
	
   	double S0_syn, S0_deg, Slt4_deg, Sgeq4_deg, MSeq0_deg, MSlt4_deg, MSgeq4_deg, DSlt4_deg, DSgeq4_deg, S0M_bind, MS0_diss, Slt4M_bind, MSlt4_diss, Sgeq4M_bind, MSgeq4_diss, Slt4D_bind, DSlt4_diss, Sgeq4D_bind, DSgeq4_diss, MSeq0cat1, MSlt4cat2, MSgeq4cat2, DSlt4cat, DSgeq4cat, prop_tot;

	double prop_this;
	double epsilon;
    
   	int nthS;
	bool reacted;
	
	double avg_S0, avg_Sprime, avg_Stot;
	double avg_M, avg_D, tot;
	avg_S0 = 0;
	avg_Sprime = 0;
	avg_Stot = 0;
	avg_M = 0;
	avg_D = 0;
	int Sinit, Dinit;
	float Minit;
	Sinit 	= (int)parameter[0]; 		
	Minit 	= (float)parameter[1];
	Dinit	= (int)parameter[2];
	for (int p=0;p<4;++p)				//loop over initial concentration of M
	{	
		vector<int>all_lengths;
		for (int k=0;k<Nsim;++k)		//loop over Nsim
		{	//Initializing variables
			time=0.0;
			//Initial values of the S0, M, D
						
			int nS0	= Sinit;		
			int nD	= Dinit;
			int nM = round(Minit);
			
			unfragmented_heap<int>Seq0;		
			unfragmented_heap<int>MSeq0;		
			unfragmented_heap<int>Slt4;
			unfragmented_heap<int>Sgeq4;
			unfragmented_heap<int>MSlt4;
			unfragmented_heap<int>MSgeq4;
			unfragmented_heap<int>DSlt4;
			unfragmented_heap<int>DSgeq4;
	
			for(int i=0;i<nS0;++i)Seq0.add(x);
			
			while (time<=(tmax+0.1*tmax))								// Run it for a time > tmax so that we can store the values at tmax
			{	S0_syn		= Q;                        				//1. Synthesis rate of S0
				S0_deg		= d1*(double)Seq0.size();      	   			//2. Degradation of S0
				Slt4_deg	= d1*(double)Slt4.size();           		//3. Degradation of S1, S2, S3
				Sgeq4_deg	= d2*(double)Sgeq4.size();          		//4. Degradation of S4, S5, ...
				MSeq0_deg	= d1*(double)MSeq0.size();          		//5. Degradation of MS0
				MSlt4_deg	= d1*(double)MSlt4.size();         			//6. Degradation of MS1, MS2, MS3
				MSgeq4_deg	= d2*(double)MSgeq4.size();         		//7. Degradation of MS4, MS5, ...
				DSlt4_deg	= d1*(double)DSlt4.size();          		//8. Degradation of DS1, DS2, DS3
				DSgeq4_deg	= d2*(double)DSgeq4.size();         		//9. Degradation of DS4, DS5, ...
				S0M_bind	= kpM*(double)nM*(double)Seq0.size();		//10. M binds to S0
				MS0_diss	= kmM*(double)MSeq0.size();         		//11. M dissociates from S0
				Slt4M_bind	= kpM*(double)nM*(double)Slt4.size();       //12. M binds to S1, S2, S3
				MSlt4_diss	= kmM*(double)MSlt4.size();        	 		//13. M dissociates from S1, S2, S3
				Sgeq4M_bind	= kpM*(double)nM*(double)Sgeq4.size();		//14. M binds to S4, S5, ...
				MSgeq4_diss	= kmM*(double)MSgeq4.size();				//15. M dissociates from S4, S5, ...
				Slt4D_bind	= kpD*(double)nD*(double)Slt4.size();		//16. D binds to S1, S2, S3
				DSlt4_diss	= kmD*(double)DSlt4.size();         		//17. D dissociates from S1, S2, S3
				Sgeq4D_bind	= kpD*(double)nD*(double)Sgeq4.size();		//18. D binds to S4, S5...
				DSgeq4_diss	= kmD*(double)DSgeq4.size();				//19. D dissociates from 
				MSeq0cat1	= kcatM2*(double)MSeq0.size();				//20. M catalyzes S0
				MSlt4cat2	= kcatM2*(double)MSlt4.size(); 				//21. M catalyzes S1, S2, S3
				MSgeq4cat2 	= kcatM2*(double)MSgeq4.size();				//22. M catalyzes S4, S5
				DSlt4cat	= kcatD*(double)DSlt4.size();				//23. D catalyzes S1, S2, S3
				DSgeq4cat	= kcatD*(double)DSgeq4.size();				//24. D catalyzes S4, S5
				
				prop_tot = S0_syn + S0_deg + Slt4_deg + Sgeq4_deg + MSeq0_deg + MSlt4_deg + MSgeq4_deg + DSlt4_deg + DSgeq4_deg + S0M_bind + MS0_diss + Slt4M_bind + MSlt4_diss + Sgeq4M_bind + MSgeq4_diss + Slt4D_bind + DSlt4_diss + Sgeq4D_bind + DSgeq4_diss + MSeq0cat1 + MSlt4cat2 + MSgeq4cat2 + DSlt4cat + DSgeq4cat;
				
				//Next time instance:
				r1 = mtrand();
				delta = -log(r1)/prop_tot;
				time = time + delta;
		
		
				//Choose the reaction
				r2 = mtrand();
				r2 = r2*prop_tot;

				reacted=0;
				
				//1.
				prop_this=prop_tot - S0_syn;
				if(r2>prop_this&&reacted==0)
				{	Seq0.add(x);
					reacted = 1;
					
				}
				
				//2.
				prop_this=prop_this - S0_deg;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)Seq0.size());
					Seq0.erase(nthS);
					reacted = 1;
	
				}

				//3.
				prop_this=prop_this - Slt4_deg;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)Slt4.size());
					Slt4.erase(nthS);
					reacted = 1;
					
				}
				
				//4.
				prop_this=prop_this - Sgeq4_deg;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)Sgeq4.size());
					Sgeq4.erase(nthS);
					reacted = 1;
				}
				
				//5.
				prop_this=prop_this - MSeq0_deg;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)MSeq0.size());
					MSeq0.erase(nthS);
					nM++;
					reacted = 1;
				
				}
				
				//6.
				prop_this=prop_this - MSlt4_deg;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)MSlt4.size());
					MSlt4.erase(nthS);
					nM++;
					reacted = 1;
				}
				
				//7.
				prop_this=prop_this - MSgeq4_deg;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)MSgeq4.size());
					MSgeq4.erase(nthS);
					nM++;
					reacted = 1;
				}
				
				//8.
				prop_this=prop_this - DSlt4_deg;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)DSlt4.size());
					DSlt4.erase(nthS);
					nD++;
					reacted = 1;
				}
				
				//9.
				prop_this=prop_this - DSgeq4_deg;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)DSgeq4.size());
					DSgeq4.erase(nthS);
					nD++;
					reacted = 1;
				}
				
				//10.
				prop_this=prop_this - S0M_bind;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)Seq0.size());
					int u = Seq0[nthS];
					Seq0.erase(nthS);
					MSeq0.add(u);
					nM--;
					reacted = 1;
				}
				
				//11.
				prop_this=prop_this - MS0_diss;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)MSeq0.size());
					int u = MSeq0[nthS];
					MSeq0.erase(nthS);
					Seq0.add(u);
					nM++;
					reacted = 1;
				}
				
				//12.
				prop_this=prop_this - Slt4M_bind;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)Slt4.size());
					int u = Slt4[nthS];
					Slt4.erase(nthS);
					MSlt4.add(u);
					nM--;
					reacted = 1;
				}
				
				//13.
				prop_this=prop_this - MSlt4_diss;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)MSlt4.size());
					int u = MSlt4[nthS];
					MSlt4.erase(nthS);
					Slt4.add(u);
					nM++;
					reacted = 1;
				}
				
				//14.
				prop_this=prop_this - Sgeq4M_bind;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)Sgeq4.size());
					int u = Sgeq4[nthS];
					Sgeq4.erase(nthS);
					MSgeq4.add(u);
					nM--;
					reacted = 1;
				}
				
				//15.
				prop_this=prop_this - MSgeq4_diss;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)MSgeq4.size());
					int u = MSgeq4[nthS];
					MSgeq4.erase(nthS);
					Sgeq4.add(u);
					nM++;
					reacted = 1;
				}
				
				//16.
				prop_this=prop_this - Slt4D_bind;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)Slt4.size());
					int u = Slt4[nthS];
					Slt4.erase(nthS);
					DSlt4.add(u);
					nD--;
					reacted = 1;
				}
				
				//17.
				prop_this=prop_this - DSlt4_diss;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)DSlt4.size());
					int u = DSlt4[nthS];
					DSlt4.erase(nthS);
					Slt4.add(u);
					nD++;
					reacted = 1;
				}
				
				//18.
				prop_this=prop_this - Sgeq4D_bind;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)Sgeq4.size());
					int u = Sgeq4[nthS];
					Sgeq4.erase(nthS);
					DSgeq4.add(u);
					nD--;
					reacted = 1;
				}
				
				//19.
				prop_this=prop_this - DSgeq4_diss;
				if(r2>prop_this&&reacted==0)
				{	nthS = (int)(mtrand()*(double)DSgeq4.size());
					int u = DSgeq4[nthS];
					DSgeq4.erase(nthS);
					Sgeq4.add(u);
					nD++;
					reacted = 1;
				}
				
				//20.
				prop_this=prop_this - MSeq0cat1;
				if(r2>prop_this&&reacted==0)
				{   nthS = (int)(mtrand()*(double)MSeq0.size());
					int u = MSeq0[nthS];
					MSeq0.erase(nthS);
					Slt4.add(u+1);
					nM++;
					reacted = 1;
				}
				
				//21.
				prop_this=prop_this - MSlt4cat2;
				if(r2>prop_this&&reacted==0)
				{   nthS = (int)(mtrand()*(double)MSlt4.size());
					int u = MSlt4[nthS];
					MSlt4.erase(nthS);
					if (u==1 || u ==2)Slt4.add(u+1);
					else Sgeq4.add(u+1);
					nM++;
					reacted = 1;
				}
				
				//22.
				prop_this=prop_this - MSgeq4cat2;
				if(r2>prop_this&&reacted==0)
				{   nthS = (int)(mtrand()*(double)MSgeq4.size());
					int u = MSgeq4[nthS];
					MSgeq4.erase(nthS);
					Sgeq4.add(u+1);
					nM++;
					reacted = 1;
				}
				
				//23.
				prop_this=prop_this - DSlt4cat;
				if(r2>prop_this&&reacted==0)
				{   nthS = (int)(mtrand()*(double)DSlt4.size());
					int u = DSlt4[nthS];
					DSlt4.erase(nthS);
					Seq0.add(x);
					nD++;
					reacted = 1;
				}
				
				//24.
				prop_this=prop_this - DSgeq4cat;
				if(r2>prop_this&&reacted==0)
				{   nthS = (int)(mtrand()*(double)DSgeq4.size());
					int u = DSgeq4[nthS];
					DSgeq4.erase(nthS);
					Seq0.add(x);
					nD++;
					reacted = 1;
				}
				
				if (time>tmax&&(time-delta)<tmax)							
				{	avg_Stot = avg_Stot + (double)Seq0.size() + (double)MSeq0.size() + (double)Slt4.size() + (double)Sgeq4.size() + (double)MSlt4.size() + (double)MSgeq4.size() + (double)DSlt4.size() + (double)DSgeq4.size();
					avg_Sprime = avg_Sprime +  ((double)Sgeq4.size() + (double)MSgeq4.size() + (double)DSgeq4.size())/((double)Seq0.size() + (double)MSeq0.size() + (double)Slt4.size() + (double)Sgeq4.size() + (double)MSlt4.size() + (double)MSgeq4.size() + (double)DSlt4.size() + (double)DSgeq4.size());
				}
			}
		}
			
		avg_Stot = avg_Stot/(double)Nsim;
		avg_Sprime = avg_Sprime/(double)Nsim;
		cout<<avg_Sprime<<"\t"<<avg_Stot<<endl;
		outFile<<avg_Sprime<<"\t"<<avg_Stot<<endl;
		Minit = Minit*10.0;
	}
}