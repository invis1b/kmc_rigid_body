//NANOROD: mc.cpp MC Class Function Definitions (Revision Date: October 27,2023)
#include "mc.h" 
void MC::Sweep()
{
    Mnew=S.M;
    
    double accept=0.0;
   
    int nsample=S.nsweep/100;
    if(nsample<1)
      nsample=1;
    
    energy=TotalEnergy();
    cout<<"Initial Energy="<<energy<<endl;
    
    WriteTemplate();
	LogProfile(0,accept);
    
    for(int i=1; i<=S.nsweep; i++)
    {
        time+=S.deltat;
        accept+=MoveMolecule();
        
        if(i%nsample==0)
        {
            accept/=double(nsample);
            LogProfile(i,accept);
            S.WriteMol2(i);
            S.WriteDump(i);
            accept=0.0;
        }
    }
	
   cout<<"Final Energy="<<energy<<endl;
   cout<<"Final Energy by recalculation="<<TotalEnergy()<<endl;
   S.WriteMol2(S.nsweep);
   S.WriteDump(S.nsweep);
}

void MC::WriteTemplate()
{
    ofstream out;
    out.open("_MC.log");
    out<<setw(12)<<"sweep"<<"\t"<<setw(12)<<"time"<<"\t"<<setw(12)<<"energy"<<"\t"<<setw(8)<<"accept"<<endl;
    out.close();
}

void MC::LogProfile(int i, double accept)
{
    ofstream out;
    out.open("_MC.log", ios::app);
    out<<setw(12)<<i<<"\t"<<setw(12)<<time<<"\t"<<setw(12)<<energy<<"\t"<<setw(8)<<accept<<endl;
    out.close();
}


//Returns acceptance fraction
double MC::MoveMolecule()
{
    
    double accept=0.0;
    vector<int> sequence_M;
    sequence_M=generateRandom(S.NMOL);//randomly generate a sequence of the N molecules
    for(int i=0; i<S.NMOL; i++)
    {
        //int i=gsl_rng_uniform_int(S.gsl_r,S.NMOL);
        index=sequence_M[i];//index of the current trial molecule
        Molecule newmolecule=S.M[index];
        newmolecule.centre=RandomTranslate(S.M[index].centre,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        newmolecule.orientation=RandomRotate(S.M[index].orientation,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        newmolecule.UpdateVertices();
        int new_hbond=-1;
        vec<hbond> new_hbondlist;
	    double r2_newbond=S.L*S.L;
        for (int j=0;j<S.NMOL;j++)//loop over all molecules;lately can be improved by celllist
        {   
            if (j!=index)
            {
                double r2_cm=min_d2(S.M[index],S.M[j],S.L);
                if (r2_cm<S.max2_cm)//check if the cm distance is in the range, can save time
                {
                    for(int k=0;k<S.M[index].N_VER;k++)
                    {
                        for(int l=0;l<S.M[j].N_VER;l++)
                        {
                            double r2_arms=min_d2(S.M[index].ver[k],S.M[j].ver[l]);
                            if (r2_arms<S.maxl2_bond&&S.M[index].vertype[k]=='A'&&S.M[j].vertype[l]=='A')
                            {
                                
                                new_hbondlist.push_back(h_bond(index,j,k,l));
                                new_hbond=1;
                            }
                        }
                    }
                    
                }
            }
        }
        if (new_hbond!=-1)
        {
            for(int m=0;m<new_hbondlist.size();m++) 
            S.M[index].hbond_list.push_back(new_hbondlist[m])
            S.M[index].vertype[new_hbondlist[m].arm1]='I'
            S.M[new_hbondlist[m].M2].hbond_list.push_back(h_bond(new_hbondlist[m].M2,new_hbondlist[m].M1,new_hbondlist[m].arm2,new_hbondlist[m].arm1))
            S.M[new_hbondlist[m].M2].vertype[new_hbondlist[m].arm2]='I'
        }
        S.M[index]=newmolecule;
    }
    return accept/double(S.NMOL);
}

bool MC::Glauber(double delta, double rand)
{
    if(1.0/(exp(delta)+1.0)>rand)
      return true;
    else
      return false;
}

double MC::TotalEnergy()
{
    double totalenergy=0.0;
    double r2;
    return totalenergy;
}
