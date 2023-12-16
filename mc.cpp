//NANOROD: mc.cpp MC Class Function Definitions (Revision Date: October 27,2023)
#include "mc.h" 
void MC::Sweep()
{
    Mnew=S.M;
    E.L=S.L;
    S.WriteDump(0);
    double accept=0.0;
   
    int nsample=S.nsweep/100;
    if(nsample<1)
      nsample=1;
    
    energy=TotalEnergy();
   
    
    WriteTemplate();
	LogProfile(0,accept);
    
    for(int i=1; i<=S.nsweep; i++)
    {
        time+=S.deltat;
        accept+=MoveMolecule();
        
        if(i%nsample==0)
        {
            
            LogProfile(i,accept);
            WriteEnergy(i);
            S.WriteMol2(i);
            S.WriteDump(i);
            S.WriteBond(i);
            S.WriteGrid(i);
            S.WriteAggregate(i);
            accept=0.0;
        }
    }
	
   
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
    if(S.M.size()<S.NMOL)
    {
        cout<<"wrongsize"<<endl;
    }
    //S.UpdateGrid();
    //check if the total particle in the celllists add to N
    int total_particles=0;

    for(int i=0;i<S.NGRID3;i++)
    {
        /*if(S.G[i].nbr.size()!=nbr_g)
        {
            cout<<"Neighbor Number wrong"<<endl;
            exit(0);
        }*/
        if(S.G[i].n!=S.G[i].plist.size())
        {
            cout<<"Number wrong"<<endl;
            exit(0);
        }
        total_particles+=S.G[i].n;
    }
    if(total_particles!=S.NMOL)
    {
            cout<<"Total number wrong"<<endl;
            exit(0);
    }
    ofstream out;
    out.open("events.txt",ios::app);
    ofstream out2;
    out2.open("Error.txt",ios::app);
    //ofstream out2;
    //out2.open("bond_energy.txt",ios::app);
    double accept=0.0;//accept events
    
    int N_break=0;
    
    
    /*//break bonds
    list<hbond>::iterator it;
    double rand=gsl_rng_uniform(S.gsl_r);
    if(rand<S.omega_B);//the relative frequency of breaking bond to frequency of diffusion step
    {
        for(it=S.H.begin();it!=S.H.end();)
        {
            hbond old_hbond=*it;
            //calculate bond_dissociation energy
            double E_dis=0;
            int free_bonds=0;
            E_dis+=S.E_1;//the basic enthalpy change of one bond
            //count # of freed bonds
            //find neighbor arms,first the one of M1, then the one of M2
            int neighborarm1=neighborarm(old_hbond.arm1);
            free_bonds+=2;
            int bonded_index1;
            vector<hbond> old_hbondlist=S.M[old_hbond.M1].hbond_list;
            for(int p=0;p<old_hbondlist.size();p++)
            {
                if(old_hbondlist[p].arm1==neighborarm1)
                {
                    free_bonds-=2;
                }
                if(old_hbondlist[p].arm1==old_hbond.arm1)
                {
                    bonded_index1=p;
                }
            } 
            int neighborarm2=neighborarm(old_hbond.arm2);
            free_bonds+=2;
            vector<hbond> bonded_neighbor_hbondlist=S.M[old_hbond.M2].hbond_list;
            int bonded_index2;
            for(int p=0;p<bonded_neighbor_hbondlist.size();p++)
            {
                if(bonded_neighbor_hbondlist[p].arm1==neighborarm2)
                {
                    free_bonds-=2;
                }
                if(bonded_neighbor_hbondlist[p].arm1==old_hbond.arm2)
                {
                    bonded_index2=p;
                }
            }                 
            E_dis+=free_bonds*S.free_bond_freeenergy;
                
            if(Arrhenius(1,E_dis,gsl_rng_uniform(S.gsl_r)))
            {
                //break bond
                out<<"Break bond"<<setw(12)<<old_hbond.M1<<setw(12)<<old_hbond.M2<<setw(12)<<old_hbond.arm1<<setw(12)<<old_hbond.arm2<<endl;
                S.M[old_hbond.M1].vertype[old_hbond.arm1]='A';
                S.M[old_hbond.M1].hbond_list[bonded_index1]=S.M[old_hbond.M1].hbond_list.back();
                S.M[old_hbond.M1].nbonds-=1;
                S.M[old_hbond.M1].hbond_list.pop_back();
                S.M[old_hbond.M2].vertype[old_hbond.arm2]='A';
                S.M[old_hbond.M2].hbond_list[bonded_index2]=S.M[old_hbond.M2].hbond_list.back();
                S.M[old_hbond.M2].nbonds-=1;
                S.M[old_hbond.M2].hbond_list.pop_back();
                it=S.H.erase(it);        
            }
            else{
                ++it;
            }
        }
    }*/
    //diffusion rotation
    //diffusion translation
    //diffusion limited bond formation
    //individual diffusion
    for(int i=0; i<S.NMOL; i++)
    {
        int index=gsl_rng_uniform_int(S.gsl_r,S.NMOL);//randomly pick one molecule
        double delta=0;//energy difference
        //int index=sequence_M[i];//index of the current trial molecule
        vector<hbond> old_hbondlist=S.M[index].hbond_list;
        Molecule newmolecule=S.M[index];//pass old molecule info to new molecule
        newmolecule.centre=RandomTranslate(S.M[index].centre,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));//translate molecule
        //update gid, put into new molecule
        XYZ image_center=image(newmolecule.centre,S.L);
        int new_gID=GridIndex_xyz(image_center,S.NGRID,S.GRIDL,S.L);
        newmolecule.gID=new_gID;
        newmolecule.orientation=RandomRotate(S.M[index].orientation,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));//rotate molecule
        newmolecule.UpdateVertices();
        //define new bondlist and a variable for whether there is new bond
        int new_hbond=-1;
        vector<hbond> new_hbondlist;
        new_hbondlist.clear();
        double r2_newbond=S.L*S.L;
        double delta_fenesum=0;
        if (newmolecule.nbonds>0)   //calculate hbond_energy
        {
            
            for(int n=0; n<newmolecule.nbonds; n++)
            {
                double delta_fene=E.hbonde_fene(newmolecule,S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2)-E.hbonde_fene(S.M[index],S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2);
                double delta_angle=E.hbonde_angle(newmolecule,S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2)-E.hbonde_angle(S.M[index],S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2);
                double delta_dihedral=E.hbonde_dihedral(newmolecule,S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2)-E.hbonde_dihedral(S.M[index],S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2);
                delta+=(delta_fene+delta_angle+delta_dihedral);
                delta_fenesum+=delta_fene;
                
            }
        }
        //add wca energy
        //add new wca
        for(int l=0;l<nbr_g;l++)
        {
            //ID of new neighboring grids
            int temp_gID=S.G[new_gID].nbr[l];
            if(temp_gID>=S.NGRID3)
            {
                cout<<new_gID<<"Error: wrong grid id"<<endl;
                cout<<temp_gID<<endl;
                exit(0);
            }
            if(S.G[temp_gID].plist.empty())
            {
                continue;
            }
            //iterate about molecules in the neighborlist
            list<int>::iterator it;
            for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
            {
                int j=*it;
                if(j!=index)
                {
                    //image distance
                    double rnew2=min_d2(S.M[j].centre,newmolecule.centre,S.L);
                    delta+=E.WCA(rnew2);
                    //find possible new bonds
                    double r2_cm=min_d2(newmolecule.centre,S.M[j].centre,S.L);
                    if (r2_cm<S.search2_cm)//check if the cm distance is in the range, can save time
                    {
                        
                        //iterate about new molecule vertices
                        for(int k=0;k<newmolecule.N_VER;k++)
                        {
                            //iterate about the neighbor vertices
                            for(int l=0;l<S.M[j].N_VER;l++)
                            {
                                double r2_arms=min_d2(newmolecule.ver[k],S.M[j].ver[l],S.L);
                                XYZ arm1=image(newmolecule.ver[k],S.L);
                                XYZ centre1=image(newmolecule.centre,S.L);
                                XYZ neighborarm1=image(newmolecule.ver[neighborarm(k)],S.L);
                                XYZ arm2=image(S.M[j].ver[l],S.L);
                                XYZ centre2=image(S.M[j].centre,S.L);
                                XYZ neighborarm2=image(S.M[j].ver[neighborarm(l)],S.L);
                                
                                double alpha=angle_vectors(real_vector(arm1-centre1,S.L),real_vector(arm2-arm1,S.L));
                                double beta=angle_vectors(real_vector(arm2-centre2,S.L),real_vector(arm1-arm2,S.L));
                                double xhi=dihedral_vectors(real_vector(neighborarm1-centre1,S.L),real_vector(centre2-centre1,S.L),real_vector(neighborarm2-centre2,S.L));
                              
                                //check type, distance between arms, angles and dihedral
                                if (newmolecule.vertype[k]=='A'&&S.M[j].vertype[l]=='A'&&r2_arms<S.searchl2_bond&&alpha<S.search_angle&&beta<S.search_angle&& abs(xhi)<S.search_xhi)
                                {
                                    
                                    
                                    bool bonded=false;//check if there is already bonds between the two molecules, only one bond can exist between two molecules
                                    /*for(int n=0;n<old_hbondlist.size();n++)
                                    {
                                        if(old_hbondlist[n].M2==j)
                                        {
                                            bonded=true;
                                        }
                                    }*/
                                    for(int n=0;n<new_hbondlist.size();n++)
                                    {
                                        if(new_hbondlist[n].M2==j)
                                        {
                                            bonded=true;
                                        }
                                    }
                                    if(bonded==false)
                                    {
                                        new_hbondlist.push_back(hbond(index,j,k,l));
                                        new_hbond=1;
                                    }
                                }
                            }
                        }
                        
                    }


                }
            }
        }
        //subtract old wca
        int old_gID=S.M[index].gID;
        for(int l=0;l<nbr_g;l++)
        {
            //ID of new neighboring grids
            int temp_gID=S.G[old_gID].nbr[l];
            list<int>::iterator it;
            for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
            {
                int j=*it;
                if(j!=index)
                {
                    //image distance
                    double rold2=min_d2(S.M[j].centre,S.M[index].centre,S.L);
                    delta-=E.WCA(rold2);
                    


                }
            }
        }
        if(Glauber(delta,gsl_rng_uniform(S.gsl_r)))
        {
            //Accept move
            if(delta_fenesum>10000)
                {
                    cout<<"toolarge"<<endl;
                }
            S.M[index]=newmolecule;
            S.Ag[newmolecule.AID].M_A[newmolecule.AsubID]=newmolecule.MOL_ID;
            accept+=1.0;
            energy+=delta;

            //update grid list
            if(new_gID!=old_gID)
            {
                S.G[old_gID].n-=1;
                S.G[old_gID].plist.remove(S.M[index].MOL_ID);
                
                S.G[new_gID].n+=1;
                S.G[new_gID].plist.push_back(S.M[index].MOL_ID);
            }
        
            //form bonds if accept
            if (new_hbond!=-1)
            {
                for(int m=0;m<new_hbondlist.size();m++) 
                {
                    //update hbond list and nbonds and vertypes of the two molecules
                    S.M[index].hbond_list.push_back(new_hbondlist[m]);
                    S.M[index].vertype[new_hbondlist[m].arm1]='I';
                    S.M[index].nbonds+=1;
                    S.M[new_hbondlist[m].M2].hbond_list.push_back(hbond(new_hbondlist[m].M2,new_hbondlist[m].M1,new_hbondlist[m].arm2,new_hbondlist[m].arm1));
                    S.M[new_hbondlist[m].M2].vertype[new_hbondlist[m].arm2]='I';
                    S.M[new_hbondlist[m].M2].nbonds+=1;
                    out<<"form a hbond"<<setw(12)<<new_hbondlist[m].M1<<setw(12)<<new_hbondlist[m].M2<<setw(12)<<new_hbondlist[m].arm1<<setw(12)<<new_hbondlist[m].arm2<<endl;
                    //update the hbondlist of the system
                    if(new_hbondlist[m].M1<new_hbondlist[m].M2)
                    {
                        S.H.push_back(new_hbondlist[m]);
                    }
                    else
                    {
                        S.H.push_back(hbond(new_hbondlist[m].M2,new_hbondlist[m].M1,new_hbondlist[m].arm2,new_hbondlist[m].arm1));
                    }
                    //update the aggregate
                    if(S.M[index].AID!=S.M[new_hbondlist[m].M2].AID)
                    {
                        int AID1=S.M[index].AID;
                        int AID2=S.M[new_hbondlist[m].M2].AID;
                        Aggregate old_aggregate1=S.Ag[AID1];
                        Aggregate old_aggregate2=S.Ag[AID2];
                        //merge two aggregate
                        Aggregate new_aggregate=old_aggregate1;
                        new_aggregate.n+=old_aggregate2.n;
                        for(int i=0;i<old_aggregate2.M_A.size();i++)
                        {
                            S.M[old_aggregate2.M_A[i]].AID=AID1;
                            S.M[old_aggregate2.M_A[i]].AsubID+=old_aggregate1.n;
                        }
                        new_aggregate.M_A.insert(new_aggregate.M_A.end(),old_aggregate2.M_A.begin(),old_aggregate2.M_A.end());
                        S.Ag[AID1]=new_aggregate;
                        S.Ag[AID2]=S.Ag.back();
                        for(int i=0;i<S.Ag[AID2].M_A.size();i++)
                        {
                            S.M[S.Ag[AID2].M_A[i]].AID=AID2;    
                        }
                        S.Ag.pop_back();

                    }
                }

            }
        }
        
        
    }
    //Aggregate diffusion
    //vector<Aggregate>::iterator ita;
    for(int a=0;a<S.Ag.size();a++)
    {
        Aggregate old_aggregate=S.Ag[a];
        if(old_aggregate.n<2)
        {
            continue;}
        if(old_aggregate.n>=2)
        {
            
            double deltaag=0;//energy difference, only WCA with other aggregate
            vector<Molecule> newmolecule_list;
            for(int i=0;i<old_aggregate.n;i++)
            {
                newmolecule_list.push_back(S.M[old_aggregate.M_A[i]]);
            }
            XYZ difference=image(S.M[old_aggregate.M_A[1]].centre-S.M[old_aggregate.M_A[0]].centre,S.L);
            Aggregate new_aggregate=old_aggregate;//first give all the information of old aggregate to new aggregate
            //calculate the diffusion step according to n
            double MCstep_ag=S.MCstep/pow(old_aggregate.n,1/6);//new MCstep according to particle number
            double MCstep_agrotate=S.MCstep_rotation/pow(old_aggregate.n,1/2);//new MCstep_rotation
            XYZ total_translate=RandomTranslatestep(MCstep_ag,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));//generate random translation for the whole aggregate
            //cout<<total_translate.x<<total_translate.x<<endl;
            quarternion total_rotate=RandomRotatestep(MCstep_agrotate,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));//random rotation operation
            XYZ cm;//center of mass
            XYZ shift=image(S.M[old_aggregate.M_A[0]].centre,S.L);//shift one particle to the center of box
            shift=XYZ(-shift.x,-shift.y,-shift.z);
            XYZ sum_xyz;
            double sum_mass=7*old_aggregate.n;
            for(int i=1;i<old_aggregate.n;i++)
            {
                XYZ centre_shift=S.M[old_aggregate.M_A[i]].centre+shift;
                XYZ image_centre=image(centre_shift,S.L);
                sum_xyz=sum_xyz+image_centre;
                for(int j=0;j<6;j++)
                {
                    XYZ ver_shift=S.M[old_aggregate.M_A[i]].ver[j]+shift;
                    XYZ image_ver=image(ver_shift,S.L);
                    sum_xyz=sum_xyz+image_ver;
                }
            }
            cm=sum_xyz/sum_mass;//generate center of mass
            //new_aggregate.cm=image(cm-shift,S.L);
            //translate and rotate all particles
            for(int i=0;i<old_aggregate.n;i++)
            {
                XYZ centre_cm_vec=image(S.M[old_aggregate.M_A[i]].centre+shift-cm,S.L);
                XYZ new_centre_cm_vec=quarterrotation(centre_cm_vec,total_rotate);
                XYZ movement=new_centre_cm_vec-centre_cm_vec;
                newmolecule_list[i].centre=newmolecule_list[i].centre+total_translate+movement;//temporarily delete the rotation
                newmolecule_list[i].orientation=quartermulti(total_rotate,newmolecule_list[i].orientation);
                newmolecule_list[i].orientation.normalize();
                XYZ image_center=image(newmolecule_list[i].centre,S.L);
                //update grid id, put in the new aggregate
                int new_gID=GridIndex_xyz(image_center,S.NGRID,S.GRIDL,S.L);
                newmolecule_list[i].gID=new_gID;
                newmolecule_list[i].UpdateVertices();
                /*for(int j=0;j<6;j++)
                {
                    XYZ ver_cm_vec=image(old_aggregate.M_A[i].ver[j]+shift-cm,S.L);
                    XYZ new_ver_cm_vec=quarterrotation(ver_cm_vec,total_rotate);
                    XYZ movement=new_ver_cm_vec-ver_cm_vec;
                    new_aggregate.M_A[i].ver[j]=new_aggregate.M_A[i].ver[j]+total_translate;//temporarily delete the rotation
                }*/
            }
            XYZ difference2=image(newmolecule_list[1].centre-newmolecule_list[0].centre,S.L);
            
            /*if(!Aresame(difference.x,difference2.x))
            {
                cout<<"error in total translate"<<endl;
                exit(2);
            }*/
            //calculate wca energy difference and find new bonds
            int new_hbond=-1;
            vector<hbond> new_hbondlist;
            new_hbondlist.clear();
            for(int i=0;i<old_aggregate.n;i++)
            {
                
                
                Molecule newmolecule=newmolecule_list[i];//pass info of new molecule in new aggregate to newmolecule
                Molecule oldmolecule=S.M[old_aggregate.M_A[i]];//pass info of old molecule in old aggregate to oldmolecule
                vector<hbond> old_hbondlist=newmolecule_list[i].hbond_list;
                
                
                //add wca energy
                //add new wca
                int new_gID=newmolecule_list[i].gID;
                for(int l=0;l<nbr_g;l++)
                {
                    //ID of new neighboring grids
                    int temp_gID=S.G[new_gID].nbr[l];
                    if(temp_gID>=S.NGRID3)
                    {
                        cout<<new_gID<<"Error: wrong grid id"<<endl;
                        cout<<temp_gID<<endl;
                        exit(0);
                    }
                    if(S.G[temp_gID].plist.empty())
                    {
                        continue;
                    }
                    //iterate about molecules in the neighborlist
                    list<int>::iterator it;
                    for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
                    {
                        int j=*it;
                        if(S.M[j].AID!=newmolecule.AID)
                        {
                            //image distance
                            double rnew2=min_d2(S.M[j].centre,newmolecule.centre,S.L);
                            deltaag+=E.WCA(rnew2);
                            //find possible new bonds
                            double r2_cm=min_d2(newmolecule.centre,S.M[j].centre,S.L);
                            if (r2_cm<S.search2_cm)//check if the cm distance is in the range, can save time
                            {
                                
                                //iterate about new molecule vertices
                                for(int k=0;k<newmolecule.N_VER;k++)
                                {
                                    //iterate about the neighbor vertices
                                    for(int l=0;l<S.M[j].N_VER;l++)
                                    {
                                        double r2_arms=min_d2(newmolecule.ver[k],S.M[j].ver[l],S.L);
                                        XYZ arm1=image(newmolecule.ver[k],S.L);
                                        XYZ centre1=image(newmolecule.centre,S.L);
                                        XYZ neighborarm1=image(newmolecule.ver[neighborarm(k)],S.L);
                                        XYZ arm2=image(S.M[j].ver[l],S.L);
                                        XYZ centre2=image(S.M[j].centre,S.L);
                                        XYZ neighborarm2=image(S.M[j].ver[neighborarm(l)],S.L);
                                        
                                        double alpha=angle_vectors(real_vector(arm1-centre1,S.L),real_vector(arm2-arm1,S.L));
                                        double beta=angle_vectors(real_vector(arm2-centre2,S.L),real_vector(arm1-arm2,S.L));
                                        double xhi=dihedral_vectors(real_vector(neighborarm1-centre1,S.L),real_vector(centre2-centre1,S.L),real_vector(neighborarm2-centre2,S.L));
                                        
                                        //check type, distance between arms, angles and dihedral
                                        
                                        if (newmolecule.vertype[k]=='A'&&S.M[j].vertype[l]=='A'&&r2_arms<S.searchl2_bond&&alpha<S.search_angle&&beta<S.search_angle&& abs(xhi)<S.search_xhi)
                                        {
                                            
                                            
                                            bool bonded=false;//check if there is already bonds between the two molecules, only one bond can exist between two molecules
                                            for(int n=0;n<new_hbondlist.size();n++)
                                            {
                                                if(new_hbondlist[n].M2==j)
                                                {
                                                    bonded=true;
                                                }
                                            }
                                            if(bonded==false)
                                            {
                                                new_hbondlist.push_back(hbond(newmolecule.MOL_ID,j,k,l));
                                                new_hbond=1;
                                            }
                                        }
                                    }
                                }
                                
                            }


                        }   
                    }
                }
                //subtract old wca
                //subtract old wca
                int old_gID=S.M[old_aggregate.M_A[i]].gID;
                for(int l=0;l<nbr_g;l++)
                {
                    //ID of new neighboring grids
                    int temp_gID=S.G[old_gID].nbr[l];
                    list<int>::iterator it;
                    for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
                    {
                        int j=*it;
                        if(S.M[j].AID!=oldmolecule.AID)
                        {
                            //image distance
                            double rold2=min_d2(S.M[j].centre,S.M[old_aggregate.M_A[i]].centre,S.L);
                            deltaag-=E.WCA(rold2);
                        }
                    }
                }       
            }
            if(Glauber(deltaag,gsl_rng_uniform(S.gsl_r)))
            {
                //Accept move
                //update aggregate
                S.Ag[a]=new_aggregate;
                for(int i=0;i<new_aggregate.n;i++)
                {
                    int old_gID=S.M[new_aggregate.M_A[i]].gID;
                    int new_gID=newmolecule_list[i].gID;
                    if(old_gID!=new_gID)
                    {
                        //update grid list
                        
                        S.G[old_gID].n-=1;
                        S.G[old_gID].plist.remove(newmolecule_list[i].MOL_ID);
                    
                        S.G[new_gID].n+=1;
                        S.G[new_gID].plist.push_back(newmolecule_list[i].MOL_ID);

                    }
                    S.M[new_aggregate.M_A[i]]=newmolecule_list[i];
                    
                }
                accept+=1.0;
                energy+=deltaag;
                //form bonds if accept
                if (new_hbond!=-1)
                {
                    for(int m=0;m<new_hbondlist.size();m++) 
                    {
                        //update hbond list and nbonds and vertypes of the two molecules
                        S.M[new_hbondlist[m].M1].hbond_list.push_back(new_hbondlist[m]);
                        S.M[new_hbondlist[m].M1].vertype[new_hbondlist[m].arm1]='I';
                        S.M[new_hbondlist[m].M1].nbonds+=1;
                        S.M[new_hbondlist[m].M2].hbond_list.push_back(hbond(new_hbondlist[m].M2,new_hbondlist[m].M1,new_hbondlist[m].arm2,new_hbondlist[m].arm1));
                        S.M[new_hbondlist[m].M2].vertype[new_hbondlist[m].arm2]='I';
                        S.M[new_hbondlist[m].M2].nbonds+=1;
                        out<<"form a hbond"<<setw(12)<<new_hbondlist[m].M1<<setw(12)<<new_hbondlist[m].M2<<setw(12)<<new_hbondlist[m].arm1<<setw(12)<<new_hbondlist[m].arm2<<endl;
                        //update the hbondlist of the system
                        if(new_hbondlist[m].M1<new_hbondlist[m].M2)
                        {
                            S.H.push_back(new_hbondlist[m]);
                        }
                        else
                        {
                            S.H.push_back(hbond(new_hbondlist[m].M2,new_hbondlist[m].M1,new_hbondlist[m].arm2,new_hbondlist[m].arm1));
                        }
                        //update the aggregate
                        if(S.M[new_hbondlist[m].M1].AID!=S.M[new_hbondlist[m].M2].AID)
                        {
                            int AID1=S.M[new_hbondlist[m].M1].AID;
                            int AID2=S.M[new_hbondlist[m].M2].AID;
                            Aggregate old_aggregate1=S.Ag[AID1];
                            Aggregate old_aggregate2=S.Ag[AID2];
                            //merge two aggregate
                            Aggregate new_aggregate=old_aggregate1;
                            new_aggregate.n+=old_aggregate2.n;
                            for(int i=0;i<old_aggregate2.M_A.size();i++)
                            {
                                S.M[old_aggregate2.M_A[i]].AID=AID1;
                                S.M[old_aggregate2.M_A[i]].AsubID+=old_aggregate1.n;
                            }
                            new_aggregate.M_A.insert(new_aggregate.M_A.end(),old_aggregate2.M_A.begin(),old_aggregate2.M_A.end());
                            S.Ag[AID1]=new_aggregate;
                            S.Ag[AID2]=S.Ag.back();
                            for(int i=0;i<S.Ag[AID2].M_A.size();i++)
                            {
                                S.M[S.Ag[AID2].M_A[i]].AID=AID2;    
                            }
                            S.Ag.pop_back();

                        }
                    }

                }
            } 
                     
        }
              
    }        
    
    
    
    
    return accept/double(S.NMOL);
}

bool MC::Glauber(double delta, double rand)
{
    if (delta>100000)
        return false;
    else
    {
        if(1.0/(exp(delta)+1.0)>rand)
            return true;
        else
            return false;
    }
    
}
bool MC::Arrhenius(double A,double delta, double rand)
{
    
    if(A*(exp(-delta))>rand)
        {
        return true;}
    else
        return false;
    
    
}
double MC::WCAEnergy()
{
    double ewca=0;
    list<int>::iterator it;
    for(int i=0;i<S.NMOL;i++)
    {
        int old_gID=S.M[i].gID;
        for(int l=0;l<nbr_g;l++)
        {
            //ID of new neighboring grids
            int temp_gID=S.G[old_gID].nbr[l];
            
            for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
            {
                int j=*it;
                if(j>i)//avoid double counting
                {
                    //image distance
                    double r2=min_d2(S.M[j].centre,S.M[i].centre,S.L);
                    ewca+=E.WCA(r2);
                    


                }
            }
        }
    }
    return ewca;
}
double MC::FENE_energy()
{
    double efene=0;
    list<hbond>::iterator it;
    for(it=S.H.begin();it!=S.H.end();it++)
    {
        hbond j=*it;
        efene+=E.hbonde_fene(S.M[j.M1],S.M[j.M2],j.arm1,j.arm2);
    }
    return efene;
}
double MC::Angle_energy()
{
    double eangle=0;
    list<hbond>::iterator it;
    for(it=S.H.begin();it!=S.H.end();it++)
    {
        hbond j=*it;
        eangle+=E.hbonde_angle(S.M[j.M1],S.M[j.M2],j.arm1,j.arm2);
    }
    return eangle;
}
double MC::Dihedral_energy()
{
    double edihedral=0;
    list<hbond>::iterator it;
    for(it=S.H.begin();it!=S.H.end();it++)
    {
        hbond j=*it;
        edihedral+=E.hbonde_dihedral(S.M[j.M1],S.M[j.M2],j.arm1,j.arm2);
    }
    return edihedral;   
}
int MC::bond_energy()
{
    double ebond=0;
    int nbond=S.H.size();
    //ebond=-nbond*S.E_1;
    return nbond;
}
int MC::bond_freeze_freenerngy()
{
    double efreeze=0;
    int nfreeze=0;
    for(int i=0;i<S.NMOL;i++)
    {
        Molecule m=S.M[i];
        for(int l=0;l<m.N_VER;l++)
        {
            if(m.vertype[l]=='I'||m.vertype[neighborarm(l)]=='I')
            nfreeze+=1;

        }
    }
    efreeze=nfreeze*S.free_bond_freeenergy;
    return nfreeze;
}
double MC::TotalEnergy()
{
    double totalenergy=0.0;
    /*for(int i;i<S.NMOL;i++)
    {
        for(int j=0;j<S.M[i].nbonds;j++)
        {
            totalenergy+=E.hbonde(S.M[i],S.M[S.M[i].hbond_list[j].M2],j);
            //totalenergy+=E.LJ(S.M[i],S.M[j]);
        }
        
    }*/
    return totalenergy;
}
double MC::WriteEnergy(int timestep)
{
    ofstream out;
    out.open("energy.txt",ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<"WCA"<<'\t'<<"FENE"<<'\t'<<"Angle"<<'\t'<<"Dihedral"<<'\t'<<"Bond"<<'\t'<<"Bond_freeze"<<endl;
    out<<WCAEnergy()<<'\t'<<FENE_energy()<<'\t'<<Angle_energy()<<'\t'<<Dihedral_energy()<<'\t'<<bond_energy()<<'\t'<<bond_freeze_freenerngy()<<endl;
    out.close();
}
