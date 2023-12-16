//NANOROD: system.cpp System Class Function Definitions (Revision Date: Oct 27, 2023)
#include "system.h"
void System::Create()
{
    
    cout<<"Creating System"<<endl;
    int i,j,k,l,n;
    //Allocate Grid
    Grid g;
    for(k=0; k<NGRID; k++)
    {
        for(j=0; j<NGRID; j++)
        {
          for(i=0; i<NGRID; i++)
          {
              g.g_index=GridIndex_index(i,j,k,NGRID);
              g.nbr.clear();
              g.plist.clear();
              g.cm.x=(double(i)+0.5)*GRIDL-0.5*L;
              g.cm.y=(double(j)+0.5)*GRIDL-0.5*L;
              g.cm.z=(double(k)+0.5)*GRIDL-0.5*L;
              for (int inei=i-1;inei<i+2;inei++)
              {
                for (int jnei=j-1;jnei<j+2;jnei++)
                {
                    for (int knei=k-1;knei<k+2;knei++)
                    {
                        int indexn=GridIndex_index(inei,jnei,knei,NGRID);
                        g.nbr.push_back(indexn);
                        
                    }
                }
              }
              g.n=0;//We fill this later now
              G.push_back(g);
             
          }
        }
    }
   
    //Fill particles Use random permutation
    bool flag;
    H.clear();//init hbond list
    int NPair=NMOL/2;
    
    int count=0;
    /*
    for(int i=0;i<NPair;i++)
    {
        Molecule m1;//initial test for 2NMOL bonded particles
        m1.MOL_ID=2*i;
        m1.AID=i;
        m1.AsubID=0;
        m1.centre=XYZ(gsl_rng_uniform(gsl_r)*L-0.5*L,gsl_rng_uniform(gsl_r)*L-0.5*L,gsl_rng_uniform(gsl_r)*L-0.5*L);
        m1.gID=GridIndex_xyz(m1.centre,NGRID,GRIDL,L);
        G[m1.gID].n+=1;
        G[m1.gID].plist.push_back(m1.MOL_ID);
        m1.UpdateVertices();
        m1.nbonds=1;
        m1.hbond_list.push_back(hbond(2*i,2*i+1,0,3));
        m1.vertype[0]='I';
        M.push_back(m1);
        Molecule m2;
        m2.MOL_ID=2*i+1;
        m2.AID=i;
        m2.AsubID=1;
        m2.centre=image((XYZ(2.26,0,1.06)+m1.centre),L);
        m2.gID=GridIndex_xyz(m2.centre,NGRID,GRIDL,L);
        G[m2.gID].n+=1;
        G[m2.gID].plist.push_back(m2.MOL_ID);
        m2.orientation=angle_to_quarternion(M_PI/3,0,0);
        m2.UpdateVertices();
        m2.nbonds=1;
        m2.hbond_list.push_back(hbond(2*i+1,2*i,3,0));
        m2.vertype[3]='I';
        M.push_back(m2);
        if(min_d2(m1.ver[0],m2.ver[3],L)>0.58*0.58)
        {
            cout<<"wrong initialization"<<endl;
            exit(3);
        }
        Aggregate A;
        A.n=2;
        A.rg2=2.2;
        A.M_A.push_back(m1.MOL_ID);
        A.M_A.push_back(m2.MOL_ID);
        A.L=L;
        Ag.push_back(A);
        H.push_back(hbond(2*i,2*i+1,0,3));

    }*/
    for(int i=0; i<NMOL; i++)
    {
        Molecule m;
        m.MOL_ID=i;
        m.AID=i;//at first each molecule is one aggregate
        m.AsubID=0;//every molecule is the only element of the aggregate
        m.centre=XYZ(gsl_rng_uniform(gsl_r)*L-0.5*L,gsl_rng_uniform(gsl_r)*L-0.5*L,gsl_rng_uniform(gsl_r)*L-0.5*L);
        m.gID=GridIndex_xyz(m.centre,NGRID,GRIDL,L);
        G[m.gID].n+=1;
        G[m.gID].plist.push_back(m.MOL_ID);
        m.orientation=RandomRotate(angle_to_quarternion(0,0,0), gsl_rng_uniform(gsl_r),gsl_rng_uniform(gsl_r),gsl_rng_uniform(gsl_r));
        m.UpdateVertices();
        M.push_back(m);
        Aggregate A;
        A.n=1;
        A.rg2=1.1;
        A.M_A.push_back(m.MOL_ID);
        A.L=L;
        Ag.push_back(A);

    }
    //init a ring

}

void System::WriteMol2(int timestep)
{
    ofstream out;
    out.open("conf.mol2",ios::app);
    out<<timestep<<endl;
    int n_ver=M[0].N_VER;
    //define filename from variable
    //char FileName[100];

    //sprintf(FileName,"_NG_%d_l0_%lf_lr_%d_t_%lf_f_%lf_e_%lf_C_%lf.mol2",NGRID,lambda0,lambdar,theta, tau_off, eta_off, cofactor);
    out<<"@<TRIPOS>MOLECULE"<<endl;
    out<<"Nanorod"<<endl;
    out<<NMOL*(n_ver+1)<<"\t"<<NMOL*n_ver+H.size()<<endl;
    out<<"SMALL"<<endl;
    out<<"NO_CHARGES"<<endl;

    out<<"@<TRIPOS>ATOM"<<endl;

    string name,type;
    
    int count=0;
    
    XYZ im,shift;
    
    for(int i=0; i<NMOL; i++)
    {
        shift=image(M[i].centre,L)-M[i].centre;
        im=M[i].centre+shift;
        out<<setw(6)<<++count<<"\t"<<"1"<<"\t"<<setw(8)<<im.x<<"\t"<<setw(8)<<im.y<<"\t"<<setw(8)<<im.z<<"\t"<<"N.3"<<endl;
        
        for(int j=0; j<n_ver; j++)
        {
            im=M[i].ver[j]+shift;
            out<<setw(6)<<++count<<"\t"<<"2"<<"\t"<<setw(8)<<im.x<<"\t"<<setw(8)<<im.y<<"\t"<<setw(8)<<im.z<<"\t"<<"C.1"<<endl;
        }
    }
          
    out<<"@<TRIPOS>BOND"<<endl;
    
    count=0;
    for(int i=0; i<NMOL; i++)
    {
        for(int j=0; j<n_ver; j++)
        {
            out<<setw(8)<<++count<<"\t"<<setw(8)<<(n_ver+1)*i+1<<"\t"<<setw(8)<<(n_ver+1)*i+j+2<<"\t"<<setw(2)<<"1"<<endl;
        }
    }
    list<hbond>::iterator it;
    for(it=H.begin();it!=H.end();it++)
    {
        hbond writebond=*it;
        out<<setw(8)<<++count<<"\t"<<setw(8)<<(n_ver+1)*writebond.M1+writebond.arm1+2<<"\t"<<setw(8)<<(n_ver+1)*writebond.M2+writebond.arm2+2<<"\t"<<setw(2)<<"1"<<endl;
            
        
        
    }
    out<<"@<TRIPOS>FF_PBC"<<endl;
    out<<setw(8)<<"v1.0"<<"\t"<<setw(8)<<"1"<<"\t"<<setw(8)<<-L/2<<"\t"<<setw(8)<<-L/2<<"\t"<<setw(8)<<-L/2<<"\t"<<setw(8)<<L/2<<"\t"<<setw(8)<<L/2<<"\t"<<setw(8)<<L/2<<endl;   
    

    out.close();
    //   cout<<"Mol2 input written in\t"<<FileName<<endl;
    return;
}
void System::WriteDump(int timestep)
{
    char FileName[200];
    sprintf(FileName,"%s_Dump.lammpstrj",Description.c_str());
    ofstream out;
    out.open(FileName,ios::app);
    
    out<<"ITEM: TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<"ITEM: NUMBER OF ATOMS"<<endl;
    int NT=NMOL*7;
    out<<NT<<endl;
    out<<"ITEM: BOX BOUNDS pp pp pp"<<endl;
    out<<-0.5*L<<"\t"<<0.5*L<<endl;
    out<<-0.5*L<<"\t"<<0.5*L<<endl;
    out<<-0.5*L<<"\t"<<0.5*L<<endl;
    out<<"ITEM: ATOMS index type x y z"<<endl;
    /*int NT=NMOL*7;
    out<<NT<<endl;
    out<<"Time="<<timestep<<endl;*/
    XYZ im;
    XYZ im_centre;
    for(int i=0;i<NMOL;i++)
    {   
        im_centre=image(M[i].centre,L);
        out<<setw(6)<<i*7+1<<"\t"<<1<<"\t"<<setw(8)<<im_centre.x<<"\t"<<setw(8)<<im_centre.y<<"\t"<<setw(8)<<im_centre.z<<endl;
        for(int j=0;j<6;j++)
        {
            im=image(M[i].ver[j],L);
            out<<setw(6)<<i*7+j+2<<"\t"<<2<<"\t"<<setw(8)<<im.x<<"\t"<<setw(8)<<im.y<<"\t"<<setw(8)<<im.z<<endl;
        }
    }
    out.close();
}
void System::WriteBond(int timestep)
{
    ofstream out;
    out.open("Bondlist.txt",ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<setw(12)<<"molecule1"<<"\t"<<setw(12)<<"molecule2"<<"\t"<<setw(12)<<"arm1"<<"\t"<<setw(8)<<"arm2"<<endl;
    list<hbond>::iterator it;
    for(it=H.begin();it!=H.end();it++)
    {
        hbond writebond=*it;
        out<<setw(12)<<writebond.M1<<setw(12)<<writebond.M2<<setw(12)<<writebond.arm1<<setw(12)<<writebond.arm2<<endl;
            
        
        
    }
    out.close();
}
void System::WriteOrientation(int timestep)
{
    ofstream out;
    out.open("Orientation.txt",ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<setw(12)<<"molecule1"<<"\t"<<setw(12)<<"molecule2"<<"\t"<<setw(12)<<"arm1"<<"\t"<<setw(8)<<"arm2"<<endl;
    for(int j=0;j<NMOL;j++)
    {
        
        out<<setw(12)<<M[j].orientation.w<<setw(12)<<M[j].orientation.x<<setw(12)<<M[j].orientation.y<<setw(12)<<M[j].orientation.z<<endl;
            
        
    
    }
    out.close();

}
void System::WriteGrid(int timestep)
{
    ofstream out;
    out.open("Grid.txt",ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    for(int i=0;i<NGRID3;i++)
    {
     out<<setw(12)<<G[i].g_index<<setw(12)<<G[i].cm.x<<setw(12)<<G[i].cm.y<<setw(12)<<G[i].cm.z<<setw(12)<<G[i].n<<setw(8)<<endl;
              out<<"neighborlist"<<endl;
              for(int l=0;l<G[i].nbr.size();l++)
              {
                out<<G[i].nbr[l]<<endl;
              }
              
    }
    out.close();
}
void System::WriteAggregate(int timestep)
{
    ofstream out;
    out.open("Aggregate.txt",ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<setw(12)<<"AID"<<setw(12)<<"moleculeid"<<endl;
    for(int i=0;i<Ag.size();i++)
    {
        Aggregate ag=Ag[i];
        for(int j=0;j<ag.n;j++)
        out<<setw(12)<<i<<setw(12)<<ag.cm.x<<setw(12)<<ag.cm.y<<setw(12)<<ag.cm.z<<setw(12)<<ag.M_A[j]<<endl;
            
        
        
    }
    out.close();
}
void System::UpdateGrid()
{
    
    int i,j,k,l,n;
    //Allocate Grid
    int g_index;
    for(k=0; k<NGRID; k++)
    {
        for(j=0; j<NGRID; j++)
        {
          for(i=0; i<NGRID; i++)
          {
              
              g_index=GridIndex_index(i,j,k,NGRID);
              
             
              G[g_index].cm.x=(double(i)+0.5)*GRIDL-0.5*L;
              G[g_index].cm.y=(double(j)+0.5)*GRIDL-0.5*L;
              G[g_index].cm.z=(double(k)+0.5)*GRIDL-0.5*L;
              for (int inei=0;inei<3;inei++)
              {
                for (int jnei=0;jnei<3;jnei++)
                {
                    for (int knei=0;knei<3;knei++)
                    {
                        int indexn=GridIndex_index(inei+i-1,jnei+j-1,knei+k-1,NGRID);
                        G[g_index].nbr[inei*9+jnei*3+knei]=indexn;
                        
                    }
                }
              }
              
             
          }
        }
    }
}


