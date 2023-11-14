//NANOROD: mc.h MC Class (Revision Date: Oct 27, 2023)
#ifndef _MC_H
#define _MC_H
#include "system.h"
#include "utils.h"
#include "molecule.h"
#include "hbond.h"
#include "energies.h"
#include "grid.h"
//#include "energies.h"
class MC
{
    public:
        System S;
        Energy E;
        double energy, time;
        int nbr_g=27;//number of number grids
        //new vector for molecule
        vector<Molecule> Mnew;
        MC(){energy=0.0; time=0.0; }
        void WriteTemplate();
        void LogProfile(int, double );
        void Sweep();
        double MoveMolecule();
        bool Glauber(double, double);
        bool MC::Arrhenius(double delta, double rand);
        double TotalEnergy();
};
#endif

