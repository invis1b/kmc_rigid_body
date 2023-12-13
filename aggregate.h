//KMC: aggregate.h Class (Added: Dec 24, 2010)
//MIXING: enegies.h Energy Class (Revision Date: Feb 2, 2012)
#ifndef _AGGREGATE_H
#define _AGGREGATE_H
#include "utils.h"
#include "molecule.h"
#include "hbond.h"
#include "xyz.h"
#include "quarternion.h"
class Aggregate
{
private:
    
    

    /*void DFS(int node, std::vector<bool>& visited, vector<Molecule>& component) {
        visited[node] = true;
        component.push_back(M_A[node]); 

        for (int neighbour : M_A[node].adjList) {
            if (!visited[neighbour]) {
                DFS(neighbour, visited, component);
            }
        }
    }*/
public:
    double L;
    int n; //number of particles
    vector<int> M_A; //vector of Molecules
    double rg2; //radius of gyration squared;
    XYZ cm;
    Aggregate(){n=0; rg2=0.0;cm=XYZ(0,0,0);}
    /*vector<vector<Molecule>> getConnectedComponents() {
        vector<bool> visited(n, false);
        vector<vector<Molecule>> connectedComponents;

        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                vector<Molecule> component;
                DFS(i, visited, component);
                connectedComponents.push_back(component);
            }
        }

        return connectedComponents;
    }*/

};
#endif