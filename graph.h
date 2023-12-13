//KMC: aggregate.h Class (Added: Dec 24, 2010)
//MIXING: enegies.h Energy Class (Revision Date: Feb 2, 2012)
#ifndef _GRAPH_H
#define _GRAPH_H
#include "utils.h"
#include "xyz.h"
#include "quarternion.h"
#include "header.h"
class Graph {
private:
    
    

    void DFS(int node, std::vector<bool>& visited, vector<vector<int>>& component) {
        visited[node] = true;
        component.push_back(adjList[node]); 
        for (int neighbour : adjList[node]) {
            if (!visited[neighbour]) {
                DFS(neighbour, visited, component);
            }
        }
    }

public:
    int V; // Number of vertices
    vector<vector<int>> adjList; // Adjacency list
    Graph(int V) : V(V), adjList(V) {}

    void addEdge(int u, int v) {
        adjList[u].push_back(v);
        adjList[v].push_back(u); // Remove this line for a directed graph
    }

    void removeEdge(int u, int v) {
        adjList[u].erase(remove(adjList[u].begin(), adjList[u].end(), v), adjList[u].end());
        adjList[v].erase(remove(adjList[v].begin(), adjList[v].end(), u), adjList[v].end());
    }

    vector<vector<vector<int>>> getConnectedComponents() {
        vector<bool> visited(V, false);
        vector<vector<vector<int>>> connectedComponents;

        for (int i = 0; i < V; ++i) {
            if (!visited[i]) {
                vector<vector<int>> component;
                DFS(i, visited, component);
                connectedComponents.push_back(component);
            }
        }

        return connectedComponents;
    }
};
#endif