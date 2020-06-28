//
//  algorithms.hpp
//  Stream_HPC
//
//  Created by changjiang GOU on 05/04/2020.
//  Copyright © 2020 Changjiang GOU. All rights reserved.
//

#ifndef algorithms_hpp
#define algorithms_hpp

#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <tuple>
#include <limits>
#include <forward_list>
#include <list>
#include <stack>
#include <queue>
#include <iostream>
#include <fstream>
#include "graph.hpp"
#include "utils.hpp"
#include "common.hpp"

double Dynamic_Get_E(Graph* graph, unsigned int task_left, unsigned int task_right, int block, int core, unordered_map<string, tuple<double,unsigned int, unsigned int>>* Register);

void InterpertRegister(unordered_map<string, tuple<double, unsigned int, unsigned int>>* Register, vector<unsigned int>* BrokenEdges, vector<unsigned int>* Blocks, vector<unsigned int>* Cores, vector<unsigned int>* Copies, unsigned int left_task, unsigned int right_task, unsigned int block, unsigned int core, bool Print);

void Traversal(Graph* graph, list<tuple<unsigned int, unsigned int, double, double>>* List_MaxS, list<tuple<unsigned int, unsigned int, double, double>>* List_Trip);

void GroupCell(Graph* graph, int block, int core, vector<unsigned int>* brokenEs);

//void BreakParaComb(Graph* graph, int block, int core, vector<unsigned int>* BrokenEdges);

bool BreakFJ_DP(Graph* graph, const int nbr_block, const int nbr_core, vector<unsigned int>* BrokenEdges);

void SaveProc(Graph* graph, unsigned int task_left, unsigned int task_right, vector<unsigned int>* brokenEs);

void SaveProc_Entry(Graph* graph, unsigned int task_left, unsigned int task_right, vector<unsigned int>* brokenEs);

double LeftPart(Graph* graph, unsigned int task_left, unsigned int task_right, int block, int core, forward_list<pair<int, int>> IdleCores, vector<unsigned int>* brokenEdges, vector<unsigned int>* blocks, vector<unsigned int>* cores, vector<unsigned int>* copies);

void BuildQGraph(Graph* graph, vector<unsigned int>* BrokenEs, Graph* Qgraph);

//bool Map_Topology(Graph* graph, vector<unsigned int>* BrokenEs, vector<unsigned int>* ProcessorsLeft, int block, Graph* Qgraph);

//bool Map_Topology(Graph* Qgraph, const int block, vector<list<unsigned int>>* Sets);

bool Map_Ranked(vector<unsigned int>* BrokenEs, vector<unsigned int>* ProcessorsLeft, int block, Graph* Qgraph);

void DFS(Graph* graph, Node* task_i, vector<unsigned int>* DFSVector);

void ConvertFormat(const char * inputfile, string outdir);

void CheckGraphs(const char * inputfile);

double WhichStructure(Graph* graph, unsigned int task_left, unsigned int task_right, char* structure);

bool MaxSpeed(Graph* graph, unsigned int nbr_blocks, unsigned int nbr_cores, vector<unsigned int>* brokenEs);

void Stastics(const char * inputfile);

double EnergyCost(Graph* Qgraph);

bool CheckMap(Graph* Qgraph);

#endif /* algorithms_hpp */
