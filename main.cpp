//
//  main.cpp
//  Stream_HPC
//
//  Created by changjiang GOU on 04/04/2020.
//  Copyright © 2020 Changjiang GOU. All rights reserved.
//

#include <iostream>
#include <string>
#include <getopt.h>
#include "utils.hpp"
#include "algorithms.hpp"
#include "common.hpp"

vector<double> BandWidth;
double period;
unsigned int Config_Core;

int main(int argc, const char * argv[]) {
//-----------------------------------
    
//    string inputfile = "apps";
//    string outdir = "../applications_newformat/";
//    ConvertFormat(inputfile.c_str(), outdir);

//    -----------------------------------

//    string inputfile = "apps";
//    CheckGraphs(inputfile.c_str());
    
//    -----------------------------------
    
//    string partition_heuristic;
//    string map_heuristic;
//
//    string result_file ="../result/general/"+partition_heuristic;
//    result_file+="_";
//    result_file+=map_heuristic;
//    result_file+="_CCR201.txt";
    
//    int opt;
//    while ((opt = getopt(argc, argv, "gbsmtr")) != -1) {
//        switch (opt) {
//            case 'g':
//                partition_heuristic = "GroupCell" ;
//                break;
//            case 'b':
//                partition_heuristic = "BreakForkJoin(Dynamic)";
//                break;
//            case 's':
//                partition_heuristic = "SaveProc(Dynamic)";
//                break;
//            case 'm':
//                partition_heuristic = "MaxSpeed";
//                break;
//            case 't':
//                map_heuristic = "MapTopology";
//                break;
//            case 'r':
//                map_heuristic = "MapRank";
//                break;
//            default: /* '?' */
//                fprintf(stderr, "Usage: %s [-g] GroupCell [-b] BreakForkJoin(Dynamic)\n",
//                        argv[0]);
//                exit(EXIT_FAILURE);
//        }
//    }
    
//    string inputfile = "period_beta_chains_101";
////    string inputfile = "period_beta_301";
////    string inputfile = "apps_test";
//
////    string result_file = "../result/general/BreakForkJoin(Dynamic)_MapTopology_CCR301.txt";
////    string result_file = "../result/general/SaveProc(Dynamic)_MapTopology_CCR301.txt";
////    string result_file = "../result/general/GroupCell_MapTopology_CCR301.txt";
////    string result_file = "../result/28May/GroupCell_MapRank_CCR101.txt";
////    string result_file = "../result/general/SaveProc(Dynamic)_MapRank_CCR301.txt";
////    string result_file = "../result/28May/BreakForkJoin(Dynamic)_MapRank_CCR301.txt";
////    string result_file = "../result/general/MaxSpeed_CCR301.txt";
//
//    string result_file = "../result/28May/DP_CCR101.txt";
//
//
//
//    string graphName;
//    string line;
//
//    unsigned int nbr_block=1;
//    vector<double> Period;
//    double energy;
//    double bandwidth_1;
//    bool map_success=false;
//    vector<unsigned int> ProcessorsLeft(nbr_block,Config_Core);
//    unordered_map<string, tuple<double, unsigned int, unsigned int>> Map;
//    vector<unsigned int> BrokenEdges;
//    vector<unsigned int> Blocks;
//    vector<unsigned int> Cores;
//    vector<unsigned int> Copies;
//    forward_list<pair<int, int>> Idle_Cores;
//    bool partition_success;
////    vector<list<unsigned int>> TasksMapped(1);
//
//    array<unsigned int, 6> Number_Blocks = {1, 2, 3, 4, 5, 6};
//    unsigned int total_core = 480;
//
//    ofstream OutFile(result_file);
//
//    for (auto nbr_block_iter=Number_Blocks.begin(); nbr_block_iter!=Number_Blocks.end(); ++nbr_block_iter) {
//        nbr_block = *nbr_block_iter;
////        TasksMapped.resize(nbr_block);
//        for (unsigned int nbr_core=1; nbr_core<=total_core/nbr_block; ++nbr_core) {
//            Config_Core = nbr_core;
//
//            ifstream OpenFile(inputfile);
//            while (getline(OpenFile,line)) {
//                BandWidth.clear();
//                Period.clear();
//                stringstream ss(line);
//                ss>>graphName;
//                while (Period.size()<7) {
//                    ss>>period;
//                    Period.push_back(period);
//                }
//                ss>>bandwidth_1;
//                BandWidth.push_back(bandwidth_1);
//                BandWidth.push_back(BandWidth[0]/16);
//
////                cout<<"graph "<<graphName<<", nbr_block "<<nbr_block<<", nbr_core "<<Config_Core<<" ";
//                OutFile<<graphName<<" "<<nbr_block<<" "<<nbr_core<<" ";
//                while (!Period.empty()) {
//                    period = Period.back();
//                    Period.pop_back();
//
//                    Graph* graph = new Graph();
//                    parse_graph(graphName.c_str(), graph);
////                    graph->Print();
//
//                    Map.clear();
//                    Idle_Cores.clear();
//                    BrokenEdges.clear();
//                    Blocks.clear();
//                    Cores.clear();
//                    Copies.clear();
//                    energy = -1;
//
////                    cout<<endl<<period<<" ";
//
////                    -----------------------------------
//                    energy = Dynamic_Get_E(graph, 1, graph->GetNodes()->size(), nbr_block, Config_Core, &Map);
////                    InterpertRegister(&Map, &BrokenEdges, &Blocks, &Cores, &Copies, 1, graph->GetNodes()->size(), nbr_block, Config_Core, true);
//
////                    -----------------------------------
////                    map_success = MaxSpeed(graph, nbr_block, Config_Core, &BrokenEdges);
//
////                    GroupCell(graph, nbr_block, Config_Core, &BrokenEdges);
////                    partition_success = true;
//
////                    partition_success = BreakFJ_DP(graph, nbr_block, Config_Core, &BrokenEdges);
//
////                    SaveProc_Entry(graph, graph->GetSource()->GetId(), graph->GetNodes()->size(), &BrokenEdges);
//
////                    cout<<"broken edges: ";
////                    for (vector<unsigned int>::iterator it=BrokenEdges.begin(); it!=BrokenEdges.end(); ++it) {
////                        cout<<*it<<" ";
////                    }
////                    cout<<endl;
//
////                    -----------------------------------
//
//
////                    if (partition_success == true) {
////                        Graph* Qgraph = new Graph();
////                        BuildQGraph(graph, &BrokenEdges, Qgraph);
//////                        Qgraph->Print();
////                        ProcessorsLeft.clear();
////                        ProcessorsLeft.assign(nbr_block, Config_Core);
////                        map_success = Map_Ranked(&BrokenEdges, &ProcessorsLeft, nbr_block, Qgraph);
////
////                        if(map_success==true){
////                            map_success = CheckMap(Qgraph);
////                        }
////
////                        if (map_success) {
////                            energy = EnergyCost(Qgraph);
////                        }
////
////                        delete Qgraph;
////                    } else {
////                        map_success = false;
////                    }
//
////                    map_success = Map_Topology(graph, &BrokenEdges, &ProcessorsLeft, nbr_block, Qgraph);
//
////                    for (unsigned int i=0; i<nbr_block; ++i) {
////                        TasksMapped[i].clear();
////                    }
////                    map_success = Map_Topology(Qgraph, nbr_block, &TasksMapped);
//
////                    for (vector<Node*>::iterator it=Qgraph->GetNodes()->begin(); it!=Qgraph->GetNodes()->end(); ++it) {
////                        cout<<"part "<<(*it)->GetId()<<" is map onto block "<<(*it)->GetBlock()<<", core "<<(*it)->GetCore()<<endl;
////                    }
//
//                    if (energy>0) {
//                        OutFile<<energy<<" ";
////                        cout<<energy<<endl;
//                    } else {
//                        OutFile<<-1<<" ";
////                        cout<<-1<<endl;
//                    }
//
//                    delete graph;
//                }
//                OutFile<<endl;
//            }
//            OpenFile.close();
//        }
//    }
//
//    OutFile.close();
//    -----------------------------------
    
//    string inputfile = "apps";
//    Stastics(inputfile.c_str());
    
//    -----------------------------------
    
//    string inputfile = "period_beta_chains_101";
//        string inputfile = "period_beta_101";
        string inputfile = "apps_test";
    
    //    string result_file = "../result/general/BreakForkJoin(Dynamic)_MapTopology_CCR301.txt";
    //    string result_file = "../result/general/SaveProc(Dynamic)_MapTopology_CCR301.txt";
    //    string result_file = "../result/general/GroupCell_MapTopology_CCR301.txt";
        string result_file = "../result/28May/Test_GroupCell_MapRank_CCR301.txt";
    //    string result_file = "../result/general/SaveProc(Dynamic)_MapRank_CCR301.txt";
//        string result_file = "../result/28May/Complet_BreakForkJoin(Dynamic)_MapRank_CCR101.txt";
//        string result_file = "../result/28May/Complet_MaxSpeed_CCR101.txt";
    
//    string result_file = "../result/28May/DP_CCR101.txt";
    
    
    
    string graphName;
    string line;
    
    unsigned int nbr_block=1;
    vector<double> Period;
    double energy;
    double bandwidth_1;
    bool map_success=false;
    vector<unsigned int> ProcessorsLeft(nbr_block,Config_Core);
    unordered_map<string, tuple<double, unsigned int, unsigned int>> Map;
    vector<unsigned int> BrokenEdges;
    vector<unsigned int> Blocks;
    vector<unsigned int> Cores;
    vector<unsigned int> Copies;
    forward_list<pair<int, int>> Idle_Cores;
    bool partition_success;
    //    vector<list<unsigned int>> TasksMapped(1);
    
    ofstream OutFile(result_file);
    
    nbr_block = 4;
    Config_Core = 128;
        //        TasksMapped.resize(nbr_block);
            
    ifstream OpenFile(inputfile);
    while (getline(OpenFile,line)) {
        BandWidth.clear();
        Period.clear();
        stringstream ss(line);
        ss>>graphName;
        while (Period.size()<7) {
            ss>>period;
            Period.push_back(period);
        }
        ss>>bandwidth_1;
        BandWidth.push_back(bandwidth_1);
        BandWidth.push_back(BandWidth[0]/16);
                
                //                cout<<"graph "<<graphName<<", nbr_block "<<nbr_block<<", nbr_core "<<Config_Core<<" ";
        OutFile<<graphName<<" "<<nbr_block<<" "<<Config_Core<<" ";
        while (!Period.empty()) {
            period = Period.back();
            Period.pop_back();
            
            Graph* graph = new Graph();
            parse_graph(graphName.c_str(), graph);
                    //                    graph->Print();
                    
            Map.clear();
            Idle_Cores.clear();
            BrokenEdges.clear();
            Blocks.clear();
            Cores.clear();
            Copies.clear();
            energy = -1;
                    
                    //                    cout<<endl<<period<<" ";
                    
                    //                    -----------------------------------
//            energy = Dynamic_Get_E(graph, 1, graph->GetNodes()->size(), nbr_block, Config_Core, &Map);
                    //                    InterpertRegister(&Map, &BrokenEdges, &Blocks, &Cores, &Copies, 1, graph->GetNodes()->size(), nbr_block, Config_Core, true);
                    
                    //                    -----------------------------------
//            map_success = MaxSpeed(graph, nbr_block, Config_Core, &BrokenEdges);
//            partition_success = true;
            
            nbr_block = 2;
            Config_Core = 8;
            GroupCell(graph, nbr_block, Config_Core, &BrokenEdges);
            partition_success = true;
            
//            partition_success = BreakFJ_DP(graph, nbr_block, Config_Core, &BrokenEdges);
            
                    //                    SaveProc_Entry(graph, graph->GetSource()->GetId(), graph->GetNodes()->size(), &BrokenEdges);
                    
                    //                    cout<<"broken edges: ";
                    //                    for (vector<unsigned int>::iterator it=BrokenEdges.begin(); it!=BrokenEdges.end(); ++it) {
                    //                        cout<<*it<<" ";
                    //                    }
                    //                    cout<<endl;
                    
                    //                    -----------------------------------
                    
                    
            if (partition_success == true) {
                Graph* Qgraph = new Graph();
                BuildQGraph(graph, &BrokenEdges, Qgraph);
                    //                        Qgraph->Print();
                ProcessorsLeft.clear();
                ProcessorsLeft.assign(nbr_block, Config_Core);
                map_success = Map_Ranked(&BrokenEdges, &ProcessorsLeft, nbr_block, Qgraph);
                    
                if(map_success==true){
                    map_success = CheckMap(Qgraph);
                }
                    
                if (map_success) {
                    energy = EnergyCost(Qgraph);
                }
                    delete Qgraph;
                } else {
                    map_success = false;
                }
            
                    //                    map_success = Map_Topology(graph, &BrokenEdges, &ProcessorsLeft, nbr_block, Qgraph);
                    
                    //                    for (unsigned int i=0; i<nbr_block; ++i) {
                    //                        TasksMapped[i].clear();
                    //                    }
                    //                    map_success = Map_Topology(Qgraph, nbr_block, &TasksMapped);
                    
                    //                    for (vector<Node*>::iterator it=Qgraph->GetNodes()->begin(); it!=Qgraph->GetNodes()->end(); ++it) {
                    //                        cout<<"part "<<(*it)->GetId()<<" is map onto block "<<(*it)->GetBlock()<<", core "<<(*it)->GetCore()<<endl;
                    //                    }
                    
                if (energy>0) {
                    OutFile<<energy<<" ";
                        //                        cout<<energy<<endl;
                } else {
                    OutFile<<-1<<" ";
                        //                        cout<<-1<<endl;
                }
                    delete graph;
                }
                OutFile<<endl;
            }
            OpenFile.close();
    
    OutFile.close();
    
    
    
}
