//
//  utils.cpp
//  Stream_HPC
//
//  Created by changjiang GOU on 04/04/2020.
//  Copyright © 2020 Changjiang GOU. All rights reserved.
//

#include "utils.hpp"

void parse_graph(const char *filename, Graph* graph){
    ifstream OpenFile(filename);
    bool nodes_cnt_read = false;
    string line;
    string last_line;
    unsigned int nb_of_nodes=0;
    
    unsigned int id;
    unsigned int predecessor;
    double ew, nw;
    
    if (!nodes_cnt_read) {
        /*Get the number of nodes*/
        
        while(getline(OpenFile,line)){
            last_line = line;
        }
        nb_of_nodes = stoi(last_line);
        OpenFile.clear();
        OpenFile.seekg(0,std::ios::beg);
        
        nodes_cnt_read = true;
        
        /*allocate space for nodes*/
        graph->AllocateNodes(nb_of_nodes);
    }
    /*skip the first line*/
    getline(OpenFile,line);
        
    /*parse actual nodes*/
    while(getline(OpenFile,line)){
        stringstream ss(line);
        ss>>id>>predecessor>>ew>>nw;
        
        graph->GetNode(id)->SetId(id);
        graph->GetNode(id)->SetWeight(nw);
        if (predecessor!=0) {
            graph->GetNode(id)->AddPredecessor(predecessor);
        }
        graph->GetNode(id)->AddInputEdgeSize(ew);
            
        if (predecessor!=0) {
            graph->GetNode(predecessor)->AddSuccessor(id);
        }
            
        if (predecessor == 0) {
            graph->SetSource(id);
        }
    }
    
    OpenFile.close();
}
