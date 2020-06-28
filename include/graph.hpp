//
//  graph.hpp
//  Stream_HPC
//
//  Created by changjiang GOU on 04/04/2020.
//  Copyright © 2020 Changjiang GOU. All rights reserved.
//

#ifndef graph_hpp
#define graph_hpp

#include <stdio.h>
#include <iostream>
#include <array>
#include "node.hpp"

class Graph {
protected:
    vector<Node*> * nodes;
    unsigned int source_id;
    unsigned int nbr_nodes;
    //unsigned int sink_id;
    
public:
    Graph(){
        source_id = 1;
        nodes = new vector<Node*>();
    }
    
    Graph(unsigned int source, unsigned int sink, unsigned int size){
        source_id = source;
        //sink_id = sink;
        nodes = new vector<Node*>();
        this->AllocateNodes(size);
    }
    
    ~Graph(){
        for (vector<Node*>::iterator it=nodes->begin(); it!=nodes->end(); ++it) {
            delete *it;
        }
        delete nodes;
    }
    
    void SetSize(unsigned int size){
        nbr_nodes = size;
    }
    
    void AllocateNodes(unsigned int size){
        this->SetSize(size);
        nodes->resize(size);
        
        unsigned int i = 0;
        for (vector<Node*>::iterator iter=nodes->begin(); iter!=nodes->end(); iter++) {
            *iter = new Node();
            (*iter)->SetId(i++);
        }
    }
    
    Node * GetSource(){
        return nodes->at(source_id-1);
    }
    
    Node * GetNode(unsigned int node_id){
        return nodes->at(node_id-1);
    }
    
    vector<Node*> * GetNodes(){
        return nodes;
    }
    
    void SetSource(unsigned int sce_id){
        source_id = sce_id;
    }
    
    void AddNode(){
        this->nodes->push_back(new Node());
        this->nodes->back()->SetId(nodes->size());
    }
    
    void MergeNode2Predeccsor(Node* node){
        //assume node has only one predecessor, and its not a fork
        Node* first_predecessor = this->GetNode(node->GetPredecessors()->front());
        vector<unsigned int>* successors;
        vector<unsigned int>* predecessors;
        
        successors = first_predecessor->GetSuccessors();
        successors->assign(node->GetSuccessors()->begin(),node->GetSuccessors()->end());
        
        for (vector<unsigned int>::iterator it=node->GetSuccessors()->begin(); it!=node->GetSuccessors()->end(); ++it) {
            predecessors = this->GetNode(*it)->GetPredecessors();
            for (vector<unsigned int>::iterator iter=predecessors->begin(); iter!=predecessors->end(); ++iter) {
                if(*iter==node->GetId()){
                    *iter = first_predecessor->GetId();
                }
            }
        }
        
        first_predecessor->SetWeight(first_predecessor->GetWeight()+node->GetWeight());
    }
    
    
    void Print(){
        cout<<"id predecessor_id edge_weight node_weight"<<endl;
        vector<double>::iterator edge_iter;
        for (unsigned int i = 1; i<=nbr_nodes; ++i) {
            if (i==source_id) {
                cout<<this->GetNode(i)->GetId()<<" 0 "<<this->GetNode(i)->GetInputEdgeSize()->front()<<" "<<this->GetNode(i)->GetWeight()<<endl;
            }
            edge_iter = this->GetNode(i)->GetInputEdgeSize()->begin();
            for (vector<unsigned int>::iterator iter=this->GetNode(i)->GetPredecessors()->begin(); iter!=this->GetNode(i)->GetPredecessors()->end(); ++iter) {
                cout<<this->GetNode(i)->GetId()<<" "<<(*iter)<<" "<<(*edge_iter)<<" "<<this->GetNode(i)->GetWeight()<<endl;
                ++edge_iter;
            }
        }
    }
};

#endif /* graph_hpp */
