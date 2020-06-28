//
//  node.hpp
//  Stream_HPC
//
//  Created by changjiang GOU on 04/04/2020.
//  Copyright © 2020 Changjiang GOU. All rights reserved.
//

#ifndef node_hpp
#define node_hpp

#include <stdio.h>
#include <vector>

using namespace std;

class Node {
protected:
    unsigned int id;
    double weight;
    unsigned int nbr_predecessor=0;
    unsigned int nbr_successor=0;
    vector<unsigned int> * predecessors;
    vector<unsigned int> * successors;
    vector<double> * input_edges_size;
    unsigned int block=0;
    unsigned int core=0;
    double speed=0;
    unsigned int copies=0;
    bool inputEdges_broken=false;
    bool visited = false;
    
public:
    Node(){
        id = 0;
        weight = 0;
        nbr_predecessor = 0;
        nbr_successor = 0;
        block = 0;
        core = 0;
        block = 0;
        inputEdges_broken=false;
        visited = false;
        predecessors = new vector<unsigned int> ();
        input_edges_size = new vector<double> ();
        successors = new vector<unsigned int> ();
    }
    
    ~Node(){
        delete predecessors;
        delete input_edges_size;
        delete successors;
    }
    
    void AddPredecessor(unsigned int predecessor_id){
        this->predecessors->push_back(predecessor_id);
        if (predecessor_id!=0) {
            nbr_predecessor = nbr_predecessor+1;
        }
    }
    
    void AddSuccessor(unsigned int successor_id){
        this->successors->push_back(successor_id);
        nbr_successor = nbr_successor+1;
    }
    
    int GetNbrPredecessor(){
        return nbr_predecessor;
    }
    
    unsigned int GetNbrSuccessor(){
        return nbr_successor;
    }
    
    unsigned int GetTheFirstSuccessor(){
        return this->successors->at(0);
    }
    
    vector<unsigned int> * GetSuccessors(){
        return successors;
    }
    
    void SetSuccessors(vector<unsigned int>* new_successors){
        this->successors = new_successors;
    }
    
    vector<unsigned int> * GetPredecessors(){
        return predecessors;
    }
    
    void ClearPredecessors(){
        this->predecessors->clear();
    }
    
    void AddInputEdgeSize(double size){
        this->input_edges_size->push_back(size);
    }
    
    vector<double> * GetInputEdgeSize(){
        return input_edges_size;
    }
    
    void SetId(unsigned int nid){
        id = nid;
    }
    
    unsigned int GetId(){
        return id;
    }
    
    void SetWeight(double w){
        weight = w;
    }
    
    double GetWeight(){
        return weight;
    }
    
    void SetBlock(unsigned int blockVal){
        block = blockVal;
    }
    
    unsigned int GetBlock(){
        return block;
    }
    
    void SetCore(unsigned int coreVal){
        core = coreVal;
    }
    
    unsigned int GetCore(){
        return core;
    }
    
    void SetSpeed(double speedVal){
        speed = speedVal;
    }
    
    double GetSpeed(){
        return speed;
    }
    
//    double GetDuration(){
//        if (speed!=0) {
//            return (double)(weight/speed);
//        } else {return -1;}
//    }
    
    void SetCopies(unsigned int copies){
        this->copies = copies;
    }
    
    unsigned int GetCopies(){
        return copies;
    }
    
    void BreakInputEdges(){
        inputEdges_broken = true;
    }
    
    bool IsInputEdgesBroken(){
        return inputEdges_broken;
    }
    
    void RestoreInputEdges(){
        inputEdges_broken = false;
    }
    
    void SetVisited(){
        visited = true;
    }
    
    void SetUnVisited(){
        visited = false;
    }
    
    bool IsVisited(){
        return visited;
    }
};

#endif /* node_hpp */
