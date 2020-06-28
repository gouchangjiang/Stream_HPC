//
//  algorithms.cpp
//  Stream_HPC
//
//  Created by changjiang GOU on 05/04/2020.
//  Copyright © 2020 Changjiang GOU. All rights reserved.
//

#include "algorithms.hpp"

///Get energy cost of running tasks between i and j, i should not be larger than j
double Dynamic_Get_E_MaxS(Graph* graph, unsigned int task_i, unsigned int task_j, int block, int core){
    double energy=-1;//-1 represents failure;
    double weight_sum=0;
    
    for (unsigned int i=task_i; i<=task_j; i++) {
        weight_sum = weight_sum + graph->GetNode(i)->GetWeight();
    }
    
    if (weight_sum>(period*SMAX)) {
        return energy;
    }
    
    if (core<1 or block<1) {//no core left
        return energy;
    }
    
    energy = Static_Power*period+ConstantC*pow(SMAX, 2)*weight_sum;
    
    return energy;
}

///Get energy cost of triplicating tasks between i and j,
///i should not larger than j
double Dynamic_Get_E_Trip(Graph* graph, unsigned int task_i, unsigned int task_j, int block, int core){
    double energy=-1; //-1 represents failure
    double weight_sum=0;
    
    for (unsigned int i=task_i; i<=task_j; ++i) {
        weight_sum = weight_sum + graph->GetNode(i)->GetWeight();
    }
    
    if (weight_sum>(period*SMAX)) {
        return energy;
    }
    
    if (core<3 or block<1) {//no core left, no block left
        energy = -1;
        return energy;
    }
    
    double speed;
    for (unsigned int i=0; i<Speed_Options; ++i) {
        speed = DiscSpeed[i];
        if (weight_sum/speed<=period) {
            break;
        }
    }
    
//    assert(graph->GetNode(task_j)->GetNbrSuccessor()==1 or graph->GetNode(task_j)->GetNbrSuccessor()==0);
    
    double ouptput_edge_size;
    if (graph->GetNode(task_j)->GetNbrSuccessor()==1) {
        ouptput_edge_size=graph->GetNode(graph->GetNode(task_j)->GetTheFirstSuccessor())->GetInputEdgeSize()->front();
    } else {
        ouptput_edge_size=0;
    }
    
    energy = 3*(Static_Power*period+ConstantC*pow(speed, 2)*weight_sum)+2*BandWidth[1]*ouptput_edge_size;
    
    return energy;
}

///Get energy cost dissipated on the communication
double Dynamic_Get_E_comm(Graph* graph, double alpha, double beta, double edge_size, unsigned int copies){
    double energy=-1;
    if (edge_size>beta*period) {
        return energy;
    }
    
    energy=copies*alpha*edge_size;
    
    return energy;
}

///Energy for the whole application, calculated in a recursive way
double Dynamic_Get_E(Graph* graph, unsigned int task_left, unsigned int task_right, int block, int core, unordered_map<string, tuple<double,unsigned int, unsigned int>>* Register){
    double energy;
    unsigned int broken_edge, best_practice_way=0;
    double energy_term_1, energy_term_2;
    
    if (block<1 or core<1) {//no available block or core
        energy = -1;
        return energy;
    }
    
    string search_key = to_string(task_left);
    search_key += '_';
    search_key += to_string(task_right);
    search_key += '_';
    search_key += to_string(block);
    search_key += '_';
    search_key += to_string(core);
    
    unordered_map<string, tuple<double,unsigned int, unsigned int>>::const_iterator got=Register->find(search_key);
    if (got!=Register->end()) {// this case has been considered already
        energy = get<0>(got->second);
        return energy;
    }
    
    double energy_maxS = Dynamic_Get_E_MaxS(graph, task_left, task_right, block, core);
    double energy_trip = Dynamic_Get_E_Trip(graph, task_left, task_right, block, core);
    energy_term_1 = energy_maxS;
    energy_term_2 = energy_trip;
    
    if (task_left==task_right) {//only one task in this part
        energy=-1;
        if (energy_term_1>0) {
            if (energy_term_2>0) {
                if (energy_term_1<=energy_term_2) {
                    best_practice_way=1;
                    energy = energy_term_1;
                } else {
                    best_practice_way=2;
                    energy = energy_term_2;
                }
            } else {
                best_practice_way=1;
                energy = energy_term_1;
            }
        } else {
            if (energy_term_2>0) {
                best_practice_way=2;
                energy = energy_term_2;
            }
        }

        Register->insert({search_key,make_tuple(energy,task_left,best_practice_way)});
        return energy;
    }
    
    double energy_rest,energy_rest_1;
    double energy_comm;
    double edge_size;
    double energy_term_3=-1,energy_term_4=-1,energy_term_5=-1,energy_term_6=-1;
    double smallest_energy=numeric_limits<double>::max();
    for (unsigned int task_middle=task_left; task_middle<task_right; task_middle++) {
        energy_maxS = Dynamic_Get_E_MaxS(graph, task_middle+1, task_right, block, core);
        energy_rest = Dynamic_Get_E(graph, task_left, task_middle, block, core-1, Register);
        edge_size = graph->GetNode(task_middle+1)->GetInputEdgeSize()->front();
        energy_comm = Dynamic_Get_E_comm(graph, Alpha[0], BandWidth[0], edge_size, 1);
        if (energy_maxS!=-1 and energy_rest!=-1 and energy_comm!=-1) {
            energy_term_3 = energy_maxS+energy_rest+energy_comm;
        } else {
            energy_term_3 = -1;
        }
        
        energy_rest_1 = Dynamic_Get_E(graph, task_left, task_middle, block-1, Config_Core, Register);
        energy_comm = Dynamic_Get_E_comm(graph, Alpha[1], BandWidth[1], edge_size, 1);
        if (energy_maxS!=-1 and energy_rest_1!=-1 and energy_comm!=-1) {
            energy_term_4 = energy_maxS+energy_rest_1+energy_comm;
        } else {
            energy_term_4 = -1;
        }
        
        energy_trip = Dynamic_Get_E_Trip(graph, task_middle+1, task_right, block, core);
        energy_rest = Dynamic_Get_E(graph, task_left, task_middle, block, core-3, Register);
        energy_comm = Dynamic_Get_E_comm(graph, Alpha[0], BandWidth[0], edge_size, 3);
        if (energy_trip!=-1 and energy_rest!=-1 and energy_comm!=-1) {
            energy_term_5 = energy_trip+energy_rest+energy_comm;
        } else {
            energy_term_5 = -1;
        }
        
//        energy_rest = Dynamic_Get_E(graph, task_left, task_middle, block-1, Config_Core, Register);
        energy_comm = Dynamic_Get_E_comm(graph, Alpha[1], BandWidth[1], edge_size, 3);
        if (energy_trip!=-1 and energy_rest_1!=-1 and energy_comm!=-1) {
            energy_term_6 = energy_trip+energy_rest_1+energy_comm;
        } else {
            energy_term_6 = -1;
        }
        
        if (energy_term_3<smallest_energy and energy_term_3>0) {
            smallest_energy = energy_term_3;
            best_practice_way = 3;
            broken_edge = task_middle+1;// edge broken is task task_middle+1 's input edge
        }
        if (energy_term_4<smallest_energy and energy_term_4>0) {
            smallest_energy = energy_term_4;
            best_practice_way = 4;
            broken_edge = task_middle+1;// edge broken is task task_middle+1 's input edge
        }
        if (energy_term_5<smallest_energy and energy_term_5>0) {
            smallest_energy = energy_term_5;
            best_practice_way = 5;
            broken_edge = task_middle+1;// edge broken is task task_middle+1 's input edge
        }
        if (energy_term_6<smallest_energy and energy_term_6>0) {
            smallest_energy = energy_term_6;
            best_practice_way = 6;
            broken_edge = task_middle+1;// edge broken is task task_middle+1 's input edge
        }
    }
    
    if (energy_term_1 < smallest_energy and energy_term_1>0) {
        smallest_energy = energy_term_1;
        best_practice_way = 1;
        broken_edge = 0; // no edge is broken
    }
    
    if (energy_term_2 < smallest_energy and energy_term_2>0) {
        smallest_energy = energy_term_2;
        best_practice_way = 2;
        broken_edge = 0; // no edge is broken
    }
    
    if (smallest_energy==numeric_limits<double>::max()) {
        smallest_energy = -1; //invalid
        best_practice_way = 0;
        broken_edge = 0; // no edge is broken
    }
    
    //cout<<search_key<<", "<<smallest_energy<<endl;
    Register->insert({search_key,make_tuple(smallest_energy,broken_edge,best_practice_way)});
    
    return smallest_energy;
}

void InterpertRegister(unordered_map<string, tuple<double, unsigned int, unsigned int>>* Register, vector<unsigned int>* BrokenEdges, vector<unsigned int>* Blocks, vector<unsigned int>* Cores, vector<unsigned int>* Copies, unsigned int left_task, unsigned int right_task, unsigned int block, unsigned int core, bool Print){
    unsigned int brokenEdge=1;
    double energy;
    unsigned int best_practice_way;
    unsigned int current_block=block, current_core=core;
    
    string search_key = to_string(left_task);
    search_key += '_';
    search_key += to_string(right_task);
    search_key += '_';
    search_key += to_string(block);
    search_key += '_';
    search_key += to_string(core);
    
    unordered_map<string, tuple<double,unsigned int, unsigned int>>::const_iterator got;
    while (brokenEdge!=0) {
        got = Register->find(search_key);
        if (got!=Register->end()) {
            energy = get<0>(got->second);
            brokenEdge = get<1>(got->second);
            best_practice_way = get<2>(got->second);
            
            Blocks->push_back(current_block);
            Cores->push_back(current_core);
            BrokenEdges->push_back(brokenEdge);
            
            switch (best_practice_way) {
                case 1:
                    Copies->push_back(1);
                    break;
                    
                case 2:
                    Copies->push_back(3);
                    break;
                    
                case 3:
                    Copies->push_back(1);
                    current_core = current_core - 1;
                    break;
                    
                case 4:
                    --current_block;
                    current_core = Config_Core;
                    Copies->push_back(1);
                    break;
                    
                case 5:
                    current_core = current_core - 3;
                    Copies->push_back(3);
                    break;
                    
                case 6:
                    --current_block;
                    current_core = Config_Core;
                    Copies->push_back(3);
                    break;
                    
                default:
                    break;
            }
        } else {
            break;
        }
        
        search_key = to_string(left_task);
        search_key += '_';
        search_key += to_string(brokenEdge-1);
        search_key += '_';
        search_key += to_string(current_block);
        search_key += '_';
        search_key += to_string(current_core);
    }
    
    if (Print==true) {
        cout<<"broken_edge, block, core, copies"<<endl;
        for (unsigned int i=0; i<BrokenEdges->size(); ++i) {
            cout<<BrokenEdges->at(i)<<" "<<Blocks->at(i)<<" "<<Cores->at(i)<<" "<<Copies->at(i)<<endl;
        }
    }
}

//double RightPart(Graph* graph, unsigned int task_left, unsigned int task_right, int block, int core, double alpha_val, double beta_val, forward_list<pair<int, int>>* IdleCores, int* num_core_used, int* id_block_used, int* id_core_used, bool* triplicated){
//    int current_block_1 = block, current_core_1 = core;
//    int current_block_2 = block, current_core_2 = core;
//    double alpha = alpha_val, beta = beta_val;
//    double edge_size, energy_term_1=-1, energy_term_2=-1, energy_comm=-1;
//    edge_size = graph->GetNode(task_left)->GetInputEdgeSize()->front();
//    if (current_block_1<1 or current_core_1<1) {
//        if (!IdleCores->empty()) {
//            current_block_1 = IdleCores->front().first;
//            current_core_1 = IdleCores->front().second;
//            alpha = Alpha[1];
//            beta = BandWidth[1];
//        }
//    }
//
//    if (current_block_1>=1 and current_core_1>=1) {
//        energy_term_1 = Dynamic_Get_E_MaxS(graph, task_left, task_right, current_block_1, current_core_1);
//        energy_comm = Dynamic_Get_E_comm(graph, alpha, beta, edge_size, 1);
//        if (energy_term_1>=0 and energy_comm>=0) {
//            energy_term_1 = energy_term_1+energy_comm;
//        }
//    }
//
//    forward_list<pair<int, int>>::iterator position_to_delete = IdleCores->before_begin();
//    auto iter=IdleCores->begin();
//    if (current_block_2<1 or current_core_2<3) {
//        for (; iter!=IdleCores->end(); ++iter, ++position_to_delete) {
//            if (iter->second>=3) {
//                current_block_2 = iter->first;
//                current_core_2 = iter->second;
//                alpha = Alpha[1];
//                beta = BandWidth[1];
//                break;
//            }
//        }
//    }
//
//    if (current_block_2>=1 and current_core_2>=3) {
//        energy_term_2 = Dynamic_Get_E_Trip(graph, task_left, task_right, current_block_2, current_core_2);
//        energy_comm = Dynamic_Get_E_comm(graph, alpha, beta, edge_size, 3);
//        if (energy_term_2>=0 and energy_comm>=0) {
//            energy_term_2 = energy_term_2+energy_comm;
//        }
//    }
//
//    *num_core_used = 0;
//    if (energy_term_1<0 and energy_term_2<0) {
//        *id_block_used = -1;
//        *id_core_used = -1;
//        return -1;
//    }
//
//    *id_block_used = block;
//    *id_core_used = core;
//
//    if (energy_term_1<0 and energy_term_2>=0) {
//        if (block<1 or core<3) {
//            if (iter->second>3) {
//                iter->second = iter->second -3;
//            } else {
//                IdleCores->erase_after(position_to_delete);
//            }
//            *num_core_used = 0;
//            *id_block_used = current_block_2;
//            *id_core_used = current_core_2;
//        } else {
//            *num_core_used = 3;
//        }
//        *triplicated = true;
//        return energy_term_2;
//    }
//
//    if (energy_term_1>=0 and energy_term_2<0) {
//        if (block<1 or core<1) {
//            if (IdleCores->front().second>1) {
//                IdleCores->front().second -= 1;
//            } else {
//                IdleCores->pop_front();
//            }
//            *num_core_used = 0;
//            *id_block_used = current_block_1;
//            *id_core_used = current_core_1;
//        } else {
//            *num_core_used = 1;
//        }
//        *triplicated = false;
//        return energy_term_1;
//    }
//
//    if (energy_term_1>=0 and energy_term_2>=0) {
//        if (energy_term_1<=energy_term_2) {
//            if (block<1 or core<1) {
//                //cout<<IdleCores->front().first<<" "<<IdleCores->front().second<<endl;
//                if (IdleCores->front().second>1) {
//                    IdleCores->front().second -= 1;
//                } else {
//                    IdleCores->pop_front();
//                }
//                *num_core_used = 0;
//                *id_block_used = current_block_1;
//                *id_core_used = current_core_1;
//            } else {
//                *num_core_used = 1;
//            }
//            *triplicated = false;
//            return energy_term_1;
//        } else {
//            if (block<1 or core<3) {
//                if (iter->second>3) {
//                    iter->second = iter->second -3;
//                } else {
//                    IdleCores->erase_after(position_to_delete);
//                }
//                *num_core_used = 0;
//                *id_block_used = current_block_2;
//                *id_core_used = current_core_2;
//            } else {
//                *num_core_used = 3;
//            }
//            *triplicated = true;
//            return energy_term_2;
//        }
//    }
//
//    return -1;
//}

double RightPart(Graph* graph, unsigned int task_left, unsigned int task_right, int block, int core, double alpha_val, double beta_val, forward_list<pair<int, int>>* IdleCores, int* num_core_used, int* id_block_used, int* id_core_used, bool* triplicated){
    int current_block_1 = block, current_core_1 = core;
    int current_block_2 = block, current_core_2 = core;
    double alpha = alpha_val, beta = beta_val;
    double edge_size, energy_term_1=-1, energy_term_2=-1, energy_comm=-1;
    edge_size = graph->GetNode(task_left)->GetInputEdgeSize()->front();

    forward_list<pair<int, int>>::iterator position_to_delete = IdleCores->before_begin();
    auto iter=IdleCores->begin();
    if (current_block_1==0) {
        if (!IdleCores->empty()) {
            current_block_1 = IdleCores->front().first;
            current_core_1 = IdleCores->front().second;
            alpha = Alpha[1];
            beta = BandWidth[1];
        }
        
        for (; iter!=IdleCores->end(); ++iter, ++position_to_delete) {
            if (iter->second>=3) {
                current_block_2 = iter->first;
                current_core_2 = iter->second;
                alpha = Alpha[1];
                beta = BandWidth[1];
                break;
            }
        }
    }
    
    energy_term_1 = Dynamic_Get_E_MaxS(graph, task_left, task_right, current_block_1, current_core_1);
    energy_comm = Dynamic_Get_E_comm(graph, alpha, beta, edge_size, 1);
    if (energy_term_1>=0 and energy_comm>=0) {
        energy_term_1 = energy_term_1+energy_comm;
    }
    
    energy_term_2 = Dynamic_Get_E_Trip(graph, task_left, task_right, current_block_2, current_core_2);
    energy_comm = Dynamic_Get_E_comm(graph, alpha, beta, edge_size, 3);
    if (energy_term_2>=0 and energy_comm>=0) {
        energy_term_2 = energy_term_2+energy_comm;
    }
    
    *num_core_used = 0;
    if (energy_term_1<0 and energy_term_2<0) {
        *id_block_used = -1;
        *id_core_used = -1;
        return -1;
    }
    
    *id_block_used = block;
    *id_core_used = core;
    if (energy_term_1<0 and energy_term_2>=0) {
        if (block==0) {
            if (iter->second>3) {
                iter->second = iter->second -3;
            } else {
                IdleCores->erase_after(position_to_delete);
            }
            *num_core_used = 0;
            *id_block_used = current_block_2;
            *id_core_used = current_core_2;
        } else {
            *num_core_used = 3;
        }
        *triplicated = true;
        return energy_term_2;
    }
    
    if (energy_term_1>=0 and energy_term_2<0) {
        if (block==0) {
            if (IdleCores->front().second>1) {
                IdleCores->front().second -= 1;
            } else {
                IdleCores->pop_front();
            }
            *num_core_used = 0;
            *id_block_used = current_block_1;
            *id_core_used = current_core_1;
        } else {
            *num_core_used = 1;
        }
        *triplicated = false;
        return energy_term_1;
    }
    
    if (energy_term_1>=0 and energy_term_2>=0) {
        if (energy_term_1<=energy_term_2) {
            if (block==0) {
                if (IdleCores->front().second>1) {
                    IdleCores->front().second -= 1;
                } else {
                    IdleCores->pop_front();
                }
                *num_core_used = 0;
                *id_block_used = current_block_1;
                *id_core_used = current_core_1;
            } else {
                *num_core_used = 1;
            }
            *triplicated = false;
            return energy_term_1;
        } else {
            if (block==0) {
                if (iter->second>3) {
                    iter->second = iter->second -3;
                } else {
                    IdleCores->erase_after(position_to_delete);
                }
                *num_core_used = 0;
                *id_block_used = current_block_2;
                *id_core_used = current_core_2;
            } else {
                *num_core_used = 3;
            }
            *triplicated = true;
            return energy_term_2;
        }
    }
    
    return -1;
}

double LeftPart(Graph* graph, unsigned int task_left, unsigned int task_right, int block, const int core, forward_list<pair<int, int>> IdleCores, vector<unsigned int>* brokenEdges, vector<unsigned int>* blocks, vector<unsigned int>* cores, vector<unsigned int>* copies){

    forward_list<pair<int, int>> original_IdleCores_copy;
    forward_list<pair<int, int>> original_IdleCores_copy1;
    forward_list<pair<int, int>> original_IdleCores_copy2;
    int num_core_used_0, num_core_used_1, num_core_used_2;
    int id_block_used_0, id_block_used_1, id_block_used_2;
    int id_core_used_0, id_core_used_1, id_core_used_2;
    vector<unsigned int> brokenEdges_1;
    vector<unsigned int> brokenEdges_2;
    vector<unsigned int> blocks_1;
    vector<unsigned int> blocks_2;
    vector<unsigned int> cores_1;
    vector<unsigned int> cores_2;
    vector<unsigned int> copies_1;
    vector<unsigned int> copies_2;
    bool triplicated_0, triplicated_1, triplicated_2;
    bool best_case_triplicated_temp, best_case_triplicated;
    unsigned int best_block_temp, best_core_temp;
    unsigned int best_block, best_core;
    
    vector<unsigned int> best_broken_Edges_temp;
    vector<unsigned int> best_blocks_temp;
    vector<unsigned int> best_cores_temp;
    vector<unsigned int> best_copies_temp;
    
    vector<unsigned int> best_broken_Edges;
    vector<unsigned int> best_blocks;
    vector<unsigned int> best_cores;
    vector<unsigned int> best_copies;
    
    int break_node = -1; //input edge is broken
    double energy_left_1, energy_right_1, energy_left_2, energy_right_2;
    
    original_IdleCores_copy.assign(IdleCores.begin(), IdleCores.end());
    double energy_toge = RightPart(graph, task_left, task_right, block, core, Alpha[0], BandWidth[0], &original_IdleCores_copy, &num_core_used_0, &id_block_used_0, &id_core_used_0, &triplicated_0);
    
    double min_E = energy_toge, min_E_temp=-1;
    best_case_triplicated_temp = triplicated_0;
    
    for (unsigned int task_middle = task_left+1; task_middle <= task_right; ++task_middle) {
//        cout<<"task middle "<<task_middle<<", task left "<<task_left<<", task right "<<task_right<<endl;
        original_IdleCores_copy1.assign(IdleCores.begin(), IdleCores.end());
        if (block<=0 and original_IdleCores_copy1.empty()) {
            energy_left_1 = -1;
            energy_right_1 = -1;
        } else if (block>=0){
            energy_right_1 = RightPart(graph, task_middle, task_right, block, core, Alpha[0], BandWidth[0], &original_IdleCores_copy1, &num_core_used_1, &id_block_used_1, &id_core_used_1, &triplicated_1);
            
            if (energy_right_1 <=0 ) {
                energy_left_1 = -1;
            } else {
                energy_left_1 = LeftPart(graph, task_left, task_middle-1, block, core-num_core_used_1, original_IdleCores_copy1, &brokenEdges_1, &blocks_1, &cores_1, &copies_1);
            }
        }
        
        original_IdleCores_copy2.assign(IdleCores.begin(), IdleCores.end());
        energy_right_2 = RightPart(graph, task_middle, task_right, block, core, Alpha[1], BandWidth[1], &original_IdleCores_copy2, &num_core_used_2, &id_block_used_2, &id_core_used_2, &triplicated_2);
        if (core > num_core_used_2 and block > 0) {//some cores on block left unused
            original_IdleCores_copy2.push_front(make_pair(block, core-num_core_used_2));
        }
        
        if (block<=1 and original_IdleCores_copy2.empty()) {
            energy_left_2 = -1;
        } else if (energy_right_2 <= 0){
            energy_left_2 = -1;
        } else if(block>=1) {
            energy_left_2 = LeftPart(graph, task_left, task_middle-1, block-1, Config_Core, original_IdleCores_copy2, &brokenEdges_2, &blocks_2, &cores_2, &copies_2);
        }
        
        if (energy_left_1>=0 and energy_right_1>=0) {
            if (energy_left_2>=0 and energy_right_2>=0) {
//                double temp = energy_left_1+energy_right_1;
//                double temp_1 = energy_left_2+energy_right_2;
//                cout<<"task left, middle, right: "<<task_left<<", "<<task_middle<<", "<<task_right<<", energy_1 "<<temp<<", energy_2 "<<temp_1<<endl;
                if ((energy_left_1+energy_right_1)<=(energy_left_2+energy_right_2)) {
                    min_E_temp = energy_left_1+energy_right_1;
                    best_broken_Edges_temp.assign(brokenEdges_1.begin(), brokenEdges_1.end());
                    best_blocks_temp.assign(blocks_1.begin(), blocks_1.end());
                    best_cores_temp.assign(cores_1.begin(), cores_1.end());
                    best_copies_temp.assign(copies_1.begin(), copies_1.end());
                    best_block_temp = id_block_used_1;
                    best_core_temp = id_core_used_1;
                    best_case_triplicated_temp = triplicated_1;
                } else {
                    min_E_temp = energy_left_2+energy_right_2;
                    best_broken_Edges_temp.assign(brokenEdges_2.begin(), brokenEdges_2.end());
                    best_blocks_temp.assign(blocks_2.begin(), blocks_2.end());
                    best_cores_temp.assign(cores_2.begin(), cores_2.end());
                    best_copies_temp.assign(copies_2.begin(), copies_2.end());
                    best_block_temp = id_block_used_2;
                    best_core_temp = id_core_used_2;
                    best_case_triplicated_temp = triplicated_2;
                }
            } else {
//                double temp_3 = energy_left_1+energy_right_1;
//                cout<<"task left, middle, right: "<<task_left<<", "<<task_middle<<", "<<task_right<<", energy_1 "<<temp_3<<endl;
                min_E_temp = energy_left_1+energy_right_1;
                best_broken_Edges_temp.assign(brokenEdges_1.begin(), brokenEdges_1.end());
                best_blocks_temp.assign(blocks_1.begin(), blocks_1.end());
                best_cores_temp.assign(cores_1.begin(), cores_1.end());
                best_copies_temp.assign(copies_1.begin(), copies_1.end());
                best_block_temp = id_block_used_1;
                best_core_temp = id_core_used_1;
                best_case_triplicated_temp = triplicated_1;
            }
        } else {
            if (energy_left_2>=0 and energy_right_2>=0) {
                min_E_temp = energy_left_2+energy_right_2;
//                cout<<"task left, middle, right: "<<task_left<<", "<<task_middle<<", "<<task_right<<", energy_2 "<<min_E_temp<<endl;
                best_broken_Edges_temp.assign(brokenEdges_2.begin(), brokenEdges_2.end());
                best_blocks_temp.assign(blocks_2.begin(), blocks_2.end());
                best_cores_temp.assign(cores_2.begin(), cores_2.end());
                best_copies_temp.assign(copies_2.begin(), copies_2.end());
                best_block_temp = id_block_used_2;
                best_core_temp = id_core_used_2;
                best_case_triplicated_temp = triplicated_2;
            }
        }
        
        if (min_E>=0 and min_E_temp>=0) {
            if (min_E_temp<min_E) {
                min_E = min_E_temp;
                break_node = task_middle;
                best_broken_Edges = best_broken_Edges_temp;
                best_blocks = best_blocks_temp;
                best_cores = best_cores_temp;
                best_copies = best_copies_temp;
                best_block = best_block_temp;
                best_core = best_core_temp;
                best_case_triplicated = best_case_triplicated_temp;
            }
        } else {
            if (min_E_temp>=0) {
                min_E = min_E_temp;
                break_node = task_middle;
                best_broken_Edges = best_broken_Edges_temp;
                best_blocks = best_blocks_temp;
                best_cores = best_cores_temp;
                best_copies = best_copies_temp;
                best_block = best_block_temp;
                best_core = best_core_temp;
                best_case_triplicated = best_case_triplicated_temp;
            }
        }
    }
    
//    if (task_right>26) {
//        cout<<"---, task left "<<task_left<<", task right "<<task_right<<", block "<<best_block<<", core "<<best_core<<endl;
//    }
    
    if (break_node!=-1) {
//        cout<<", break node "<<break_node;
        best_broken_Edges.push_back(break_node);
        best_blocks.push_back(best_block);
        best_cores.push_back(best_core);
        best_copies.push_back(best_case_triplicated);
        
        brokenEdges->assign(best_broken_Edges.begin(), best_broken_Edges.end());
        blocks->assign(best_blocks.begin(), best_blocks.end());
        cores->assign(best_cores.begin(), best_cores.end());
        copies->assign(best_copies.begin(), best_copies.end());
    } else {
        min_E = energy_toge;
        blocks->push_back(id_block_used_0);
        cores->push_back(id_core_used_0);
        copies->push_back(triplicated_0);
    }
    
//    cout<<", energy cost "<<min_E<<endl;
    
    return min_E;
}

double GetStructureEndNode(Graph* graph, unsigned int StartNode, unsigned int* endNode){
    double weight=graph->GetNode(StartNode)->GetWeight();
    stack<unsigned int> waiting_nodes;
    waiting_nodes.push(StartNode);
    Node* currendNode;
    unsigned int end_node_temp;
    double weight_temp;
    
    while (!waiting_nodes.empty()) {
        currendNode = graph->GetNode(waiting_nodes.top());
//        cout<<"current node is "<<currendNode->GetId()<<endl;
        waiting_nodes.pop();
        for (vector<unsigned int>::iterator iter = currendNode->GetSuccessors()->begin(); iter!=currendNode->GetSuccessors()->end(); ++iter) {
            if (graph->GetNode(*iter)->GetNbrSuccessor()>1) {//this node is a fork node
//                cout<<"node "<<(*iter)<<" is a fork node"<<endl;
                weight_temp = GetStructureEndNode(graph, (*iter), &end_node_temp);
                weight+=weight_temp;
//                cout<<"its join node is "<<end_node_temp<<endl;
//                cout<<"push join node "<<end_node_temp<<" into waiting stack"<<endl;
                waiting_nodes.push(end_node_temp);
            }
            
            if (graph->GetNode(*iter)->GetNbrSuccessor()<=1) {
                if (graph->GetNode(*iter)->GetPredecessors()->size()<=1) {
//                    cout<<"push node "<<(*iter)<<" into waiting stack"<<endl;
                    waiting_nodes.push((*iter));
//                    cout<<"plus weight of node "<<*iter<<endl;
                    weight+=graph->GetNode(*iter)->GetWeight();
                } else {//this node is a join node
                    *endNode = (*iter);
                }
            }
        }
    }
    
//    cout<<"plus weight of node "<<*endNode<<endl;
    weight+=graph->GetNode(*endNode)->GetWeight();
    
    return weight;
}

double GetStructureBeginNode(Graph* graph, unsigned int endNode, unsigned int* startNode){
    double weight=graph->GetNode(endNode)->GetWeight();
    stack<unsigned int> waiting_nodes;
    waiting_nodes.push(endNode);
    Node* currendNode;
    unsigned int start_node_temp;
    double weight_temp;
    
    while (!waiting_nodes.empty()) {
        currendNode = graph->GetNode(waiting_nodes.top());
//        cout<<"current node is "<<currendNode->GetId()<<endl;
        waiting_nodes.pop();
        for (vector<unsigned int>::iterator iter = currendNode->GetPredecessors()->begin(); iter!=currendNode->GetPredecessors()->end(); ++iter) {
//            cout<<*iter<<endl;
            if (graph->GetNode(*iter)->GetPredecessors()->size()>1) {//this node is a join node
//                cout<<"node "<<(*iter)<<" is a join node"<<endl;
                weight_temp = GetStructureBeginNode(graph, (*iter), &start_node_temp);
                weight+=weight_temp;
//                cout<<"its fork node is "<<start_node_temp<<endl;
//                cout<<"push fork node "<<start_node_temp<<" into waiting stack"<<endl;
                waiting_nodes.push(start_node_temp);
            }
            
            if (graph->GetNode(*iter)->GetPredecessors()->size()<=1) {
                if (graph->GetNode(*iter)->GetNbrSuccessor()<=1) {
//                    cout<<"push node "<<(*iter)<<" into waiting stack"<<endl;
                    waiting_nodes.push((*iter));
//                    cout<<"plus weight of node "<<*iter<<endl;
                    weight+=graph->GetNode(*iter)->GetWeight();
                } else {//this node is a fork node
                    *startNode = (*iter);
                }
            }
        }
    }
    
//    cout<<"plus weight of node "<<*startNode<<endl;
    weight+=graph->GetNode(*startNode)->GetWeight();
    
    return weight;
}

double GetEMaxS(const double weight){
    double energy = Static_Power*period+ConstantC*pow(SMAX, 2)*weight;
//    cout<<"weight "<<weight<<", running at maxS costs energy "<<energy<<endl;
    
    return energy;
}

double GetETrip(const double weight, double ouptput_edge_size){
    double speed;
    for (unsigned int i=0; i<Speed_Options; ++i) {
        speed = DiscSpeed[i];
        if (weight/speed<=period) {
            break;
        }
    }
    
    double energy = 3*(Static_Power*period+ConstantC*pow(speed, 2)*weight)+2*Alpha[0]*ouptput_edge_size;
//    cout<<"weight "<<weight<<", running at speed "<<speed<<", three copies costs energy "<<energy<<endl;
    
    return energy;
}

///Given a fork node and join node
///return if the output edge of fork node or the input edge of join node is too big
bool EdgeTooLarge(Graph* graph, unsigned int fork, unsigned int join){
    bool tooLarge = false;
    vector<unsigned int>* successors;
    
    successors = graph->GetNode(fork)->GetSuccessors();
    for (vector<unsigned int>::iterator iter=successors->begin(); iter!=successors->end(); ++iter) {
        if (graph->GetNode((*iter))->GetInputEdgeSize()->front()>period*BandWidth[0]) {
            tooLarge = true;
            break;
        }
    }
    
    for (vector<double>::iterator iter=graph->GetNode(join)->GetInputEdgeSize()->begin(); iter!=graph->GetNode(join)->GetInputEdgeSize()->end(); ++iter) {
        if ((*iter)>period*BandWidth[0]) {
            tooLarge = true;
            break;
        }
    }
    
    return tooLarge;
}

///Given a fork node, return if it has an edge that is too big
void EncountFork(Graph* graph, unsigned int fork, bool* taken_whole, double* part_weight, unsigned int* join){
    unsigned int join_node;
    double part_weight_temp = GetStructureEndNode(graph, fork, &join_node);
    *part_weight = part_weight_temp;
    *join = join_node;
    bool edge_too_large = EdgeTooLarge(graph, fork, join_node);
    
    if (edge_too_large == true) {//the fork-join has to be taken as a whole
        *taken_whole = true;
    } else {
        *taken_whole = false;
    }
}


/// break all edges except those that are too large
/// put parts into maxlist or triplicationList
void Traversal(Graph* graph, list<tuple<unsigned int, unsigned int, double, double>>* List_MaxS, list<tuple<unsigned int, unsigned int, double, double>>* List_Trip){
    stack<unsigned int> waiting_nodes;
    waiting_nodes.push(graph->GetSource()->GetId());
    unsigned int current_node;
    unsigned int part_start_node, part_end_node, join_node;
    double part_weight, part_weight_temp;
    vector<unsigned int>* successors;
    double output_edge_size, input_edge_size;
    bool taken_whole;
    double beta_used = BandWidth[0];
    vector<double>::iterator input_edge_size_iter;
    vector<unsigned int>::iterator predecessor_iter;
    
    //waiting_nodes stores part_start_node
    while (!waiting_nodes.empty()) {
        current_node = waiting_nodes.top();
        part_start_node = current_node;
        waiting_nodes.pop();
        part_weight = 0;
        
        while (true) {//try to find part_end_node, part_weight, output_edge_size
            part_end_node=current_node;
            part_weight += graph->GetNode(current_node)->GetWeight();
            successors = graph->GetNode(current_node)->GetSuccessors();
            if(successors->size()==0){
                break;
            }
            if (successors->size()==1) {//a chain
                current_node = successors->front();
                output_edge_size = graph->GetNode(current_node)->GetInputEdgeSize()->front();
//                if (output_edge_size<=BandWidth[0]*period) {
//                    break;
//                }
                if (output_edge_size<=beta_used*period) {
                    break;
                }
//                part_weight+=graph->GetNode(current_node)->GetWeight();
            } else if (successors->size()>1) {//a fork node
                EncountFork(graph, current_node, &taken_whole, &part_weight_temp, &join_node);
                if (taken_whole==true) {//the fork-join has to be taken as a whole
                    part_weight += part_weight_temp;
                    part_weight -= graph->GetNode(current_node)->GetWeight();
                    current_node = join_node;
                    part_weight -= graph->GetNode(current_node)->GetWeight();
                } else {
                    output_edge_size = 0;
                    successors = graph->GetNode(current_node)->GetSuccessors();
                    for (vector<unsigned int>::iterator iter=successors->begin(); iter!=successors->end(); ++iter) {
                        predecessor_iter = graph->GetNode((*iter))->GetPredecessors()->begin();
                        input_edge_size_iter = graph->GetNode((*iter))->GetInputEdgeSize()->begin();
                        while ((*predecessor_iter)!=current_node) {
                            ++predecessor_iter;
                            ++input_edge_size_iter;
                        }
                        output_edge_size+=(*input_edge_size_iter);
                    }
                        
//                  cout<<"add node "<<join_node<<" into waiting stack"<<endl;
                    waiting_nodes.push(join_node);//in case that join_node is inserted repeatedly, we inser it only here and no other place
                    break;
                }
            }
        }
        
        input_edge_size = 0;
        for (vector<double>::iterator it=graph->GetNode(part_start_node)->GetInputEdgeSize()->begin(); it!=graph->GetNode(part_start_node)->GetInputEdgeSize()->end(); ++it) {
            input_edge_size=input_edge_size+(*it);
        }
        
        if (GetEMaxS(part_weight)<=GetETrip(part_weight, output_edge_size)) {
//            cout<<"add "<<part_start_node<<"--"<<part_end_node<<" into MaxSpeed list"<<endl;
            List_MaxS->push_front(make_tuple(part_start_node,part_end_node,part_weight,input_edge_size));
        } else {
//            cout<<"add "<<part_start_node<<"--"<<part_end_node<<" into Triplication list"<<endl;
            List_Trip->push_front(make_tuple(part_start_node,part_end_node,part_weight,input_edge_size));
        }
        
        for (vector<unsigned int>::iterator iter=successors->begin(); iter!=successors->end(); ++iter) {
            if (graph->GetNode((*iter))->GetPredecessors()->size()<=1) {//avoid add join nodes
//                cout<<"add "<<(*iter)<<" in to waiting stack"<<endl;
                waiting_nodes.push((*iter));
            }
        }
    }
}


bool compare_nonIncreasing(const tuple<unsigned int, unsigned int, double, double>& first, const tuple<unsigned int, unsigned int, double, double>& second){
    return get<3>(first) > get<3>(second);
}

///node_id is the start node of a part
///return the part whose start node is node_id
list<tuple<unsigned int, unsigned int, double, double>>::iterator FindByStartId(list<tuple<unsigned int, unsigned int, double, double>>* parts, unsigned int node_id){
    list<tuple<unsigned int, unsigned int, double, double>>:: iterator iterator_found = parts->end();
    for (list<tuple<unsigned int, unsigned int, double, double>>::iterator it=parts->begin(); it!=parts->end(); it++){
        if (get<0>((*it))==node_id) {
            iterator_found = it;
            break;
        }
    }
    
    return iterator_found;
}

///node_id is the end node of a part
///return the part whose end node is node_id
list<tuple<unsigned int, unsigned int, double, double>>::iterator FindByEndId(list<tuple<unsigned int, unsigned int, double, double>>* parts, unsigned int node_id){
    list<tuple<unsigned int, unsigned int, double, double>>:: iterator iterator_found = parts->end();
    for (list<tuple<unsigned int, unsigned int, double, double>>::iterator it=parts->begin(); it!=parts->end(); it++){
        if (get<1>((*it))==node_id) {
            iterator_found = it;
            break;
        }
    }
    
    return iterator_found;
}

///start node is a frok node,
///end node is a join node,
///keep start node and end_node, and delete other parts between start node and end_node
bool CheckAllMaxS(Graph* graph, unsigned int start_node, unsigned int end_node, double * total_w, bool clean, list<tuple<unsigned int, unsigned int, double, double>> *parts){
    queue<unsigned int> waiting_nodes;
    waiting_nodes.push(end_node);
    Node* currentNode;
    bool all_in = true;
    list<tuple<unsigned int, unsigned int, double, double>>:: iterator iter_found;
    double total_weight=0;
    unsigned int temp_start_node;
    iter_found = FindByStartId(parts, end_node);
    if (iter_found != parts->end()) {
//        cout<<"count part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
        total_weight = get<2>(*iter_found);
    } else {
        all_in = false;
        return all_in;
//        cout<<"CheckAllMaxS error 11"<<endl;
    }
    
    while (!waiting_nodes.empty()) {
        currentNode = graph->GetNode(waiting_nodes.front());
        waiting_nodes.pop();
        
        //consider fork part only once
        if (currentNode->GetPredecessors()->size()>1) {//this is a join node
            GetStructureBeginNode(graph, currentNode->GetId(), &temp_start_node);
//            cout<<", a join node"<<endl;
            if (temp_start_node != start_node) {
                if (graph->GetNode(temp_start_node)->GetCopies()==1) {
                    iter_found = FindByEndId(parts, temp_start_node);
                    if (iter_found!=parts->end()) {
                        waiting_nodes.push(get<0>(*iter_found));//count fork node of the fork-join
//                        cout<<"count part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
                        total_weight += get<2>(*iter_found);
                        if (clean==true) {
//                            cout<<" delete "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
                            parts->erase(iter_found);
                        }
                    } else {
                        cout<<"CheckAllMaxS error 15: "<<temp_start_node<<" does not found."<<endl;
                    }
                } else {
                    all_in = false;
                    return all_in;
                }
            }
        }

        for (vector<unsigned int>::iterator iter = currentNode->GetPredecessors()->begin(); iter!=currentNode->GetPredecessors()->end(); ++iter){
            if (graph->GetNode(*iter)->GetCopies()==1) {//this part is in MaxS speed option.
                iter_found = FindByEndId(parts, *iter);//iter_found contains the begin node and end node of the part
                if (iter_found!=parts->end()) {// this part has been found
                    if (graph->GetNode(get<1>(*iter_found))->GetSuccessors()->size()<=1) {//do not count fork node
                        if (get<0>(*iter_found)!=start_node and get<1>(*iter_found)!=start_node) {
                            waiting_nodes.push(get<0>(*iter_found));
//                            cout<<"count part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
                            total_weight += get<2>(*iter_found);
                            if (clean == true) {
//                                cout<<" delete "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<" from MaxS list"<<endl;
                                parts->erase(iter_found);
                            }
                        }
                    }
                }
            } else {
                all_in = false;
                return all_in;
            }
        }
    }
    
    iter_found = FindByEndId(parts, start_node);
    if (iter_found != parts->end()) {
//        cout<<"count part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
        total_weight += get<2>(*iter_found);
    }
//    } else {
//        cout<<"error 12"<<endl;
//    }
    
    *total_w = total_weight;
    return all_in;
}

///start node is a frok node,
///end node is a join node,
///keep start node and end_node, and delete other parts between start node and end_node
double CheckSumWeight(Graph* graph, unsigned int start_node, unsigned int end_node, bool clean, list<tuple<unsigned int, unsigned int, double, double>>* MaxS, list<tuple<unsigned int, unsigned int, double, double>>* Triplication){
    queue<unsigned int> waiting_nodes;
    waiting_nodes.push(end_node);
    Node* currentNode;
    list<tuple<unsigned int, unsigned int, double, double>>:: iterator iter_found;
    double sum_weight = 0;
    unsigned int temp_start_node;
    
    if (graph->GetNode(end_node)->GetCopies()==1) {//this part is in MaxS speed option.
        iter_found = FindByStartId(MaxS, end_node);
        if (iter_found != MaxS->end()) {
            sum_weight = get<2>(*iter_found);
//            cout<<"part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<"( "<<get<2>(*iter_found)<<" )"<<endl;
        } else {
//            cout<<"error 13"<<endl;
        }
    } else {
        iter_found = FindByStartId(Triplication, end_node);
        if (iter_found != Triplication->end()) {
            sum_weight = get<2>(*iter_found);
//            cout<<"part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<"( "<<get<2>(*iter_found)<<" )"<<endl;
        } else {
//            cout<<"error 14"<<endl;
        }
    }

    while (!waiting_nodes.empty()) {
        currentNode = graph->GetNode(waiting_nodes.front());
        waiting_nodes.pop();
        
        if (currentNode->GetPredecessors()->size()>1) {//this is a join node
            GetStructureBeginNode(graph, currentNode->GetId(), &temp_start_node);
//            cout<<", a join node"<<endl;
            if (temp_start_node != start_node ) {
                if (graph->GetNode(temp_start_node)->GetCopies()==1) {
                    iter_found = FindByEndId(MaxS, temp_start_node);
                    if (iter_found!=MaxS->end()) {
                        waiting_nodes.push(get<0>(*iter_found));//count fork node of the fork-join
                        sum_weight += get<2>(*iter_found);
//                        cout<<"part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<"( "<<get<2>(*iter_found)<<" )"<<endl;
                        if (clean==true) {
                            MaxS->erase(iter_found);
//                            cout<<" delete "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
                        }
                    } else {
                        cout<<"CheckSumWeight error 15: "<<temp_start_node<<" does not found."<<endl;
                    }
                } else {
                    iter_found = FindByEndId(Triplication, temp_start_node);
                    if (iter_found!=Triplication->end()) {
                        waiting_nodes.push(get<0>(*iter_found));//count fork node of the fork-join
                        sum_weight += get<2>(*iter_found);
//                        cout<<"part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<"( "<<get<2>(*iter_found)<<" )"<<endl;
                        if (clean==true) {
                            Triplication->erase(iter_found);
//                            cout<<" delete "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
                        }
                    } else {
                        cout<<"CheckSumWeight error 16: "<<temp_start_node<<" does not found."<<endl;
                    }
                }
            }
        }
        
        for (vector<unsigned int>::iterator iter = currentNode->GetPredecessors()->begin(); iter!=currentNode->GetPredecessors()->end(); ++iter){
            if (*iter!=start_node) {
                if (graph->GetNode(*iter)->GetSuccessors()->size()<=1) {//avoid accounting fork node repeatedly
                    if (graph->GetNode(*iter)->GetCopies()==1) {//this part is in MaxS speed option.
                        iter_found = FindByEndId(MaxS, *iter);
                        if (iter_found!=MaxS->end()) {// this part has been found
                            waiting_nodes.push(get<0>(*iter_found));
                            sum_weight += get<2>(*iter_found);
//                            cout<<"part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<"( "<<get<2>(*iter_found)<<" )"<<endl;
                            if (clean==true) {
//                                cout<<" delete "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
                                MaxS->erase(iter_found);
                            }
                        } else {
                            cout<<"error 9: part "<<*iter<<" is not found."<<endl;
                        }
                    } else {// this part is in Triplication
                        iter_found = FindByEndId(Triplication, *iter);
                        if (iter_found!=Triplication->end()) {// this part has been found
                            waiting_nodes.push(get<0>(*iter_found));
                            sum_weight += get<2>(*iter_found);
//                            cout<<"part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<"( "<<get<2>(*iter_found)<<" )"<<endl;
                            if (clean==true) {
//                                cout<<" delete "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
                                Triplication->erase(iter_found);
                            }
                        } else {
                            cout<<"error 10"<<endl;
                        }
                    }
                }
            }
        }
    }
    
    if (graph->GetNode(start_node)->GetCopies()==1) {//this part is in MaxS speed option.
        iter_found = FindByEndId(MaxS, start_node);
        if (iter_found != MaxS->end()) {
            sum_weight += get<2>(*iter_found);
//            cout<<"part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<"( "<<get<2>(*iter_found)<<" )"<<endl;
        }
    } else {
        iter_found = FindByEndId(Triplication, start_node);
        if (iter_found != Triplication->end()) {
            sum_weight += get<2>(*iter_found);
//            cout<<"part "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<"( "<<get<2>(*iter_found)<<" )"<<endl;
        }
    }
    
    return sum_weight;
}

void GroupCell(Graph* graph, int block, int core, vector<unsigned int>* brokenEs){
    list<tuple<unsigned int, unsigned int, double, double>> MaxSpeedList;//start_node, end_node, weight, edge_size
    list<tuple<unsigned int, unsigned int, double, double>> TriplicationList;
    Traversal(graph, &MaxSpeedList, &TriplicationList);
    
    MaxSpeedList.sort(compare_nonIncreasing);//by an non-increasing order of input edge’s weight;
    
//    cout<<"Parts in MaxSpeed list"<<endl;
    for (list<tuple<unsigned int, unsigned int, double, double>>::iterator it=MaxSpeedList.begin(); it!=MaxSpeedList.end(); ++it) {
        graph->GetNode(get<0>((*it)))->SetCopies(1);
        graph->GetNode(get<1>((*it)))->SetCopies(1);
//        cout<<get<0>(*it)<<"-"<<get<1>(*it)<<" ( "<<get<2>(*it)<<","<<get<3>(*it)<<")"<<endl;
    }
//    cout<<endl;
    
//    cout<<"Parts in Triplication list"<<endl;
    for (list<tuple<unsigned int, unsigned int, double, double>>::iterator it=TriplicationList.begin(); it!=TriplicationList.end(); ++it) {
        graph->GetNode(get<0>((*it)))->SetCopies(3);
        graph->GetNode(get<1>((*it)))->SetCopies(3);
//        cout<<get<0>(*it)<<"-"<<get<1>(*it)<<" ( "<<get<2>(*it)<<","<<get<3>(*it)<<")"<<endl;
    }
//    cout<<endl;
    
    long size_before=0;
    Node* predecessor;
    list<tuple<unsigned int, unsigned int, double, double>>:: iterator iterator_found;
    list<tuple<unsigned int, unsigned int, double, double>>:: iterator iterator_found_end_node;
    double temp_size;
    unsigned int start_node, end_node;
    bool all_in;
    while (MaxSpeedList.size()!=size_before) {
        size_before = MaxSpeedList.size();
        for (list<tuple<unsigned int, unsigned int, double, double>>::iterator search=MaxSpeedList.begin(); search!=MaxSpeedList.end();) {
//            cout<<"try to merge "<<get<0>(*search)<<"-"<<get<1>(*search);
            if (graph->GetNode(get<0>((*search)))->GetPredecessors()->size()>1) {//this is a join node
                GetStructureBeginNode(graph, get<0>((*search)) , &start_node);
                all_in = CheckAllMaxS(graph, start_node, get<0>((*search)), &temp_size, false, &MaxSpeedList);
                if (all_in==true and temp_size<=period*SMAX) {
                    iterator_found = FindByEndId(&MaxSpeedList, start_node);
//                    cout<<", "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found))<<" to "<<get<0>((*iterator_found))<<"-"<<get<1>((*search))<<endl;
                    get<1>((*iterator_found))=get<1>((*search));
                    get<2>((*iterator_found))=temp_size;
                    all_in = CheckAllMaxS(graph, start_node, get<0>((*search)), &temp_size, true, &MaxSpeedList);//clean
//                    cout<<"delete "<<get<0>((*search))<<"-"<<get<1>((*search))<<endl;
                    search = MaxSpeedList.erase(search);
                } else {
//                    cout<<endl;
                    ++search;
                }
            } else {
                if (get<0>(*search)==graph->GetSource()->GetId()) {//the source node has no predecessor
                    ++search;
                } else {
                    predecessor = graph->GetNode(graph->GetNode(get<0>((*search)))->GetPredecessors()->front());
                    if (predecessor->GetSuccessors()->size()>1) {//its predecessor is a fork node
                        temp_size = GetStructureEndNode(graph, predecessor->GetId(), &end_node);
                        all_in = CheckAllMaxS(graph, predecessor->GetId(), end_node, &temp_size, false, &MaxSpeedList);
                        if (all_in==true and temp_size<=period*SMAX) {
                            iterator_found = FindByEndId(&MaxSpeedList, predecessor->GetId());
                            iterator_found_end_node = FindByStartId(&MaxSpeedList, end_node);
//                            cout<<", "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found))<<" to "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found_end_node))<<endl;
                            get<1>((*iterator_found))=get<1>((*iterator_found_end_node));
                            get<2>((*iterator_found))=temp_size;
                            all_in = CheckAllMaxS(graph, predecessor->GetId(), end_node, &temp_size, true, &MaxSpeedList);
//                            cout<<"delete "<<get<0>((*iterator_found_end_node))<<"-"<<get<1>((*iterator_found_end_node))<<endl;
                            MaxSpeedList.erase(iterator_found_end_node);
                            search = MaxSpeedList.begin();
                        } else {
                            ++search;
                        }
                    } else {// it's a chain
                        iterator_found = FindByEndId(&MaxSpeedList, predecessor->GetId());
                        if (iterator_found!=MaxSpeedList.end()) {
                            if ((get<2>((*iterator_found))+get<2>((*search)))<=period*SMAX) {
//                                cout<<", "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found))<<" to "<<get<0>((*iterator_found))<<"-"<<get<1>((*search))<<endl;
                                get<1>((*iterator_found))=get<1>((*search));
                                get<2>((*iterator_found))=get<2>((*iterator_found))+get<2>((*search));
//                                cout<<"delete "<<get<0>((*search))<<"-"<<get<1>((*search))<<endl;
                                search = MaxSpeedList.erase(search);
                            } else {
                                ++search;
                            }
                        } else {
                            ++search;
                        }
                    }
                }
            }
        }
    }
    
    unsigned int number_cores = block*core;
    unsigned int current_node;
    unsigned int new_copies, temp_node_one;
    if ((MaxSpeedList.size()+3*TriplicationList.size())>number_cores) {
        TriplicationList.sort(compare_nonIncreasing);//by an non-increasing order of input edge’s weight;
        for (list<tuple<unsigned int, unsigned int, double, double>>::iterator part=TriplicationList.begin(); part!=TriplicationList.end();) {
            current_node = get<0>(*part);
//            cout<<"try to merge "<<current_node<<"-"<<get<1>(*part);
            if (current_node != graph->GetSource()->GetId()) {//source node has no predecessor
                if (graph->GetNode(current_node)->GetPredecessors()->size()>1) {// this is a join node
                    GetStructureBeginNode(graph, current_node , &start_node);
                    temp_size = CheckSumWeight(graph, start_node, current_node, false, &MaxSpeedList, &TriplicationList);
                    if (temp_size<=period*SMAX) {
                        if (graph->GetNode(start_node)->GetCopies()==1) {//start_node is in MaxS
                            iterator_found = FindByEndId(&MaxSpeedList, start_node);
                            if (iterator_found!=MaxSpeedList.end()) {
//                                cout<<", "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found))<<" to "<<get<0>((*iterator_found))<<"-"<<get<1>(*part)<<endl;
                                get<1>(*iterator_found)=get<1>(*part);
                                get<2>(*iterator_found)=temp_size;
                                new_copies = 1;
                                temp_node_one = get<1>(*iterator_found);
                            } else {
                                cout<<"error 5"<<endl;
                            }
                        } else {//start_node is in Triplication
                            iterator_found = FindByEndId(&TriplicationList, start_node);
                            if (iterator_found!=TriplicationList.end()) {
//                                cout<<", "<<get<0>(*iterator_found)<<"-"<<get<1>(*iterator_found)<<" to "<<get<0>((*iterator_found))<<"-"<<get<1>(*part)<<endl;
                                get<1>(*iterator_found)=get<1>(*part);
                                get<2>(*iterator_found)=temp_size;
                                new_copies = 3;
                                temp_node_one = get<1>(*iterator_found);
                            } else {
                                cout<<"error 6"<<endl;
                            }
                        }
                        CheckSumWeight(graph, start_node, current_node, true, &MaxSpeedList, &TriplicationList);
//                        cout<<" delete "<<get<0>(*part)<<"-"<<get<1>(*part)<<endl;
                        part = TriplicationList.erase(part);
                        graph->GetNode(temp_node_one)->SetCopies(new_copies);
                        graph->GetNode(current_node)->SetCopies(new_copies);
                    } else {
//                        cout<<endl;
                        ++part;
                    }
                } else {
                    predecessor = graph->GetNode(graph->GetNode(current_node)->GetPredecessors()->front());
                    if (predecessor->GetSuccessors()->size()>1) {//its predecessor is a fork node
                        GetStructureEndNode(graph, predecessor->GetId(), &end_node);
                        temp_size = CheckSumWeight(graph, predecessor->GetId(), end_node, false, &MaxSpeedList, &TriplicationList);
                        if (temp_size<=period*SMAX) {
                            if (graph->GetNode(predecessor->GetId())->GetCopies()==1) {//in max list
                                iterator_found = FindByEndId(&MaxSpeedList, predecessor->GetId());
                            } else {
                                iterator_found = FindByEndId(&TriplicationList, predecessor->GetId());
                            }
                                
                            if (graph->GetNode(end_node)->GetCopies()==1) {//in max list
                                iterator_found_end_node = FindByStartId(&MaxSpeedList, end_node);
                            } else {
                                iterator_found_end_node = FindByStartId(&TriplicationList, end_node);
                                graph->GetNode(get<1>(*iterator_found_end_node))->SetCopies(1);
                            }
                            
//                            cout<<", "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found))<<" to "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found_end_node))<<endl;
                            
                            get<1>((*iterator_found))=get<1>(*iterator_found_end_node);
                            get<2>((*iterator_found))=temp_size;

                            CheckSumWeight(graph, predecessor->GetId(), end_node, true, &MaxSpeedList, &TriplicationList);
//                            cout<<" delete "<<get<0>(*iterator_found_end_node)<<"-"<<get<1>(*iterator_found_end_node)<<endl;
                            
                            if (graph->GetNode(end_node)->GetCopies()==1) {//in max list
                                MaxSpeedList.erase(iterator_found_end_node);
                                part = TriplicationList.begin();
                            } else {
                                part = TriplicationList.erase(iterator_found_end_node);
                            }
                            
                            if (graph->GetNode(predecessor->GetId())->GetCopies()==1) {//in max list
                                graph->GetNode(get<1>((*iterator_found)))->SetCopies(1);
                            } else {
                                graph->GetNode(get<1>(*iterator_found))->SetCopies(3);
                            }
                        } else {
//                            cout<<endl;
                            ++part;
                        }
                    } else {// it's a chain
                        if (graph->GetNode(predecessor->GetId())->GetCopies()==1) {
                            iterator_found = FindByEndId(&MaxSpeedList, predecessor->GetId());
                        } else {
                            iterator_found = FindByEndId(&TriplicationList, predecessor->GetId());
                        }
                            
                        if ((get<2>((*iterator_found))+get<2>(*part))<=period*SMAX) {
//                            cout<<", "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found))<<" to "<<get<0>((*iterator_found))<<"-"<<get<1>((*part))<<endl;
                            get<1>((*iterator_found))=get<1>((*part));
                            get<2>((*iterator_found))=get<2>((*iterator_found))+get<2>((*part));
                            
                            if (graph->GetNode(predecessor->GetId())->GetCopies()==1) {
                                graph->GetNode(get<1>(*part))->SetCopies(1);
                            } else {
                                graph->GetNode(get<1>(*part))->SetCopies(3);
                            }
                                
//                            cout<<" delete "<<get<0>(*part)<<"-"<<get<1>(*part)<<endl;
                            part = TriplicationList.erase(part);
                        } else {
//                            cout<<endl;
                            ++part;
                        }
                    }
            }
        } else {
            ++part;
        }
            
            if ((MaxSpeedList.size()+3*TriplicationList.size())<=number_cores) {
                break;
            }
            
            if (TriplicationList.empty()) {
                break;
            }
        }
    }
    
    double output_edge_size = 0;
    vector<unsigned int>* successors;
    
    vector<double>::iterator input_edge_size_iter;
    vector<unsigned int>::iterator predecessor_iter;

    for (list<tuple<unsigned int, unsigned int, double, double>>::iterator it=TriplicationList.begin(); it!=TriplicationList.end();) {
        output_edge_size = 0;
        successors = graph->GetNode(get<1>(*it))->GetSuccessors();
        for (vector<unsigned int>::iterator iter=successors->begin(); iter!=successors->end(); ++iter) {
            predecessor_iter = graph->GetNode((*iter))->GetPredecessors()->begin();
            input_edge_size_iter = graph->GetNode((*iter))->GetInputEdgeSize()->begin();
            while ((*predecessor_iter)!=get<1>(*it)) {
                ++predecessor_iter;
                ++input_edge_size_iter;
            }
            
            output_edge_size+=(*input_edge_size_iter);
        }
        
        if (GetEMaxS(get<2>(*it))<=GetETrip(get<2>(*it), output_edge_size)) {
            MaxSpeedList.push_back(*it);
            graph->GetNode(get<0>(*it))->SetCopies(1);
            graph->GetNode(get<1>(*it))->SetCopies(1);
            it = TriplicationList.erase(it);
        } else {
            brokenEs->push_back(get<0>(*it));
            ++it;
        }
    }
    
    for (list<tuple<unsigned int, unsigned int, double, double>>::iterator it=MaxSpeedList.begin(); it!=MaxSpeedList.end(); ++it) {
        brokenEs->push_back(get<0>(*it));
    }
    
    return;
}

///given a node, it returns the chain that starts from this node, the chain may end at a fork node
///if the given node is a fork, it returns the fork node itself
///if the given node is a join node, the chain starts from here
unsigned int AddPart(Graph* graph, list<tuple<unsigned int, unsigned int, double, double>>* targetList, Node* node){
    //node should be a start node of a chain
    Node* current_node = node;
    double weight = current_node->GetWeight();
    double input_edge_size = 0;
    unsigned int end_node = current_node->GetId();
    
    vector<double>* input_edges = current_node->GetInputEdgeSize();
    for (vector<double>::iterator iter=input_edges->begin(); iter!=input_edges->end(); ++iter) {
        input_edge_size += (*iter);
    }
    
    vector<unsigned int>* successors=current_node->GetSuccessors();
    while (successors->size()==1) {
        current_node = graph->GetNode(successors->front());
        if (current_node->GetPredecessors()->size()==1) {//current_node is not a join node
            weight+=current_node->GetWeight();
            successors = current_node->GetSuccessors();
            end_node = current_node->GetId();
        } else { // current_node is a fork node or join node
            break;
        }
    }
    
    targetList->push_back(make_tuple(node->GetId(), end_node, weight, input_edge_size));
    
    if (current_node->GetSuccessors()->size()>1) {//current_node is a fork node
        return current_node->GetId();
    }
    
    return node->GetId();
}

///Break output edges of fork nodes,
///input edges of join nodes,
///parts generated are stored in Parts
void BreakForkJoins(Graph* graph, list<tuple<unsigned int, unsigned int, double, double>>* Parts){
    queue<unsigned int> waiting_nodes;
    waiting_nodes.push(graph->GetSource()->GetId());
    unsigned int current_node, return_node, join_node;
    vector<unsigned int>* successors;
    
    while (!waiting_nodes.empty()) {
        current_node = waiting_nodes.front();
        waiting_nodes.pop();
        if (graph->GetNode(current_node)->GetSuccessors()->size()<=1) {
            //a chain node or a join
            return_node = AddPart(graph, Parts, graph->GetNode(current_node));
            if (return_node!=current_node) {// return_node is a fork node
                successors = graph->GetNode(return_node)->GetSuccessors();
                for (vector<unsigned int>::iterator it=successors->begin(); it!=successors->end(); ++it) {
                    waiting_nodes.push(*it);
                }
                
                if (graph->GetNode(return_node)->GetSuccessors()->size()>1) {//the return node is a fork node
                    GetStructureEndNode(graph, return_node, &join_node);
                    waiting_nodes.push(join_node);
                }
            }
        } else {// a fork node
            Parts->push_back(make_tuple(current_node,current_node,graph->GetNode(current_node)->GetWeight(),graph->GetNode(current_node)->GetInputEdgeSize()->front()));
            successors = graph->GetNode(current_node)->GetSuccessors();
            for (vector<unsigned int>::iterator iter = successors->begin(); iter!=successors->end(); ++iter) {
                waiting_nodes.push(*iter);
            }
            GetStructureEndNode(graph, current_node, &join_node);
            waiting_nodes.push(join_node);
        }
    }
}

///end node may be a fork node
///keep start node and end_node, and delete other parts between start node and end_node
double CheckSumWeight(Graph* graph, unsigned int start_node, unsigned int end_node, bool clean, list<tuple<unsigned int, unsigned int, double, double>>* Parts){
    queue<unsigned int> waiting_nodes;
    waiting_nodes.push(end_node);
    Node* currentNode;
    list<tuple<unsigned int, unsigned int, double, double>>:: iterator iter_found;
    double sum_weight = 0;
    unsigned int temp_start_node;
    
    iter_found = FindByStartId(Parts, end_node);
    if (iter_found != Parts->end()) {
        sum_weight = get<2>(*iter_found);
    } else {
        cout<<"CheckSumWeight error 15"<<endl;
    }
    
    while (!waiting_nodes.empty()) {
        currentNode = graph->GetNode(waiting_nodes.front());
        waiting_nodes.pop();
//        cout<<"current node "<<currentNode->GetId();
        
        if (currentNode->GetPredecessors()->size()>1) {//this is a join node
            GetStructureBeginNode(graph, currentNode->GetId(), &temp_start_node);
//            cout<<", a join node"<<endl;
            if (temp_start_node != start_node ) {
                iter_found = FindByEndId(Parts, temp_start_node);
                if (iter_found!=Parts->end()) {
                    waiting_nodes.push(get<0>(*iter_found));//count fork node of the fork-join
                    sum_weight += get<2>(*iter_found);
                    if (clean) {
                        Parts->erase(iter_found);
//                        cout<<" delete "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
                    }
//                    cout<<", add node "<<temp_start_node<<" into queue."<<endl;
                } else {
                    cout<<"CheckSumWeight error 14: "<<temp_start_node<<" does not found."<<endl;
                }
            }
        }
        
        for (vector<unsigned int>::iterator iter=currentNode->GetPredecessors()->begin(); iter!=currentNode->GetPredecessors()->end(); ++iter){
            if (*iter!=start_node) {
                if (graph->GetNode(*iter)->GetSuccessors()->size()<=1) {// do not count fork node here, otherwise it will be consider repeatedly
                    iter_found = FindByEndId(Parts, *iter);
                    if (iter_found!=Parts->end()) {// this part has been found
//                        cout<<"count node "<<*iter<<"push "<<get<0>(*iter_found)<<" into queue."<<endl;
                        waiting_nodes.push(get<0>(*iter_found));
                        sum_weight += get<2>(*iter_found);
                        if (clean==true) {
//                            cout<<" delete "<<get<0>(*iter_found)<<"-"<<get<1>(*iter_found)<<endl;
                            Parts->erase(iter_found);
                        }
                    } else {
                        cout<<"CheckSumWeight error 16: "<<*iter<<" does not found."<<endl;
                    }
                }
            }
        }
    }
    
    iter_found = FindByEndId(Parts, start_node);
    if (iter_found != Parts->end()) {
        sum_weight += get<2>(*iter_found);
    }
    
    return sum_weight;
}

void AddPartby2Ends(Graph* graph, unsigned int start_node, unsigned int end_node, list<tuple<unsigned int, unsigned int, double, double>>* targetList){
    unsigned int current_node = end_node;
    double part_weight = 0;
    vector<unsigned int>* predecessors;
    
    while (current_node!=start_node) {
        part_weight += graph->GetNode(current_node)->GetWeight();
        predecessors = graph->GetNode(current_node)->GetPredecessors();
        current_node = predecessors->front();
    }
    
    part_weight += graph->GetNode(start_node)->GetWeight();
    
    double input_edge_size = 0;
    vector<double>* input_edges_size = graph->GetNode(start_node)->GetInputEdgeSize();
    for (vector<double>::iterator it=input_edges_size->begin(); it!=input_edges_size->end(); ++it) {
        input_edge_size += (*it);
    }
    
    targetList->push_back(make_tuple(start_node,end_node,part_weight,input_edge_size));//start_node, end_node, weight, edge_weight
//    cout<<start_node<<"-"<<end_node<<"("<<part_weight<<"), ";
}

bool compare_nonDecreasingWeight(const tuple<unsigned int, unsigned int, double, double>& first, const tuple<unsigned int, unsigned int, double, double>& second){
    return get<2>(first) < get<2>(second);
}

bool CheckFeasibility(Graph* graph, unsigned int task_left, unsigned int task_right){
    bool feasible = true;
    bool too_heavy_on_one_core = false;
    bool all_edges_too_big = true;
    double total_weight = 0;
    double edge_size = 0;

    for (int i = task_left; i <= task_right; ++i){
        total_weight += graph->GetNode(i)->GetWeight();
    }

    if (total_weight>(SMAX*period)){
        too_heavy_on_one_core = true;
    }

    if (too_heavy_on_one_core==true){
        for (int i=task_left+1; i<=task_right; ++i) {
            edge_size = graph->GetNode(i)->GetInputEdgeSize()->front();
            if (edge_size<=(BandWidth[0]*period)) {
                all_edges_too_big = false;
                break;
            }
        }
    }
    
    if (too_heavy_on_one_core==true and all_edges_too_big==true) {
        feasible = false;
    }

    return feasible;
}

//old version of BreakFJ-DP of our paper
void BreakParaComb(Graph* graph, int block, int core, vector<unsigned int>* BrokenEdges){
    list<tuple<unsigned int, unsigned int, double, double>> Parts;//start_node, end_node, weight, input_edge_size
    BreakForkJoins(graph, &Parts);
    
//    cout<<"Partitions: "<<endl;
//    for (list<tuple<unsigned int, unsigned int, double, double>>::iterator it=Parts.begin(); it!=Parts.end(); ++it) {
//        cout<<get<0>(*it)<<"-"<<get<1>(*it)<<" ("<<get<2>(*it)<<", "<<get<3>(*it)<<")"<<endl;
//    }
//    cout<<endl;
    
    list<tuple<unsigned int, unsigned int, double, double>> Parts_Large;
    //parts whose size are larger than P_t*S_max;
    for (list<tuple<unsigned int, unsigned int, double, double>>::iterator it=Parts.begin(); it!=Parts.end(); ++it) {
        if (get<2>(*it)>period*SMAX) {
//            cout<<"part "<<get<0>(*it)<<"-"<<get<1>(*it)<<" is too large"<<endl;
            Parts_Large.push_back(*it);
        }
    }
    
    Parts_Large.sort(compare_nonDecreasingWeight);//by a non-decreasing weight
    tuple<unsigned int, unsigned int, double, double> Part;
    unsigned int num_block = 1, num_core = 1;
    forward_list<pair<int, int> > IdleCores;
    
    vector<unsigned int> brokenEs;
    vector<unsigned int> Blocks;
    vector<unsigned int> Cores;
    vector<unsigned int> Copies;
    
    unordered_map<string, tuple<double, unsigned int, unsigned int>> Register;
    
    double energy=-1;
    vector<unsigned int>::reverse_iterator iter_break_edge;
    vector<unsigned int>::reverse_iterator iter_break_edge_plus1;
    list<tuple<unsigned int, unsigned int, double, double>>:: iterator iterator_found;
    list<tuple<unsigned int, unsigned int, double, double>>:: iterator iterator_found_end_node;
    unsigned int end_node, start_node, part_start_node;
    bool feasible;
    while (!Parts_Large.empty()) {
        Part = Parts_Large.front();
        Parts_Large.pop_front();
        
        num_core = 1;
        energy = -1;
        
        feasible = CheckFeasibility(graph, get<0>(Part), get<1>(Part));
        if (feasible==true) {
            do {
                num_core++;
                Register.clear();
                if (num_core>(get<1>(Part)-get<0>(Part)+1)) {
                    break;
                }
                energy = Dynamic_Get_E(graph, get<0>(Part), get<1>(Part), num_block, num_core, &Register);
            } while (energy==-1);
//            cout<<"break part "<<get<0>(Part)<<"-"<<get<1>(Part)<<", into ";
        }
        
        if (energy!=-1) {
            brokenEs.clear();
            Blocks.clear();
            Cores.clear();
            Copies.clear();
            InterpertRegister(&Register, &brokenEs, &Blocks, &Cores, &Copies, get<0>(Part), get<1>(Part), num_block, num_core, false);
            
            //delte the original part
            start_node = get<0>(Part);
            end_node = get<1>(Part);
            iterator_found = FindByEndId(&Parts, end_node);
            if (iterator_found!=Parts.end()) {
                Parts.erase(iterator_found);
            } else {
                cout<<"BreakParaComb error 1"<<endl;
            }
            
            //create new parts
            ///the input edge of a node is broken, so the end_node is brokenEs -1
            brokenEs.pop_back();
            AddPartby2Ends(graph, start_node, brokenEs.back()-1, &Parts);
            
            iter_break_edge = brokenEs.rbegin();
            iter_break_edge_plus1 = iter_break_edge+1;
            for (; iter_break_edge_plus1!=brokenEs.rend();) {
                AddPartby2Ends(graph, *iter_break_edge, (*iter_break_edge_plus1)-1, &Parts);
                ++iter_break_edge;
                ++iter_break_edge_plus1;
            }
            
            AddPartby2Ends(graph, *iter_break_edge, end_node, &Parts);
            //        cout<<endl;
        }
    }
    
    long number_cores = block*core;
    double part_size;
    unsigned long size_before;
    Node* predecessor;
    while (Parts.size()>number_cores) {//if there is more parts than cores
        size_before = Parts.size();
        Parts.sort(compare_nonDecreasingWeight);// sort parts by a non-decreasing order of weight
        
        for (list<tuple<unsigned int, unsigned int, double, double>>::iterator part=Parts.begin(); part!=Parts.end();){
            part_start_node = get<0>(*part);
//            cout<<"try to merge "<<part_start_node<<"-"<<get<1>(*part);
            if (part_start_node != graph->GetSource()->GetId()) {//source node has no predecessor
                if (graph->GetNode(part_start_node)->GetPredecessors()->size()>1) {// this is a join node
                    GetStructureBeginNode(graph, part_start_node , &start_node);
                    part_size = CheckSumWeight(graph, start_node, part_start_node, false, &Parts);
                    if (part_size<=period*SMAX) {
                        iterator_found = FindByEndId(&Parts, start_node);
                        if (iterator_found!=Parts.end()) {
//                            cout<<", "<<get<0>(*iterator_found)<<"-"<<get<1>(*iterator_found)<<" to "<<get<0>((*iterator_found))<<"-"<<get<1>(*part)<<endl;
                            get<1>(*iterator_found)=get<1>(*part);
                            get<2>(*iterator_found)=part_size;
                        } else {
                            cout<<"error 16"<<endl;
                        }
                        CheckSumWeight(graph, start_node, part_start_node, true, &Parts);
//                        cout<<" delete "<<get<0>(*part)<<"-"<<get<1>(*part)<<endl;
                        part = Parts.erase(part);
                        break;
                    } else {
//                        cout<<endl;
                        ++part;
                    }
                } else {
                    predecessor = graph->GetNode(graph->GetNode(part_start_node)->GetPredecessors()->front());
                    if (predecessor->GetSuccessors()->size()>1) {//its predecessor is a fork node
                        GetStructureEndNode(graph, predecessor->GetId(), &end_node);
                        part_size = CheckSumWeight(graph, predecessor->GetId(), end_node, false, &Parts);
                        if (part_size<=period*SMAX) {
                            iterator_found = FindByEndId(&Parts, predecessor->GetId());
                            iterator_found_end_node = FindByStartId(&Parts, end_node);
//                            cout<<", "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found))<<" to "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found_end_node))<<endl;
                            
                            get<1>((*iterator_found))=get<1>(*iterator_found_end_node);
                            get<2>((*iterator_found))=part_size;

                            CheckSumWeight(graph, predecessor->GetId(), end_node, true, &Parts);
//                            cout<<" delete "<<get<0>(*iterator_found_end_node)<<"-"<<get<1>(*iterator_found_end_node)<<endl;
                            
                            Parts.erase(iterator_found_end_node);
                            part = Parts.begin();
                            break;
                        } else {
                            ++part;
                        }
                    } else {// it's a chain
                        iterator_found = FindByEndId(&Parts, predecessor->GetId());
                            
                        if ((get<2>((*iterator_found))+get<2>(*part))<=period*SMAX) {
                            get<1>((*iterator_found))=get<1>((*part));
                            get<2>((*iterator_found))=get<2>((*iterator_found))+get<2>((*part));
//                            cout<<", "<<get<0>((*iterator_found))<<"-"<<get<1>((*iterator_found))<<" to "<<get<0>((*iterator_found))<<"-"<<get<1>((*part))<<endl;
//                            cout<<" delete "<<get<0>(*part)<<"-"<<get<1>(*part)<<endl;
                            part = Parts.erase(part);
                            break;
                        } else {
                            ++part;
                        }
                    }
            }
        } else {
            ++part;
        }
    }
        
        if (size_before == Parts.size()) {
            break;
        }
    }
    
    vector<unsigned int>* successors;
    double output_edge_size;
    vector<double>::iterator input_edge_size_iter;
    vector<unsigned int>::iterator predecessor_iter;
    for (list<tuple<unsigned int, unsigned int, double, double>>::iterator it=Parts.begin(); it!=Parts.end(); ++it) {
        part_size = get<2>((*it));
        successors = graph->GetNode(get<1>(*it))->GetSuccessors();
        output_edge_size = 0;
        for (vector<unsigned int>::iterator iter=successors->begin(); iter!=successors->end(); ++iter) {
            predecessor_iter = graph->GetNode((*iter))->GetPredecessors()->begin();
            input_edge_size_iter = graph->GetNode((*iter))->GetInputEdgeSize()->begin();
            while ((*predecessor_iter)!=get<1>(*it)) {
                ++predecessor_iter;
                ++input_edge_size_iter;
            }
            
            output_edge_size+=(*input_edge_size_iter);
        }
        
        if (GetEMaxS(part_size)<=GetETrip(part_size, output_edge_size)){
            graph->GetNode(get<0>((*it)))->SetCopies(1);
            graph->GetNode(get<1>((*it)))->SetCopies(1);
        } else {
            graph->GetNode(get<0>((*it)))->SetCopies(3);
            graph->GetNode(get<1>((*it)))->SetCopies(3);
        }
        
        BrokenEdges->push_back(get<0>(*it));
    }
    
    return;
}

bool BreakFJ_DP(Graph* graph, const int nbr_block, const int nbr_core, vector<unsigned int>* BrokenEdges){
    list<tuple<unsigned int, unsigned int, double, double>> Parts;//start_node, end_node, weight, input_edge_size
    BreakForkJoins(graph, &Parts);
    
    tuple<unsigned int, unsigned int, double, double> Part;
    forward_list<pair<int, int> > IdleCores;
    
    vector<unsigned int> local_brokenEs;
    vector<unsigned int> Blocks;
    vector<unsigned int> Cores;
    vector<unsigned int> Copies;
    vector<unsigned int>::iterator it_copies;
    double energy = 0;
    unordered_map<string, tuple<double, unsigned int, unsigned int>> Register;
    BrokenEdges->clear();
    
    while (!Parts.empty()) {
        Part = Parts.front();
        Parts.pop_front();
        
        energy = Dynamic_Get_E(graph, get<0>(Part), get<1>(Part), nbr_block, nbr_core, &Register);

        if (energy==-1) {
            return false;
        } else {
            local_brokenEs.clear();
            Blocks.clear();
            Cores.clear();
            Copies.clear();
            InterpertRegister(&Register, &local_brokenEs, &Blocks, &Cores, &Copies, get<0>(Part), get<1>(Part), nbr_block, nbr_core, false);
            
            it_copies = Copies.begin();
            for (vector<unsigned int>::iterator it=local_brokenEs.begin(); it!=local_brokenEs.end(); ++it, ++it_copies) {
                if (*it==0) {//0 means the first node is broken, then we add the start node of this part
                    BrokenEdges->push_back(get<0>(Part));
                    graph->GetNode(get<0>(Part))->SetCopies(*it_copies);
                } else {
                    BrokenEdges->push_back(*it);
                    graph->GetNode(*it)->SetCopies(*it_copies);
                }
            }
        }
    }
    
    return true;
}

bool SaveProc_chain(Graph* graph, unsigned int task_left, unsigned int task_right, vector<unsigned int>* brokenEdges){
//    cout<<"SaveProc_Chain: "<<task_left<<" to "<<task_right<<endl;
    int nbr_core = 1;
    double energy;
    int block = 1;
    bool success = true;
    forward_list<pair<int, int> > IdleC;
    vector<unsigned int> brokenE;
    vector<unsigned int> blocks;
    vector<unsigned int> cores;
    vector<unsigned int> copies;
    unordered_map<string, tuple<double, unsigned int, unsigned int>> Map;
    
    bool feasible = CheckFeasibility(graph, task_left, task_right);
    
    if (feasible==true) {
        do {
            nbr_core++;
            if (nbr_core> (task_right-task_left+1)) {
                success = false;
                break;
            }
            
            Map.clear();
            energy = Dynamic_Get_E(graph, task_left, task_right, block, nbr_core, &Map);
        } while (energy==-1);
    } else {
        success = false;
    }
    
    if (success == true) {
        InterpertRegister(&Map, &brokenE, &blocks, &cores, &copies, task_left, task_right, block, nbr_core, false);
        brokenE.pop_back();
        
        brokenEdges->assign(brokenE.begin(), brokenE.end());
    }
    
    return success;
}

unsigned int GetChainEndNode(Graph* graph, unsigned int start_node, double* total_weight){
    unsigned int current_node = start_node;
    double weight = graph->GetNode(current_node)->GetWeight();
    vector<unsigned int>* successors = graph->GetNode(current_node)->GetSuccessors();
    while (successors->size()==1) {
        if (graph->GetNode(successors->front())->GetPredecessors()->size()<=1 and graph->GetNode(successors->front())->GetSuccessors()->size()<=1) {
            current_node = successors->front();
            successors = graph->GetNode(current_node)->GetSuccessors();
            weight+=graph->GetNode(current_node)->GetWeight();
        } else {
            break;
        }
    }
    
    *total_weight = weight;
    
    return current_node;
}

//include task_left, but not task_right
void GetNodesBetween(Graph* graph, unsigned int task_left, unsigned int task_right, unsigned int* end_node){
    stack<unsigned int> waiting_nodes;
    waiting_nodes.push(task_left);
    Node* currendNode;
    unsigned int end_node_temp;
    
    while (!waiting_nodes.empty()) {
        currendNode = graph->GetNode(waiting_nodes.top());
//        cout<<"current node is "<<currendNode->GetId()<<endl;
        waiting_nodes.pop();
        for (vector<unsigned int>::iterator iter = currendNode->GetSuccessors()->begin(); iter!=currendNode->GetSuccessors()->end(); ++iter) {
            if ((*iter)!=task_right) {
                if (graph->GetNode(*iter)->GetNbrSuccessor()>1) {//this node is a fork node
                    //cout<<"node "<<(*iter)<<" is a fork node"<<endl;
                    GetStructureEndNode(graph, (*iter), &end_node_temp);
                    //cout<<"its join node is "<<end_node_temp<<endl;
                    //cout<<"push join node "<<end_node_temp<<" into waiting stack"<<endl;
                    waiting_nodes.push(end_node_temp);
                }
                
                if (graph->GetNode(*iter)->GetNbrSuccessor()<=1) {
                     waiting_nodes.push((*iter));
                }
            } else {
                *end_node = currendNode->GetId();
                break;
            }
        }
    }
    return;
}

void SaveProc_forkJoin(Graph* graph, unsigned int task_left, unsigned int task_right, vector<unsigned int>* brokenEs){
    //break output edges of task_left, a fork node
    //break input edges of task_right, a join node
//    cout<<"SaveProc-frokjoin: "<<task_left<<" to "<<task_right<<endl;
    vector<unsigned int>* successors = graph->GetNode(task_left)->GetSuccessors();
    for (vector<unsigned int>::iterator it=successors->begin(); it!=successors->end(); ++it) {
        brokenEs->push_back(*it);
    }
    brokenEs->push_back(task_right);

    vector<unsigned int> brokenE_temp;
    unsigned int end_node;
    for (vector<unsigned int>::iterator iter=successors->begin(); iter!=successors->end(); ++iter) {
        GetNodesBetween(graph, *iter, task_right, &end_node);
        brokenE_temp.clear();
//        cout<<"SaveProc: "<<*iter<<"-"<<end_node<<endl;
        SaveProc(graph, *iter, end_node, &brokenE_temp);
        brokenEs->insert(brokenEs->end(), brokenE_temp.begin(), brokenE_temp.end());
    }
    
    return;
}

///the structure to be break is supposed to be a series combination of chains and fork-joins, which may contain fork-joins
///if the
void BreakStructures(Graph* graph, unsigned int task_left, unsigned int task_right, list<tuple<unsigned int, unsigned int, double, double>>* Parts, vector<bool>* Label){
    queue<unsigned int> waiting_nodes;
    waiting_nodes.push(task_left);
    unsigned int current_node, return_node, join_node, start_node;
    vector<unsigned int>* successors;
    double edge_size, part_weight, weight;
    vector<unsigned int> empty_successor;
    empty_successor.clear();
    
    start_node = task_left;
    part_weight = 0;
    bool a_new_structue = true;
    while (!waiting_nodes.empty()) {
        current_node = waiting_nodes.front();
        waiting_nodes.pop();
        
//        cout<<"current node "<<current_node<<endl;
        if (graph->GetNode(current_node)->GetPredecessors()->size()<=1 and graph->GetNode(current_node)->GetSuccessors()->size()<=1) {
            //a chain node
            if (a_new_structue==false) {
                start_node = current_node;
                part_weight = 0;
                a_new_structue=true;
            }
            
            return_node = GetChainEndNode(graph, current_node, &weight);
            part_weight+=weight;
            if (return_node == task_right) {
                successors = &empty_successor;
            } else {
                successors = graph->GetNode(return_node)->GetSuccessors();
            }
            
            if (successors->size()==0) {
                edge_size = graph->GetNode(start_node)->GetInputEdgeSize()->front();
                Parts->push_back(make_tuple(start_node,return_node, part_weight, edge_size));
                Label->push_back(false);
                a_new_structue = false;
            }
        } else if (graph->GetNode(current_node)->GetSuccessors()->size()>1) {// a fork node
            weight = GetStructureEndNode(graph, current_node, &join_node);
            
            if (weight>period*SMAX) {
                if (a_new_structue==true and graph->GetNode(start_node)->GetSuccessors()->size()==1) {
                    edge_size = graph->GetNode(start_node)->GetInputEdgeSize()->front();
                    Parts->push_back(make_tuple(start_node,graph->GetNode(current_node)->GetPredecessors()->front(), part_weight, edge_size));
                    Label->push_back(false);
                }
                
                edge_size = graph->GetNode(current_node)->GetInputEdgeSize()->front();
                Parts->push_back(make_tuple(current_node,join_node, weight, edge_size));
                Label->push_back(true);
                a_new_structue = false;
            } else {
                part_weight+=weight;
                
                if (a_new_structue==false) {
                    start_node = current_node;
                    part_weight = weight;
                    a_new_structue=true;
                }
            }
            
            if (join_node == task_right) {
                successors = &empty_successor;
            } else {
                successors = graph->GetNode(join_node)->GetSuccessors();
            }
        }
        
        for (vector<unsigned int>::iterator iter = successors->begin(); iter!=successors->end(); ++iter) {
            waiting_nodes.push(*iter);
        }
    }
    
    return;
}

bool checkIfChain(Graph* graph, unsigned int task_left, unsigned int task_right){
    bool is_a_chain = true;
    
    Node* current_node = graph->GetNode(task_left);
    while (current_node->GetSuccessors()->size()==1) {
        if (current_node->GetId()==task_right) {
            break;
        }
        current_node=graph->GetNode(current_node->GetSuccessors()->front());
    }
    
    if (current_node->GetSuccessors()->size()>1) {
        is_a_chain = false;
    }
    
    return is_a_chain;
}

void AddNodetoGraph(Graph* new_graph, Node* original_node, double new_weight){
    new_graph->AddNode();
    unsigned int new_graph_node_id = new_graph->GetNodes()->size();
    new_graph->GetNode(new_graph_node_id)->SetWeight(new_weight);
    if (new_graph_node_id>1) {
        new_graph->GetNode(new_graph_node_id)->AddPredecessor(new_graph_node_id-1);
        new_graph->GetNode(new_graph_node_id-1)->AddSuccessor(new_graph_node_id);
    }
    new_graph->GetNode(new_graph_node_id)->AddInputEdgeSize(original_node->GetInputEdgeSize()->front());
}

//this function will change the stucture of the graph
void HandleMixStruc(Graph* graph, unsigned int task_left, unsigned int task_right, vector<unsigned int>* BrokenEdges){
    unsigned int end_node;
    double part_weight;
    unordered_map<unsigned int, unsigned int> Map;//new_graph_node_id to old_graph_node_id
    Graph* new_graph = new Graph();
    Node* current_node = graph->GetNode(task_left);
    double new_node_weight;
    
    while (current_node->GetId()!=task_right) {
//        cout<<"original node: "<<current_node->GetId()<<endl;
        if (current_node->GetSuccessors()->size()==1) {//a chain node
            new_node_weight = current_node->GetWeight();
            AddNodetoGraph(new_graph, current_node, new_node_weight);
            Map.insert(make_pair(new_graph->GetNodes()->size(), current_node->GetId()));
            current_node = graph->GetNode(current_node->GetSuccessors()->front());
        } else {//a fork node
            part_weight = GetStructureEndNode(graph, current_node->GetId(), &end_node);
            new_node_weight = part_weight;
            AddNodetoGraph(new_graph, current_node, new_node_weight);
            Map.insert(make_pair(new_graph->GetNodes()->size(), current_node->GetId()));
            
            if (end_node==task_right) {
                break;
            } else {
                if (graph->GetNode(end_node)->GetSuccessors()->size()==1) {//avoid the case that a join node is also a fork node
                    current_node = graph->GetNode(graph->GetNode(end_node)->GetSuccessors()->front());
                }
            }
        }
    }
    
    if (current_node->GetId()==task_right and current_node->GetPredecessors()->size()==1) {//a join node has already been considered
        new_node_weight = current_node->GetWeight();
        AddNodetoGraph(new_graph, current_node, new_node_weight);
        Map.insert(make_pair(new_graph->GetNodes()->size(),current_node->GetId()));
//        cout<<"original node: "<<current_node->GetId()<<endl;
    }
    new_graph->SetSize(new_graph->GetNodes()->size());
    
    
    vector<unsigned int> broken_Edges_temp;
    SaveProc_chain(new_graph, 1, new_graph->GetNodes()->size(), &broken_Edges_temp);
    
    for (vector<unsigned int>::iterator it=broken_Edges_temp.begin(); it!=broken_Edges_temp.end(); ++it) {
        BrokenEdges->push_back(Map[*it]);
    }
    
    delete new_graph;
}

void SaveProc_mix(Graph* graph, unsigned int task_left, unsigned int task_right, vector<unsigned int>* brokenEs){
//    cout<<"SaveProc_Mix: "<<task_left<<" to "<<task_right<<endl;
    list<tuple<unsigned int, unsigned int, double, double>> Parts;//start_node, end_node, weight, input_edge_weight
    vector<bool> Labels;
    
    BreakStructures(graph, task_left, task_right, &Parts, &Labels);
    
//    cout<<"Parts: "<<endl;
    for (list<tuple<unsigned int, unsigned int, double, double>>::iterator it=Parts.begin(); it!=Parts.end(); ++it) {
        graph->GetNode(get<0>((*it)))->SetCopies(1);
        graph->GetNode(get<1>((*it)))->SetCopies(1);
//        cout<<get<0>(*it)<<"-"<<get<1>(*it)<<" ( "<<get<2>(*it)<<", "<<get<3>(*it)<<")"<<endl;
    }
    
    vector<unsigned int> brokenEs_temp;
    list<tuple<unsigned int, unsigned int, double, double>>::iterator it_parts=Parts.begin();
    vector<bool>::iterator it_labels = Labels.begin();
    for (; it_parts!=Parts.end(); ) {
        brokenEs_temp.clear();
        if ((*it_labels)==true) {// a fork-join structure which is larger than period*Smax
            SaveProc_forkJoin(graph, get<0>((*it_parts)), get<1>((*it_parts)), &brokenEs_temp);
            brokenEs->insert(brokenEs->end(), brokenEs_temp.begin(), brokenEs_temp.end());
            it_parts = Parts.erase(it_parts);
            it_labels = Labels.erase(it_labels);
        } else {
            ++it_parts; ++it_labels;
        }
    }
    
    it_parts = Parts.begin();
    it_labels = Labels.begin();
    for (; it_parts!=Parts.end();it_parts++, it_labels++) {
        brokenEs_temp.clear();
        if (get<2>(*it_parts)>period*SMAX) {
            if (checkIfChain(graph, get<0>(*it_parts), get<1>(*it_parts))) {
                SaveProc_chain(graph, get<0>(*it_parts), get<1>(*it_parts), &brokenEs_temp);
            } else {
                HandleMixStruc(graph, get<0>(*it_parts), get<1>(*it_parts), &brokenEs_temp);
            }
            brokenEs->insert(brokenEs->end(), brokenEs_temp.begin(), brokenEs_temp.end());
        }
    }
    
    return;
}

double WhichStructure(Graph* graph, unsigned int task_left, unsigned int task_right, char* structure){
    double weight=0;
    queue<unsigned int> waiting_nodes;
    Node* currendNode;
    unsigned int end_node_temp;
    double weight_temp=0;
    vector<unsigned int>* successors;
    
    if (graph->GetNode(task_left)->GetSuccessors()->size()<=1) {
        //it maybe a chain or a mix, assume it is a chain
        bool isChain = checkIfChain(graph, task_left, task_right);
        if (isChain == true) {
            *structure = 'c';
        } else {
            *structure = 'm';
        }
        waiting_nodes.push(task_left);
    } else {
        //it is a parallel combination or a mix, assume it is a parallel combination
        weight_temp = GetStructureEndNode(graph, task_left, &end_node_temp);
        weight+=weight_temp;
        if (end_node_temp==task_right) {
            *structure = 'p';
        } else {
            *structure = 'm';
            successors = graph->GetNode(end_node_temp)->GetSuccessors();
            for (vector<unsigned int>::iterator iter=successors->begin(); iter!=successors->end(); ++iter) {
                if ((*iter)==task_right) {
                    weight+=graph->GetNode(*iter)->GetWeight();
                    return weight;
                } else {
                    waiting_nodes.push(*iter);
                }
            }
        }
    }
    
    while (!waiting_nodes.empty()) {
        currendNode = graph->GetNode(waiting_nodes.front());
        waiting_nodes.pop();
        
        if(currendNode->GetId()==task_right){
            weight = weight + currendNode->GetWeight();
            break;
        }
        
        if (currendNode->GetSuccessors()->size()>1) {//this is a frok node
            weight_temp = GetStructureEndNode(graph, currendNode->GetId(), &end_node_temp);
            successors = graph->GetNode(end_node_temp)->GetSuccessors();
        } else {
            //a chain structure
            weight_temp = currendNode->GetWeight();
            successors = currendNode->GetSuccessors();
        }
        weight+=weight_temp;
        
        for (vector<unsigned int>::iterator iter=successors->begin(); iter!=successors->end(); ++iter) {
            if (graph->GetNode((*iter))->GetPredecessors()->size()<=1) {//avoid add join nodes
                waiting_nodes.push((*iter));
            }
        }
    }
    
    return weight;
}

void SaveProc(Graph* graph, unsigned int task_left, unsigned int task_right, vector<unsigned int>* brokenEs){
    char structure;
    
    double weight = WhichStructure(graph, task_left, task_right, &structure);
    
    if (weight>period*SMAX) {
        switch (structure) {
            case 'c':
                SaveProc_chain(graph, task_left, task_right, brokenEs);
                break;
                
            case 'p':
                SaveProc_forkJoin(graph, task_left, task_right, brokenEs);
                break;
                
            case 'm':
                SaveProc_mix(graph, task_left, task_right, brokenEs);
                break;
                
            default:
                break;
        }
    }
    
    return;
}

void SaveProc_Entry(Graph* graph, unsigned int task_left, unsigned int task_right, vector<unsigned int>* brokenEs){
    SaveProc(graph, task_left, task_right, brokenEs);
    for (vector<unsigned int>::iterator it=brokenEs->begin(); it!=brokenEs->end(); ++it) {
        graph->GetNode(*it)->SetCopies(1);
    }
    graph->GetSource()->SetCopies(1);
}

///given original graph node, it returns predecessors' original graph node that correspond to Q grap node
void FindQPredecessors(unsigned int node, Graph* graph, unordered_map<unsigned int, unsigned int>* ID2QID, vector<unsigned int>* output){
    unsigned int current_id;
    vector<unsigned int>* predecessors;
    vector<unsigned int> waitingNodes;
    waitingNodes.push_back(node);
    unsigned int start_node;
    vector<unsigned int>::iterator iter;
    
    while (!waitingNodes.empty()) {
        current_id = waitingNodes.back();
        waitingNodes.pop_back();
        
        do {
            predecessors = graph->GetNode(current_id)->GetPredecessors();
            
            if (predecessors->size()>1) {//current node is a join node
                if (!graph->GetNode(current_id)->IsInputEdgesBroken()) {//not broken
                    GetStructureBeginNode(graph, current_id, &start_node);
                    current_id = start_node;
                } else {//broken
                    current_id = predecessors->front();
                    iter = predecessors->begin();
                    iter++;
                    for (; iter!=predecessors->end(); ++iter) {
                        if (graph->GetNode(*iter)->IsInputEdgesBroken()) {
                            output->push_back(*iter);
                        } else {
                            waitingNodes.push_back(*iter);
                        }
                    }
                }
            } else if(predecessors->size()==1) {
                current_id = predecessors->front();
            }
            
            if (graph->GetNode(current_id)->IsInputEdgesBroken()) {
                output->push_back(current_id);
            }
        } while (!graph->GetNode(current_id)->IsInputEdgesBroken());
    }
}

void BuildQGraph(Graph* graph, vector<unsigned int>* BrokenEs, Graph* Qgraph){
    bool include_one = false;
    for (vector<unsigned int>::iterator iter_edge=BrokenEs->begin(); iter_edge!=BrokenEs->end(); ++iter_edge) {
        if ((*iter_edge)==1) {
            include_one = true;
        }
    }

    unsigned long QGraph_size;
    if (include_one==false) {
//        cout<<"error in buildQgraph: root edge should always be broken"<<endl;
        QGraph_size = BrokenEs->size()+1;
    } else {
        QGraph_size = BrokenEs->size();
    }

    unordered_map<unsigned int, unsigned int> Qid_Id;
    unordered_map<unsigned int, unsigned int> Id_Qid;

    for (vector<unsigned int>::iterator it=BrokenEs->begin(); it!=BrokenEs->end(); ++it) {
        graph->GetNode(*it)->BreakInputEdges();
    }
    graph->GetSource()->BreakInputEdges();

    vector<unsigned int> Qid;
    vector<unsigned int> Qpredecessors;
    vector<unsigned int> QInputEdgesWeight;
    vector<double> QNodeWeight;

    //root
    unsigned int current_node;
    vector<unsigned int>* successors;
    vector<unsigned int> Nodes;
    double Qnode_weight = 0;
    queue<unsigned int> NextNodes;

    unsigned int current_Qnode=1;
    double Qinputedge=graph->GetSource()->GetInputEdgeSize()->front();
    Qpredecessors.push_back(0);
    QInputEdgesWeight.push_back(Qinputedge);

    NextNodes.push(graph->GetSource()->GetId());

    Qid_Id.insert(make_pair(1, graph->GetSource()->GetId()));
    Id_Qid.insert(make_pair(graph->GetSource()->GetId(),1));

    ///NextNodes stores node whose input edges are broken, Nodes stores nodes who start from a node with broken input edges
    while (!NextNodes.empty()) {
        Nodes.push_back(NextNodes.front());
//        cout<<"original_id to Q_id: "<<NextNodes.front()<<"-"<<current_Qnode<<endl;
        Qid_Id.insert(make_pair(current_Qnode, NextNodes.front()));
        Id_Qid.insert(make_pair(NextNodes.front(),current_Qnode));
        NextNodes.pop();
        Qnode_weight = 0;
        Qid.push_back(current_Qnode);
        while (!Nodes.empty()) {//calculate Qnode weight
            current_node = Nodes.back();
//            cout<<"   current node "<<current_node<<endl;
            Nodes.pop_back();
            Qnode_weight+=graph->GetNode(current_node)->GetWeight();
            successors = graph->GetNode(current_node)->GetSuccessors();
            for (vector<unsigned int>::iterator it=successors->begin(); it!=successors->end(); ++it) {
                if (graph->GetNode((*it))->GetPredecessors()->size()>1) {//a join node
                    if (current_node==graph->GetNode((*it))->GetPredecessors()->front()) {//make sure that join only count once
                        if (!graph->GetNode((*it))->IsInputEdgesBroken()) {
                            Nodes.push_back((*it));
                        } else {
                            NextNodes.push((*it));
                        }
                    }
                } else {
                    if (!graph->GetNode((*it))->IsInputEdgesBroken()) {
                        Nodes.push_back((*it));
                    } else {
                        NextNodes.push((*it));
                    }
                }
            }
        }
        QNodeWeight.push_back(Qnode_weight);

        current_Qnode++;
    }

    vector<unsigned int> predecessors;
    vector<unsigned int>::iterator Qid_iter=Qid.begin();
    ++Qid_iter;
    vector<double>* inputEdgesSize;
    vector<double>::iterator iter_Qweight = QNodeWeight.begin();
    for (; Qid_iter!=Qid.end(); ++Qid_iter, ++iter_Qweight) {
        current_Qnode = *Qid_iter;
//        cout<<"Qnode "<<current_Qnode<<", its node "<<Qid_Id[current_Qnode]<<endl;

        inputEdgesSize = graph->GetNode(Qid_Id[current_Qnode])->GetInputEdgeSize();
        for (vector<double>::iterator iter_weight=inputEdgesSize->begin(); iter_weight!=inputEdgesSize->end(); ++iter_weight) {
            if (iter_weight!=inputEdgesSize->begin()) {
                Qid_iter = Qid.insert(Qid_iter, current_Qnode);
                iter_Qweight = QNodeWeight.insert(iter_Qweight, *iter_Qweight);
                ++Qid_iter;
                ++iter_Qweight;
            }
            QInputEdgesWeight.push_back(*iter_weight);
        }

        predecessors.clear();
        FindQPredecessors(Qid_Id[current_Qnode], graph, &Id_Qid, &predecessors);
        for (vector<unsigned int>::iterator iter_node=predecessors.begin(); iter_node!=predecessors.end(); ++iter_node) {
            Qpredecessors.push_back(Id_Qid[*iter_node]);
        }
    }

//    cout<<"id predecessor_id"<<endl;
//    vector<unsigned int>::iterator it_temp = Qid.begin();
//    vector<unsigned int>::iterator it_temp_pred = Qpredecessors.begin();
//    for (; it_temp!=Qid.end(); ++it_temp, ++it_temp_pred) {
//        cout<<*it_temp<<" "<<*it_temp_pred<<endl;
//    }

    Qgraph->AllocateNodes(QGraph_size);
    for (unsigned int i=0; i<Qid.size(); ++i) {
//        cout<<Qid[i]<<" "<<Qpredecessors[i]<<" "<<QInputEdgesWeight[i]<<" "<<QNodeWeight[i]<<endl;
        Qgraph->GetNode(Qid[i])->SetId(Qid[i]);
        Qgraph->GetNode(Qid[i])->SetWeight(QNodeWeight[i]);
        Qgraph->GetNode(Qid[i])->SetCopies(graph->GetNode(Qid_Id[Qid[i]])->GetCopies());
        Qgraph->GetNode(Qid[i])->AddPredecessor(Qpredecessors[i]);
        if (Qpredecessors[i]!=0) {
            Qgraph->GetNode(Qpredecessors[i])->AddSuccessor(Qid[i]);
        }
        Qgraph->GetNode(Qid[i])->AddInputEdgeSize(QInputEdgesWeight[i]);

        Qgraph->GetNode(Qid[i])->SetBlock(graph->GetNode(Qid_Id[Qid[i]])->GetBlock());
        Qgraph->GetNode(Qid[i])->SetCore(graph->GetNode(Qid_Id[Qid[i]])->GetCore());
    }

    Qgraph->SetSource(1);

    return;
}

void DFS(Graph* graph, Node* task_i, vector<unsigned int>* DFSVector){
    vector<unsigned int>* successors = task_i->GetSuccessors();
    for (vector<unsigned int>::iterator iter=successors->begin(); iter!=successors->end(); ++iter) {
        if (!graph->GetNode(*iter)->IsVisited()) {
            DFS(graph, graph->GetNode(*iter), DFSVector);
        }
    }
    DFSVector->push_back(task_i->GetId());
    task_i->SetVisited();
    return;
}

bool Map_same_place_parallel_part(Graph* Qgraph, unsigned int part){
    bool success = false;
    vector<Node*> parallel_parts;
    parallel_parts.clear();
    unsigned int predecessor, successor;
    vector<unsigned int>* successors;
    
//    part must be a chain, otherwise it does not have a parallel part
    if (Qgraph->GetNode(part)->GetPredecessors()->size()==1 and Qgraph->GetNode(part)->GetSuccessors()->size()==1) {
        predecessor = Qgraph->GetNode(part)->GetPredecessors()->front();
        if(predecessor == 0){//part is the source part, it does not have a parallel part
            return false;
        }
        successor = Qgraph->GetNode(part)->GetSuccessors()->front();
        if (Qgraph->GetNode(predecessor)->GetSuccessors()->size()>1 and Qgraph->GetNode(successor)->GetPredecessors()->size()>1) {//it is a branch between fork and join
            successors = Qgraph->GetNode(predecessor)->GetSuccessors();
            for (vector<unsigned int>::iterator it=successors->begin(); it!=successors->end(); ++it) {
                if ((*it)!=part) {
                    parallel_parts.push_back(Qgraph->GetNode(*it));
                }
            }
        } else {
            return false;
        }
    } else {
        return false;
    }
    
    Node* mini_size_part;
    if (parallel_parts.empty()) {
        return false;
    } else {
        double smallest_weight=numeric_limits<double>::max();
        for (vector<Node*>::iterator it=parallel_parts.begin(); it!=parallel_parts.end(); ++it) {
            if ((*it)->GetId()!=part) {
                if ((*it)->GetWeight()<smallest_weight) {
                    smallest_weight = (*it)->GetWeight();
                    mini_size_part = *it;
                }
            }
        }
    }
    
    double sum_weight = Qgraph->GetNode(part)->GetWeight() + mini_size_part->GetWeight();
    double sum_output_weight = 0;
    vector<double>::iterator iter = Qgraph->GetNode(successor)->GetInputEdgeSize()->begin();
    for (vector<unsigned int>::iterator it=Qgraph->GetNode(successor)->GetPredecessors()->begin(); it!=Qgraph->GetNode(successor)->GetPredecessors()->end(); ++it, ++iter) {
        if (*it==part or *it==mini_size_part->GetId()) {
            sum_output_weight += (*iter);
        }
    }
    
    unsigned int block_id, core_id, copies;
    if (sum_weight<=period*SMAX) {
        success = true;
        block_id = mini_size_part->GetBlock();
        core_id = mini_size_part->GetCore();
        copies = mini_size_part->GetCopies();
        
        if (GetEMaxS(sum_weight) <= GetETrip(sum_weight, sum_output_weight)) {
            if (copies==3) {
                mini_size_part->SetCopies(1);
            }
            Qgraph->GetNode(part)->SetBlock(block_id);
            Qgraph->GetNode(part)->SetCore(core_id);
            Qgraph->GetNode(part)->SetCopies(1);
        } else {
            Qgraph->GetNode(part)->SetBlock(block_id);
            Qgraph->GetNode(part)->SetCore(core_id);
            Qgraph->GetNode(part)->SetCopies(3);
        }
    } else {
        return false;
    }
    
    return success;
}

//simple version of Map_Topology
//bool Map_Topology(vector<unsigned int>* BrokenEs, vector<unsigned int>* ProcessorsLeft, int block, Graph* Qgraph){
//    vector<unsigned int> DFSVector;
//    for (vector<Node*>::iterator iter=Qgraph->GetNodes()->begin(); iter!=Qgraph->GetNodes()->end(); ++iter) {
//        (*iter)->SetUnVisited();
//    }
//    DFS(Qgraph, Qgraph->GetSource(), &DFSVector);
//
//    unsigned int current_part, first_predecessor;
//    int temp_block;
//    unsigned int copies;
//    bool success;
//    double inputEdge_size = 0;
//    while (!DFSVector.empty()) {
//        current_part = DFSVector.back();
//        DFSVector.pop_back();
////        cout<<"current_part "<<current_part;
//        inputEdge_size = 0;
//
//        if (Qgraph->GetNode(current_part)->GetBlock()>0) {//already mapped by MapRanked
////            cout<<", already mapped by MapRanked."<<endl;
//        } else {
//            temp_block = block;
//            first_predecessor = 0;
//            if (Qgraph->GetNode(current_part)->GetPredecessors()->size()>0) {
//                first_predecessor = Qgraph->GetNode(current_part)->GetPredecessors()->front();
//            }
//
//            if (first_predecessor!=0) {
//                if (Qgraph->GetNode(first_predecessor)->GetBlock()!=0) {//already mapped
//                    temp_block = Qgraph->GetNode(first_predecessor)->GetBlock();
//                    inputEdge_size = Qgraph->GetNode(current_part)->GetInputEdgeSize()->front();
//                }
//            }
//
//            success = false;
//            copies = Qgraph->GetNode(current_part)->GetCopies();
////            if (ProcessorsLeft->at(temp_block-1)<copies) {
////                if (inputEdge_size > (BandWidth[1]*period)) {
////                    success = Map_same_place_parallel_part(Qgraph, current_part);
////                } else {
////                    temp_block = block;
////                }
////            }
//
//            if (ProcessorsLeft->at(temp_block-1)<copies) {
//                temp_block = block;
//            }
//
//            if (success == false) {
//                while (ProcessorsLeft->at(temp_block-1)<copies) {
//                    temp_block--;
//                    if (temp_block<=0) {
//                        //try to allocate it to the same core as its parallel part
//                        success = Map_same_place_parallel_part(Qgraph, current_part);
//                        if (success==false) {
//                            copies = 1;
//                            temp_block = block;
//                            while (ProcessorsLeft->at(temp_block-1)<copies) {
//                                temp_block--;
//                                if (temp_block<=0) {
//                                    return false;
//                                }
//                            }
//                            Qgraph->GetNode(current_part)->SetCopies(1);
//                        } else {
//                            break;
//                        }
//                    }
//                }
//            }
//
//            if (success == false) {
//                if (copies==3 and ProcessorsLeft->at(temp_block-1)>=3) {
////                    cout<<", map onto block "<<temp_block<<", core "<<ProcessorsLeft->at(temp_block-1)<<", "<<ProcessorsLeft->at(temp_block-1)-1<<", "<<ProcessorsLeft->at(temp_block-1)-2<<endl;
//                    Qgraph->GetNode(current_part)->SetBlock(temp_block);
//                    Qgraph->GetNode(current_part)->SetCore(ProcessorsLeft->at(temp_block-1));
//                    ProcessorsLeft->at(temp_block-1)-=3;
//                } else if (copies==1 and ProcessorsLeft->at(temp_block-1)>=1) {
////                    cout<<", map onto block "<<temp_block<<", core "<<ProcessorsLeft->at(temp_block-1)<<endl;
//                    Qgraph->GetNode(current_part)->SetBlock(temp_block);
//                    Qgraph->GetNode(current_part)->SetCore(ProcessorsLeft->at(temp_block-1));
//                    ProcessorsLeft->at(temp_block-1)-=1;
//                } else {
//                    return false;
//                }
//            }
//        }
//    }
//
//    return true;
//}

//MapRank with simple version of MapTopology
//bool Map_Ranked(vector<unsigned int>* BrokenEs, vector<unsigned int>* ProcessorsLeft, int block, Graph* Qgraph){
//    vector<unsigned int> DFSVector;
//    DFS(Qgraph, Qgraph->GetSource(), &DFSVector);
//
//    int temp_block = block;
//    vector<tuple<unsigned int, unsigned int, double>> Vec;//store parts that are connected by a large edge
//    double max_comm_size = BandWidth[1]*period;
//    vector<unsigned int>::iterator it_pred;
//    Node* current_node;
//    for (vector<unsigned int>::reverse_iterator it=DFSVector.rbegin(); it!=DFSVector.rend(); ++it) {
//        current_node = Qgraph->GetNode(*it);
//        it_pred = current_node->GetPredecessors()->begin();
//        for (vector<double>::iterator iter=current_node->GetInputEdgeSize()->begin(); iter!=current_node->GetInputEdgeSize()->end(); ++iter,++it_pred) {
//            if ((*iter)>max_comm_size) {
//                Vec.push_back(make_tuple(*it_pred,*it,*iter));//left_part, right_part, edge_size
//            }
//        }
//    }
//
//    //put neighbour parts of Vec into the same vector
//    unordered_map<unsigned int, unsigned int> NodeId2VecId;
//    unsigned int part_left, part_right, vec_id;
//    unordered_map<unsigned int, unsigned int>::const_iterator it_part_left;
//    unordered_map<unsigned int, unsigned int>::const_iterator it_part_right;
//    unsigned int new_vec_id=0;
//    for (vector<tuple<unsigned int, unsigned int, double>>::iterator it=Vec.begin(); it!=Vec.end(); ++it) {
//        part_left = get<0>((*it));
//        part_right = get<1>((*it));
//        //        cout<<"part_left -- part_right "<<part_left<<" - "<<part_right<<endl;
//        it_part_left = NodeId2VecId.find(part_left);
//        it_part_right = NodeId2VecId.find(part_right);
//
//        if (it_part_right!=NodeId2VecId.end() and it_part_left==NodeId2VecId.end()) {//right_part is already in a vector
//            vec_id = (*it_part_right).second;
//            NodeId2VecId.insert(make_pair(part_left, vec_id));
//        }
//
//        if (it_part_left!=NodeId2VecId.end() and it_part_right==NodeId2VecId.end()) {//left_part is already in a vector
//            vec_id = (*it_part_left).second;
//            NodeId2VecId.insert(make_pair(part_right, vec_id));
//        }
//
//        if (it_part_left==NodeId2VecId.end() and it_part_right==NodeId2VecId.end()) {//neither is in a vector
//            NodeId2VecId.insert(make_pair(part_left, new_vec_id));
//            NodeId2VecId.insert(make_pair(part_right, new_vec_id));
//            ++new_vec_id;
//        }
//    }
//
//    vector<list<unsigned int>> L(new_vec_id);//L is a vector of vector, each vector stores parts that should be mapped onto
//    //the same block
//    for (unordered_map<unsigned int, unsigned int>::iterator it=NodeId2VecId.begin(); it!=NodeId2VecId.end(); ++it) {
//        vec_id = (*it).second;
//        L[vec_id].push_back((*it).first);
//    }
//
//    vector<unsigned int> current_vector;
//    unsigned int num_procs_needed = 0;
//    int min_w_node_id;
//    while (!L.empty()) {
//        current_vector.assign(L.back().begin(), L.back().end());
//        L.pop_back();
//
//        num_procs_needed = 0;
//        for (vector<unsigned int>::iterator it=current_vector.begin(); it!=current_vector.end(); ++it) {
//            num_procs_needed+=Qgraph->GetNode(*it)->GetCopies();
//        }
//
//        min_w_node_id = -1;
//        while (num_procs_needed>Config_Core) {//parts in a vector demand too many cores
//            for (vector<unsigned int>::iterator it=current_vector.begin(); it!=current_vector.end(); ++it) {
//                if (Qgraph->GetNode(*it)->GetCopies()==3) {
//                    min_w_node_id = *it;
//                    break;
//                }
//            }
//
//            for (vector<unsigned int>::iterator it=current_vector.begin(); it!=current_vector.end(); ++it) {
//                if (Qgraph->GetNode(*it)->GetCopies()==3) {
//                    if (Qgraph->GetNode(*it)->GetWeight()<Qgraph->GetNode(min_w_node_id)->GetWeight()) {
//                        min_w_node_id = *it;
//                    }
//                }
//            }
//
//            if(min_w_node_id != -1){
//                if (Qgraph->GetNode(min_w_node_id)->GetCopies()==3) {
//                    Qgraph->GetNode(min_w_node_id)->SetCopies(1);//reduce requested cores by 2
//                    num_procs_needed = num_procs_needed - 2;
//                } else {
//                    return false;
//                }
//            } else {
//                return false;
//            }
//        }
//
//        if (ProcessorsLeft->at(temp_block-1)<num_procs_needed) {
//            temp_block--;
//        }
//
//        if (temp_block<1) {
//            return false;
//        }
//
//        for (vector<unsigned int>::iterator it=current_vector.begin(); it!=current_vector.end(); ++it) {
//            //            cout<<"current_part "<<*it;
//            Qgraph->GetNode(*it)->SetBlock(temp_block);
//            Qgraph->GetNode(*it)->SetCore(ProcessorsLeft->at(temp_block-1));
//            if (Qgraph->GetNode(*it)->GetCopies()==3) {
////                cout<<", map onto block "<<temp_block<<", core "<<ProcessorsLeft->at(temp_block-1)<<", "<<ProcessorsLeft->at(temp_block-1)-1<<", "<<ProcessorsLeft->at(temp_block-1)-2<<endl;
//                ProcessorsLeft->at(temp_block-1) = ProcessorsLeft->at(temp_block-1)-3;
//            } else {
////                cout<<", map onto block "<<temp_block<<", core "<<ProcessorsLeft->at(temp_block-1)<<endl;
//                ProcessorsLeft->at(temp_block-1) = ProcessorsLeft->at(temp_block-1)-1;
//            }
//        }
//    }
//
//    bool success = Map_Topology(BrokenEs, ProcessorsLeft, block, Qgraph);
//
//    return success;
//}

void UpdateQgraph(Graph* Qgraph){
    Node* current_part;
    unsigned int first_predecessor;

    vector<unsigned int> DFSVector;
    for (vector<Node*>::iterator iter=Qgraph->GetNodes()->begin(); iter!=Qgraph->GetNodes()->end(); ++iter) {
        (*iter)->SetUnVisited();
    }
    DFS(Qgraph, Qgraph->GetSource(), &DFSVector);

    while (!DFSVector.empty()) {
        current_part = Qgraph->GetNode(DFSVector.back());
        DFSVector.pop_back();

        if (!current_part->GetPredecessors()->empty()) {
            first_predecessor = current_part->GetPredecessors()->front();
            if (first_predecessor!=0) {//first_predecessor exists
                if (Qgraph->GetNode(first_predecessor)->GetCore()==current_part->GetCore()) {
                    if (Qgraph->GetNode(first_predecessor)->GetBlock()==current_part->GetBlock()) {
                        //these two parts are mapped onto the same core, we are gonna merge them into one part
                        //current_part will be deleted
                        Qgraph->MergeNode2Predeccsor(current_part);
//                        cout<<"merge part "<<current_part->GetId()<<" to part "<<first_predecessor<<endl;
                    }
                }
            }
        }
    }
}

//complex version of MapTopology
bool Map_Topology(Graph* Qgraph, const int block, vector<list<unsigned int>>* Sets){
    vector<unsigned int> DFSVector;
    for (vector<Node*>::iterator iter=Qgraph->GetNodes()->begin(); iter!=Qgraph->GetNodes()->end(); ++iter) {
        (*iter)->SetUnVisited();
    }
    DFS(Qgraph, Qgraph->GetSource(), &DFSVector);

    Node* current_part;
    unsigned int first_predecessor;
    double inputEdge_size=0;
    int current_block;
    unsigned int copies = 0;
    int index_h;
    Node* tempt_part;
    bool all_input_edge_small=false;
    double part_size, total_size, output_edge_size;
    double part_weight;
    vector<unsigned int>* successors;
    vector<double>::iterator input_edge_size_iter;
    vector<unsigned int>::iterator predecessor_iter;

//    vector<list<unsigned int>> Sets(block);//Sets is a vector of list, each list contains parts of a block
    vector<Node*> Head_parts_combined;
    bool map_success = true;
    unsigned int part_move;

    while (!DFSVector.empty()) {
        current_part = Qgraph->GetNode(DFSVector.back());
        DFSVector.pop_back();
        current_block = block;
        copies = current_part->GetCopies();

//        cout<<"current_part "<<current_part->GetId();

        if (current_part->GetBlock()>0) {//already mapped by MapRanked
//            cout<<", already mapped by MapRanked."<<endl;
        } else {
            if(!current_part->GetPredecessors()->empty()){
                first_predecessor = current_part->GetPredecessors()->front();
                if (first_predecessor!=0) {
                    if (Qgraph->GetNode(first_predecessor)->GetBlock()!=0) {//already mapped
                        current_block = Qgraph->GetNode(first_predecessor)->GetBlock();
                        inputEdge_size = current_part->GetInputEdgeSize()->front();
                    }
                }
            }

            if (Sets->at(current_block-1).size()>(Config_Core-copies) and inputEdge_size>(BandWidth[1]*period)) {
                index_h = 0;
                list<unsigned int>::reverse_iterator rit=Sets->at(current_block-1).rbegin();

                for (; rit!=Sets->at(current_block-1).rend(); ++rit) {
                    tempt_part = Qgraph->GetNode(*rit);
                    all_input_edge_small = true;
                    for (vector<double>::iterator size_it=tempt_part->GetInputEdgeSize()->begin(); size_it!=tempt_part->GetInputEdgeSize()->end(); ++size_it) {
                        if ((*size_it)>BandWidth[1]*period) {
                            all_input_edge_small = false;
                            break;
                        }
                    }

                    if (all_input_edge_small==true) {
                        if (tempt_part->GetId()==Qgraph->GetSource()->GetId()) {
                            all_input_edge_small = false;
                        }
                        break;
                    }
                }

                if (all_input_edge_small==true and current_block>1) {
//                    cout<<"move part between "<<*rit<<" and part "<<Sets[current_block-1].back()<<" to block "<<current_block-1<<endl;
                    part_move = *rit;
                    for (list<unsigned int>::reverse_iterator riter=Sets->at(current_block-1).rbegin(); riter!=rit; ++riter) {
                        Sets->at(current_block-2).push_back(*riter);
                        Qgraph->GetNode(*riter)->SetBlock(current_block-1);
                    }
                    Sets->at(current_block-2).push_back(part_move);
                    Qgraph->GetNode(part_move)->SetBlock(current_block-1);

                    while (Sets->at(current_block-1).back()!=part_move) {
                        Sets->at(current_block-1).pop_back();
                    }

                    while (Sets->at(current_block-1).back()==part_move) {
                        Sets->at(current_block-1).pop_back();
                    }

                    current_block = current_block - 1;
                } else {
                    while (Sets->at(current_block-1).size()>(Config_Core-copies)) {
                        part_size = numeric_limits<double>::max();
                        for (list<unsigned int>::iterator it=Sets->at(current_block-1).begin(); it!=Sets->at(current_block-1).end(); ++it) {
                            if (Qgraph->GetNode(*it)->GetNbrPredecessor()==1) {
                                if (Qgraph->GetNode((Qgraph->GetNode(*it)->GetPredecessors()->front()))->GetSuccessors()->size()==1) {
                                    //its predecessor is not a fork
                                    if (Qgraph->GetNode(*it)->GetWeight()<part_size) {
                                        tempt_part = Qgraph->GetNode(*it);
                                        part_size = Qgraph->GetNode(*it)->GetWeight();
                                    }
                                }
                            }
                        }

                        if (part_size!=numeric_limits<double>::max()) {//map tempt_part to the place where its predecessor is mapped
//                            cout<<"merge part "<<tempt_part->GetId()<<" and its predecessor "<<tempt_part->GetPredecessors()->front()<<endl;

                            for (list<unsigned int>::iterator it=Sets->at(current_block-1).begin(); it!=Sets->at(current_block-1).end();) {
                                if (*it==tempt_part->GetId()) {
                                    it = Sets->at(current_block-1).erase(it);
                                } else {++it;}
                            }

                            total_size = part_size + Qgraph->GetNode(tempt_part->GetPredecessors()->front())->GetWeight();

                            output_edge_size = 0;
                            successors = tempt_part->GetSuccessors();
                            part_weight = current_part->GetWeight();
                            for (vector<unsigned int>::iterator iter=successors->begin(); iter!=successors->end(); ++iter) {
                                predecessor_iter = Qgraph->GetNode((*iter))->GetPredecessors()->begin();
                                input_edge_size_iter = Qgraph->GetNode((*iter))->GetInputEdgeSize()->begin();
                                while ((*predecessor_iter)!=tempt_part->GetId()) {
                                    ++predecessor_iter;
                                    ++input_edge_size_iter;
                                }
                                output_edge_size+=(*input_edge_size_iter);
                            }

                            if (GetEMaxS(total_size)<=GetETrip(total_size, output_edge_size)) {
                                if (Qgraph->GetNode(tempt_part->GetPredecessors()->front())->GetCopies()==3) {
                                    /*remove two copies from Sets[current_block]*/
                                    for (list<unsigned int>::iterator it=Sets->at(current_block-1).begin(); it!=Sets->at(current_block-1).end();) {
                                        if (*it==tempt_part->GetId()) {
                                            it = Sets->at(current_block-1).erase(it);
                                            it = Sets->at(current_block-1).erase(it);
                                            break;
                                        } else {++it;}
                                    }
                                }

                                current_part->SetCopies(1);
                                tempt_part->SetCopies(1);
                            } else {
                                current_part->SetCopies(3);
                                tempt_part->SetCopies(3);
                            }

                            /*we label temp_part, current_part will be mapped later when temp_part will no longer change*/
                            Head_parts_combined.push_back(Qgraph->GetNode(tempt_part->GetPredecessors()->front()));
                        } else {
                            map_success = false;
                            break;
                        }
                    }
                }
            }

            if (map_success==true) {
                while (Sets->at(current_block-1).size()>(Config_Core-copies) and current_block>1) {
                    current_block--;
                }
                
                if (copies==3 and (Config_Core - Sets->at(current_block-1).size())>=3) {
//                    cout<<", map onto block "<<current_block<<", and core "<<Sets[current_block-1].size()+1<<" "<<Sets[current_block-1].size()+2<<" "<<Sets[current_block-1].size()+3<<endl;
                    current_part->SetBlock(current_block);
                    Sets->at(current_block-1).push_back(current_part->GetId());
                    Sets->at(current_block-1).push_back(current_part->GetId());
                    Sets->at(current_block-1).push_back(current_part->GetId());
                } else if(copies==1 and ((Config_Core - Sets->at(current_block-1).size())>0)){
//                    cout<<", map onto block "<<current_block<<", and core "<<Sets[current_block-1].size()+1<<endl;
                    current_part->SetBlock(current_block);
                    Sets->at(current_block-1).push_back(current_part->GetId());
                } else {
                    map_success = false;
                    break;
                }
            } else {
                break;
            }
        }
    }

    current_block = block;
    for (; current_block>0; --current_block) {
        for (list<unsigned int>::iterator it=Sets->at(current_block-1).begin(); it!=Sets->at(current_block-1).end(); ++it) {
//            cout<<"part "<<*it<<", block "<<current_block<<", core "<<(distance(Sets[current_block-1].begin(), it)+1)<<endl;
            Qgraph->GetNode(*it)->SetBlock(current_block);
            Qgraph->GetNode(*it)->SetCore(distance(Sets->at(current_block-1).begin(), it)+1);
        }
    }

    bool stop = false;
    while (stop == false) {
        stop = true;
        for (vector<Node*>::iterator it=Head_parts_combined.begin(); it!=Head_parts_combined.end(); ++it) {
            if((*it)->GetCore()>0) {
//                cout<<"map "<<(*it)->GetSuccessors()->front()<<" as "<<(*it)->GetId()<<endl;
                current_part = Qgraph->GetNode((*it)->GetSuccessors()->front());
                current_part->SetBlock((*it)->GetBlock());
                current_part->SetCore((*it)->GetCore());
            } else {
                stop = false;
            }
        }
    }

    UpdateQgraph(Qgraph);

    return map_success;
}

//MapRank with complex version of MapTopology
bool Map_Ranked(vector<unsigned int>* BrokenEs, vector<unsigned int>* ProcessorsLeft, int block, Graph* Qgraph){
    vector<unsigned int> DFSVector;
    DFS(Qgraph, Qgraph->GetSource(), &DFSVector);

    vector<list<unsigned int>> Sets(block);//Sets is a vector of list, each list contains parts of a block

    int temp_block = block;
    vector<tuple<unsigned int, unsigned int, double>> Vec;//store parts that are connected by a large edge
    double max_comm_size = BandWidth[1]*period;
    vector<unsigned int>::iterator it_pred;
    Node* current_node;
    for (vector<unsigned int>::reverse_iterator it=DFSVector.rbegin(); it!=DFSVector.rend(); ++it) {
        current_node = Qgraph->GetNode(*it);
        it_pred = current_node->GetPredecessors()->begin();
        for (vector<double>::iterator iter=current_node->GetInputEdgeSize()->begin(); iter!=current_node->GetInputEdgeSize()->end(); ++iter,++it_pred) {
            if ((*iter)>max_comm_size) {
                Vec.push_back(make_tuple(*it_pred,*it,*iter));//left_part, right_part, edge_size
            }
        }
    }

    //put neighbour parts of Vec into the same vector
    unordered_map<unsigned int, unsigned int> NodeId2VecId;
    unsigned int part_left, part_right, vec_id;
    unordered_map<unsigned int, unsigned int>::const_iterator it_part_left;
    unordered_map<unsigned int, unsigned int>::const_iterator it_part_right;
    unsigned int new_vec_id=0;
    for (vector<tuple<unsigned int, unsigned int, double>>::iterator it=Vec.begin(); it!=Vec.end(); ++it) {
        part_left = get<0>((*it));
        part_right = get<1>((*it));
//        cout<<"part_left -- part_right "<<part_left<<" - "<<part_right<<endl;
        it_part_left = NodeId2VecId.find(part_left);
        it_part_right = NodeId2VecId.find(part_right);

        if (it_part_right!=NodeId2VecId.end() and it_part_left==NodeId2VecId.end()) {//right_part is already in a vector
            vec_id = (*it_part_right).second;
            NodeId2VecId.insert(make_pair(part_left, vec_id));
        }

        if (it_part_left!=NodeId2VecId.end() and it_part_right==NodeId2VecId.end()) {//left_part is already in a vector
            vec_id = (*it_part_left).second;
            NodeId2VecId.insert(make_pair(part_right, vec_id));
        }

        if (it_part_left==NodeId2VecId.end() and it_part_right==NodeId2VecId.end()) {//neither is in a vector
            NodeId2VecId.insert(make_pair(part_left, new_vec_id));
            NodeId2VecId.insert(make_pair(part_right, new_vec_id));
            ++new_vec_id;
        }
    }

    vector<list<unsigned int>> L(new_vec_id);//L is a vector of vector, each vector stores parts that should be mapped onto
    //the same block
    for (unordered_map<unsigned int, unsigned int>::iterator it=NodeId2VecId.begin(); it!=NodeId2VecId.end(); ++it) {
        vec_id = (*it).second;
        L[vec_id].push_back((*it).first);
    }

    vector<unsigned int> current_vector;
    unsigned int num_procs_needed = 0;
    int min_w_node_id;
    while (!L.empty()) {
        current_vector.assign(L.back().begin(), L.back().end());
        L.pop_back();

        num_procs_needed = 0;
        for (vector<unsigned int>::iterator it=current_vector.begin(); it!=current_vector.end(); ++it) {
            num_procs_needed+=Qgraph->GetNode(*it)->GetCopies();
        }

        min_w_node_id = -1;
        while (num_procs_needed>Config_Core) {//parts in a vector demand too many cores
            for (vector<unsigned int>::iterator it=current_vector.begin(); it!=current_vector.end(); ++it) {
                if (Qgraph->GetNode(*it)->GetCopies()==3) {
                    min_w_node_id = *it;
                    break;
                }
            }

            for (vector<unsigned int>::iterator it=current_vector.begin(); it!=current_vector.end(); ++it) {
                if (Qgraph->GetNode(*it)->GetCopies()==3) {
                    if (Qgraph->GetNode(*it)->GetWeight()<Qgraph->GetNode(min_w_node_id)->GetWeight()) {
                        min_w_node_id = *it;
                    }
                }
            }

            if(min_w_node_id != -1){
                if (Qgraph->GetNode(min_w_node_id)->GetCopies()==3) {
                    Qgraph->GetNode(min_w_node_id)->SetCopies(1);//reduce requested cores by 2
                    num_procs_needed = num_procs_needed - 2;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }

        if (ProcessorsLeft->at(temp_block-1)<num_procs_needed) {
            temp_block--;
        }

        if (temp_block<1) {
            return false;
        }

        for (vector<unsigned int>::iterator it=current_vector.begin(); it!=current_vector.end(); ++it) {
//            cout<<"current_part "<<*it;
            Sets[temp_block-1].push_back(*it);

            Qgraph->GetNode(*it)->SetBlock(temp_block);
//            Qgraph->GetNode(*it)->SetCore(ProcessorsLeft->at(temp_block-1));
            if (Qgraph->GetNode(*it)->GetCopies()==3) {
//                cout<<", map onto block "<<temp_block<<", core "<<ProcessorsLeft->at(temp_block-1)<<", "<<ProcessorsLeft->at(temp_block-1)-1<<", "<<ProcessorsLeft->at(temp_block-1)-2<<endl;
                Sets[temp_block-1].push_back(*it);
                Sets[temp_block-1].push_back(*it);
                ProcessorsLeft->at(temp_block-1) = ProcessorsLeft->at(temp_block-1)-3;
            } else {
//                cout<<", map onto block "<<temp_block<<", core "<<ProcessorsLeft->at(temp_block-1)<<endl;
                ProcessorsLeft->at(temp_block-1) = ProcessorsLeft->at(temp_block-1)-1;
            }
        }
    }

//    bool success = Map_Topology(graph, BrokenEs, ProcessorsLeft, block, Qgraph);
    bool success = Map_Topology(Qgraph, block, &Sets);

    return success;
}

void PrintNewFormat(Graph* graph, string out){
    vector<unsigned int> new2original;
    unordered_map<unsigned int, unsigned int> original2new;
    DFS(graph, graph->GetSource(), &new2original);
    
    unsigned int new_id=1;
    vector<unsigned int>* predecessors;
    vector<double>::iterator inputedge_size_it;
    unsigned int original_node;
    ofstream OutFile(out.c_str());
    
    OutFile<<"id predecessor_id input_edge_weight node_weight"<<endl;
    OutFile<<1<<" "<<0<<" "<<graph->GetNode(1)->GetInputEdgeSize()->front()<<" "<<graph->GetNode(1)->GetWeight()<<endl;
    while (!new2original.empty()) {
        original_node = new2original.back();
        cout<<"original node "<<original_node<<", new node "<<new_id<<endl;
        new2original.pop_back();
        
        original2new.insert(make_pair(original_node, new_id));
        predecessors = graph->GetNode(original_node)->GetPredecessors();
        inputedge_size_it = graph->GetNode(original_node)->GetInputEdgeSize()->begin();
        
        for (vector<unsigned int>::iterator iter=predecessors->begin(); iter!=predecessors->end(); ++iter,++inputedge_size_it) {
            if (*iter==0) {
                OutFile<<new_id<<" "<<0<<" "<<graph->GetNode(original_node)->GetInputEdgeSize()->front()<<" "<<graph->GetNode(original_node)->GetWeight()<<endl;
//                cout<<new_id<<" "<<0<<" "<<graph->GetNode(original_node)->GetInputEdgeSize()->front()<<" "<<graph->GetNode(original_node)->GetWeight()<<endl;
            } else {
                OutFile<<new_id<<" "<<original2new[(*iter)]<<" "<<*inputedge_size_it<<" "<<graph->GetNode(original_node)->GetWeight()<<endl;
//                cout<<new_id<<" "<<original2new[(*iter)]<<" "<<*inputedge_size_it<<" "<<graph->GetNode(original_node)->GetWeight()<<endl;
            }
        }
        
        new_id++;
    }
    OutFile.close();
}

void ConvertFormat(const char * inputfile, string outdir){
    ifstream OpenFile(inputfile);
    string graphName;
    string newfilename;
    do {
        OpenFile>>graphName;
        cout<<"Transfer format of "<<graphName<<endl;
        Graph* graph = new Graph();
        parse_graph(graphName.c_str(), graph);
        newfilename = outdir;
        PrintNewFormat(graph, newfilename.append(graphName));
        delete graph;
    } while (OpenFile.good());
    OpenFile.close();
}

bool Check(Graph* graph){
    bool no_joinfork_same_node = true;
    vector<Node*>* Nodes=graph->GetNodes();
    
    for (vector<Node*>::iterator it=Nodes->begin(); it!=Nodes->end(); ++it) {
        if ((*it)->GetPredecessors()->size()>1 and (*it)->GetSuccessors()->size()>1) {
            no_joinfork_same_node = false;
            break;
        }
    }
    
    return no_joinfork_same_node;
}

void CheckGraphs(const char * inputfile){
    ifstream OpenFile(inputfile);
    string graphName;
    bool result = true;
    do {
        OpenFile>>graphName;
        Graph* graph = new Graph();
        parse_graph(graphName.c_str(), graph);
        result = Check(graph);
        delete graph;
        cout<<"Check "<<graphName<<", no fork join same node "<<result<<endl;
    } while (OpenFile.good());
    OpenFile.close();
}

void GetStastic(Graph* graph, unsigned int* nbr_edges, double* max_node_w, double* min_node_w, double* max_edge_w, double* min_edge_w, unsigned int* max_degree, double* weight_total_node, double* weight_total_edge){
    vector<Node*>* Nodes = graph->GetNodes();
    
    *nbr_edges = 0;
    *max_node_w = graph->GetSource()->GetWeight();
    *min_node_w = graph->GetSource()->GetWeight();
    *max_edge_w = graph->GetSource()->GetInputEdgeSize()->front();
    *min_edge_w = graph->GetSource()->GetInputEdgeSize()->front();
    *max_degree = graph->GetSource()->GetInputEdgeSize()->size();
    *weight_total_node = 0;
    *weight_total_edge = 0;
    
    vector<double>* inputEdgeSize;
    for (vector<Node*>::iterator currentNode_iter=Nodes->begin(); currentNode_iter!=Nodes->end(); ++currentNode_iter) {
        *nbr_edges+=(*currentNode_iter)->GetInputEdgeSize()->size();
        *weight_total_node+=(*currentNode_iter)->GetWeight();
        
        if ((*currentNode_iter)->GetWeight()>(*max_node_w)) {
            *max_node_w = (*currentNode_iter)->GetWeight();
        }
        
        if (*min_node_w == 0 and (*currentNode_iter)->GetWeight()>0) {
            *min_node_w = (*currentNode_iter)->GetWeight();
        }
        
        if ((*currentNode_iter)->GetWeight()<(*min_node_w) and (*currentNode_iter)->GetWeight()>0){
            *min_node_w = (*currentNode_iter)->GetWeight();
        }
        
        inputEdgeSize = (*currentNode_iter)->GetInputEdgeSize();
        for (vector<double>::iterator edge_size_iter=inputEdgeSize->begin(); edge_size_iter!=inputEdgeSize->end(); ++edge_size_iter) {
            if (*min_edge_w == 0 and (*edge_size_iter)>0) {
                *min_edge_w = *edge_size_iter;
            }
            
            if ((*edge_size_iter)<(*min_edge_w) and (*edge_size_iter)>0) {
                *min_edge_w = *edge_size_iter;
            }
            
            if ((*edge_size_iter)>(*max_edge_w)) {
                *max_edge_w = *edge_size_iter;
            }
            
            *weight_total_edge+=(*edge_size_iter);
        }
        
        if (inputEdgeSize->size()>(*max_degree)) {
            *max_degree = inputEdgeSize->size();
        }
    }

    return;
}

void Stastics(const char * inputfile){
    ifstream OpenFile(inputfile);
    string graphName;

    cout<<"name nbr_nodes nbr_edges max_node_weight min_node_weight max_edge_size min_edge_size max_degree node_weight edge_weight"<<endl;
    
    unsigned int nbr_nodes=0;
    unsigned int nbr_edges=0;
    double max_node_w=0;
    double min_node_w=0;
    double max_edge_w=0;
    double min_edge_w=0;
    unsigned int max_degree=0;
    double total_weight_node=0;
    double total_weight_edge=0;
    
    do {
        OpenFile>>graphName;
        Graph* graph = new Graph();
        parse_graph(graphName.c_str(), graph);
        
        GetStastic(graph, &nbr_edges, &max_node_w, &min_node_w, &max_edge_w, &min_edge_w, &max_degree, &total_weight_node, &total_weight_edge);
        nbr_nodes = graph->GetNodes()->size();
        cout<<graphName<<" "<<nbr_nodes<<" "<<nbr_edges<<" "<<max_node_w<<" "<<min_node_w<<" "<<max_edge_w<<" "<<min_edge_w<<" "<<max_degree<<" "<<total_weight_node<<" "<<total_weight_edge<<endl;
        
        delete graph;
    } while (OpenFile.good());
    OpenFile.close();
}

double GetEnergyCom(Node* part_pred, Node* part_succ){
    double alpha;
    if (part_pred->GetBlock()==part_succ->GetBlock()) {
        alpha = Alpha[0];
    } else {
        alpha = Alpha[1];
    }
    
    double edge_size;
    vector<double>::iterator iter_edge_size = part_succ->GetInputEdgeSize()->begin();
    for (vector<unsigned int>::iterator it=part_succ->GetPredecessors()->begin(); it!=part_succ->GetPredecessors()->end(); ++it, ++iter_edge_size) {
        if (*it == part_pred->GetId()) {
            edge_size = *iter_edge_size;
            break;
        }
    }
    
    double energy = alpha*edge_size;
    
    return energy;
}

double EnergyCost(Graph* Qgraph){
    double energy=0;
    
    queue<Node*> Waiting_Parts;
    Waiting_Parts.push(Qgraph->GetSource());
    Node* current_part;
    vector<unsigned int>* successors;
    double energy_part;
    double output_weight = 0;
    double energy_comm;
    
//    calculate energy cost of running current_part and of sending outfiles to next_part
    while (!Waiting_Parts.empty()) {
        current_part = Waiting_Parts.front();
//        cout<<"calculating energy cost on part "<<current_part->GetId();
        Waiting_Parts.pop();
        successors = current_part->GetSuccessors();
        output_weight = 0;
        for (vector<unsigned int>::iterator it=successors->begin(); it!=successors->end(); ++it) {
            output_weight+=Qgraph->GetNode((*it))->GetInputEdgeSize()->front();
        }
        
        if (current_part->GetCopies()==1) {//energy cost on computing
            energy_part = GetEMaxS(current_part->GetWeight());
        } else if (current_part->GetCopies()==3){
            energy_part = GetETrip(current_part->GetWeight(), output_weight);
        } else {
            cout<<"EnergyCost error 1"<<endl;
        }
        
        for (vector<unsigned int>::iterator it=successors->begin(); it!=successors->end(); ++it) {
            //energy cost on sending output files
            energy_comm = GetEnergyCom(current_part, Qgraph->GetNode(*it));
            energy_part += energy_comm;
        }
//        cout<<" and on sending out, "<<energy_part<<endl;
        
        energy += energy_part;
        
        if (successors->size()==1) {
            if (Qgraph->GetNode(successors->front())->GetPredecessors()->size()>1) {//successor includes a join node
                if (current_part->GetId()==Qgraph->GetNode(successors->front())->GetPredecessors()->front()) {
                    //make sure we only put join node into Waitiong_Parts once
                    Waiting_Parts.push(Qgraph->GetNode(successors->front()));
                }
            } else {
                Waiting_Parts.push(Qgraph->GetNode(successors->front()));
            }
        } else {
            for (vector<unsigned int>::iterator it=successors->begin(); it!=successors->end(); ++it) {
                Waiting_Parts.push(Qgraph->GetNode(*it));
            }
        }
    }
    
    return energy;
}

///check if the map is valid
bool CheckMap(Graph* Qgraph){
    queue<Node*> Waiting_Parts;
    Waiting_Parts.push(Qgraph->GetSource());
    Node* current_part;
    vector<unsigned int>* successors;
    double input_edge_weight = 0;
    bool success = true;
    vector<double>::iterator it_edge_size;
    double beta;
    
    while (!Waiting_Parts.empty()){
        current_part = Waiting_Parts.front();
//        cout<<"check part "<<current_part->GetId()<<endl;
        Waiting_Parts.pop();
        successors = current_part->GetSuccessors();
        
        if (current_part->GetWeight()>SMAX*period) {
            success = false;
//            cout<<", fail because of its size."<<endl;
            break;
        }
        
        if (current_part->GetBlock()<=0 or current_part->GetCore()<=0) {
            success = false;
//            cout<<", fail because of no core assigned."<<endl;
            break;
        }
        
        if(current_part->GetId()!=Qgraph->GetSource()->GetId()){
            it_edge_size = current_part->GetInputEdgeSize()->begin();
            for (vector<unsigned int>::iterator it_node=current_part->GetPredecessors()->begin(); it_node!=current_part->GetPredecessors()->end();
                 ++it_node,++it_edge_size) {
                input_edge_weight = *it_edge_size;
                
                if (current_part->GetBlock()==Qgraph->GetNode(*it_node)->GetBlock()) {
                    if(current_part->GetCore()==Qgraph->GetNode(*it_node)->GetCore()){
                        beta = numeric_limits<double>::max();
                    }
                    beta = BandWidth[0];
                } else {
                    beta = BandWidth[1];
                }
                
                if (input_edge_weight>beta*period) {
                    success = false;
//                    cout<<", fail because of input communication."<<endl;
                    break;
                }
            }
            
            if (success==false) {
                break;
            }
        }
        
        if (successors->size()==1) {
            if (Qgraph->GetNode(successors->front())->GetPredecessors()->size()>1) {//successor includes a join node
                if (current_part->GetId()==Qgraph->GetNode(successors->front())->GetPredecessors()->front()) {
                    //make sure we only put join node into Waitiong_Parts once
                    Waiting_Parts.push(Qgraph->GetNode(successors->front()));
                }
            } else {
                Waiting_Parts.push(Qgraph->GetNode(successors->front()));
            }
        } else {
            for (vector<unsigned int>::iterator it=successors->begin(); it!=successors->end(); ++it) {
                Waiting_Parts.push(Qgraph->GetNode(*it));
            }
        }
    }
    
    return success;
}

void LabelNodesinForkJoin(Graph* graph, unsigned int StartNode){
    stack<unsigned int> waiting_nodes;
    waiting_nodes.push(StartNode);
    Node* currendNode;
    unsigned int end_node_temp;
    int block = graph->GetNode(StartNode)->GetBlock();
    int core = graph->GetNode(StartNode)->GetCore();
    
    while (!waiting_nodes.empty()) {
        currendNode = graph->GetNode(waiting_nodes.top());
        currendNode->SetBlock(block);
        currendNode->SetCore(core);
        currendNode->SetCopies(1);
//        cout<<"label node "<<currendNode->GetId()<<endl;
        
        waiting_nodes.pop();
        for (vector<unsigned int>::iterator iter = currendNode->GetSuccessors()->begin(); iter!=currendNode->GetSuccessors()->end(); ++iter) {
            if (graph->GetNode(*iter)->GetNbrSuccessor()>1) {//this node is a fork node
                GetStructureEndNode(graph, (*iter), &end_node_temp);
                graph->GetNode(*iter)->SetBlock(block);
                graph->GetNode(*iter)->SetCore(core);
                graph->GetNode(*iter)->SetCopies(1);
                LabelNodesinForkJoin(graph, (*iter));
                waiting_nodes.push(end_node_temp);
            }
            
            if (graph->GetNode(*iter)->GetNbrSuccessor()<=1) {
                if (graph->GetNode(*iter)->GetPredecessors()->size()<=1) {
                    waiting_nodes.push((*iter));
                }// else {this node is a join node}
            }
        }
    }
    
    return;
}

//return if it's success
bool NextCore(int* current_core, int* current_block, double* current_workload){
    bool success = true;
    *current_workload = 0;
    
    if (*current_core>1) {
        *current_core = *current_core -1;
    } else if(*current_block>1){
        *current_block = *current_block - 1;
        *current_core = Config_Core;
    } else {
        success = false;
    }
    
    return success;
}

bool Is_successor(Graph* graph, unsigned int predecessor, unsigned int successor){
    bool is_succ = false;
    for (vector<unsigned int>::iterator iter=graph->GetNode(predecessor)->GetSuccessors()->begin();
         iter!=graph->GetNode(predecessor)->GetSuccessors()->end(); ++iter) {
        if((*iter)==successor){
            is_succ = true;
            break;
        }
    }
    return is_succ;
}

//return if mapping is success, but not sure it's valid
//return edges broken
bool MaxSpeed(Graph* graph, unsigned int nbr_blocks, unsigned int nbr_cores, vector<unsigned int>* brokenEs){
    bool success = true;
    vector<unsigned int> vectorL;
    unsigned int current_node;
    int current_core = nbr_cores;
    int current_block = nbr_blocks;
    double workload_on_current_core = 0;
    double node_weight;
    double maximum_load = period*SMAX;
    brokenEs->clear();
    unsigned int end_node;
    bool is_successor=true;
    
    DFS(graph, graph->GetSource(), &vectorL);
    unsigned int node_mapped_before=vectorL.back();
    while (!vectorL.empty()) {
        current_node = vectorL.back();
        vectorL.pop_back();
        
        if (graph->GetNode(current_node)->GetBlock()==0) {//not map yet
//            if (workload_on_current_core!=0) {
            if (node_mapped_before != current_node) {
                is_successor = Is_successor(graph, node_mapped_before, current_node);
                if (is_successor==false or graph->GetNode(node_mapped_before)->GetSuccessors()->size()>1) {
                    //node_mapped_before is a fork node, map this node onto a new core
                    success = NextCore(&current_core, &current_block, &workload_on_current_core);
                    if (success == false) {
                        break;
                    }
                    brokenEs->push_back(current_node);
                }
            }
            
            if (graph->GetNode(current_node)->GetSuccessors()->size()>1) {//this is a fork node
                node_weight = GetStructureEndNode(graph, current_node, &end_node);
                if ((node_weight+workload_on_current_core) > maximum_load){
                    //a whole fork-join is too much for this core
                    
                    node_weight = graph->GetNode(current_node)->GetWeight();
                    if ((node_weight+workload_on_current_core) > maximum_load) {
                        success = NextCore(&current_core, &current_block, &workload_on_current_core);
                        if (success == false) {
                            break;
                        }
                        brokenEs->push_back(current_node);
                    }
                    
                    workload_on_current_core += node_weight;
                    graph->GetNode(current_node)->SetBlock(current_block);
                    graph->GetNode(current_node)->SetCore(current_core);
                    graph->GetNode(current_node)->SetCopies(1);
                    node_mapped_before = current_node;
//                    cout<<"map node "<<current_node<<" onto block "<<current_block<<", core "<<current_core<<endl;
                } else {//we can map the whole fork-join onto this core
                    workload_on_current_core += node_weight;
                    graph->GetNode(current_node)->SetBlock(current_block);
                    graph->GetNode(current_node)->SetCore(current_core);
                    graph->GetNode(current_node)->SetCopies(1);
                    graph->GetNode(end_node)->SetBlock(current_block);
                    graph->GetNode(end_node)->SetCore(current_core);
                    graph->GetNode(end_node)->SetCopies(1);
                    node_mapped_before = end_node;
//                    cout<<"map nodes between "<<current_node<<" and "<<end_node<<" onto block "<<current_block<<", core "<<current_core<<endl;
                    LabelNodesinForkJoin(graph, current_node);
                }
            } else if (graph->GetNode(current_node)->GetPredecessors()->size()>1){//this is a join node
//                if (workload_on_current_core != 0) {
                    success = NextCore(&current_core, &current_block, &workload_on_current_core);
                    if (success==false) {
                        break;
                    }
                    brokenEs->push_back(current_node);
//                }
                node_weight = graph->GetNode(current_node)->GetWeight();
                workload_on_current_core += node_weight;
                graph->GetNode(current_node)->SetBlock(current_block);
                graph->GetNode(current_node)->SetCore(current_core);
                graph->GetNode(current_node)->SetCopies(1);
                node_mapped_before = current_node;
//                cout<<"map node "<<current_node<<" onto block "<<current_block<<", core "<<current_core<<endl;
            } else {// a chain node
                node_weight = graph->GetNode(current_node)->GetWeight();
                if ((node_weight+workload_on_current_core) > maximum_load) {
                    success = NextCore(&current_core, &current_block, &workload_on_current_core);
                    if (success == false) {
                        break;
                    }
                    
                    brokenEs->push_back(current_node);
                }
                
                workload_on_current_core += node_weight;
                graph->GetNode(current_node)->SetBlock(current_block);
                graph->GetNode(current_node)->SetCore(current_core);
                graph->GetNode(current_node)->SetCopies(1);
                node_mapped_before = current_node;
//                cout<<"map node "<<current_node<<" onto block "<<current_block<<", core "<<current_core<<endl;
            }
        }
    }
    
    return success;
}
