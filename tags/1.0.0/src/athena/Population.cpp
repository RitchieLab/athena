#include "Population.h"
#include "Terminals.h"

int floatLT::order_adjust=0;

///
/// Assignment operator
///
Population& Population::operator=(const Population& other){
    if(this != &other){
      copy(other);   
    }
    return *this;
}


///
/// copy population
///
void Population::copy(const Population& other){
    SolutionTreeNode* current_node;
    for(current_node = other.solutions.GetFirst(); current_node != other.solutions.GetLast();
      current_node = current_node->GetNext()){
        // need to create a clone operator  
        Solution* sol = current_node->GetData()->clone();
        solutions.Insert(current_node->GetKey(), sol);
    }
    currnode = solutions.GetFirst();
    convertScores = other.convertScores;
}


///
/// Outputs fitness of the solutions in tree in order
/// 
void Population::output_tree(){
    SolutionTreeNode* current_node;
    for(current_node = solutions.GetFirst(); current_node != solutions.GetLast();
      current_node = current_node->GetNext()){
        // need to create a clone operator 
       cout << "tree solution fitness=" << 
           current_node->GetData()->fitness() << endl;
    }
}


///
/// Copy constructor
///
Population::Population(const Population& p){
    copy(p);
}


///
/// Destructor frees any Solutions
///
Population::~Population(){
    clear();
}


///
/// Clears the tree and deletes the Solution pointers
///
void Population::clear(){
    
    if(solutions.GetCount() == 0)
        return;

    SolutionTreeNode* current_node;
    for(current_node = solutions.GetFirst(); current_node != solutions.GetLast();
      current_node = current_node->GetNext()){
        delete current_node->GetData();
    }
    delete current_node->GetData();
    solutions.Clear();
}


///
/// converts scores of population
/// @param train Dataset
/// @param test Dataset
///
void Population::convert_scores(Dataset* train, Dataset* test){

    SolutionTreeNode* node;
    for(node = solutions.GetFirst(); node != solutions.GetLast();
      node = node->GetNext()){
      node->GetData()->adjust_score_out(train, test);
    }
}

///
/// converts scores of population
/// @param train Dataset
///
void Population::convert_scores(Dataset* train){

    SolutionTreeNode* node;
    for(node = solutions.GetFirst(); node != solutions.GetLast();
      node = node->GetNext()){
      
      node->GetData()->adjust_score_out(train);
    }

}




