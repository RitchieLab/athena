// ScoreHolder.cpp

#include "ScoreHolder.h"

float ScoreHolder::noValue = -1.0;

ScoreHolder::ScoreHolder(){
	maximumCount = 1000000;
	initialize();
}


void ScoreHolder::initialize(){
// cout << "COUNT WAS " << count << endl;
	count = 0;
	holder.clear();
	none.sc = noValue;
}

void ScoreHolder::clear(){
	initialize();
}


///
/// Find and return pointer to node with the score.
///
ScoreHolder::ScoreNode* ScoreHolder::getScore(vector<IndividualTerm*>& terms){

	std::map<IndividualTerm*, ScoreNode>* h = &holder;
	ScoreNode insertNode;

	std::map<IndividualTerm*, ScoreNode>::iterator nodeIter;
	ScoreNode* returnNode;

	for(vector<IndividualTerm*>::iterator iter=terms.begin(); iter != terms.end();
		++iter){
		// if node not found create and insert the node
		if((nodeIter=h->find(*iter)) == h->end()){
			if(count < maximumCount){
// cout << "NOT FOUND insert with score=" << insertNode.sc << endl;
				h->insert(std::pair<IndividualTerm*, ScoreNode>(*iter, insertNode));
				nodeIter = h->find(*iter);
// cout << "score insert check=" << nodeIter->second.sc << endl;
// 				returnNode = &(nodeIter->second);
				count++;
			}
			else{
				// can't return a valid node because reached maximum size
				none.sc = noValue;
				return &none;
			}
		}
// 		if(iter+1 != terms.end())
	 		h = &(nodeIter->second.scores);
	}

// nodeIter=h->find(terms.back());
// if(nodeIter == h->end())
// cout << "NO MATCH " << endl;
// cout << "return score is " << nodeIter->second.sc << endl;
	return &(nodeIter->second);
	// return node 
// 	return &(h->find(terms.back())->second);
}




