#include "CycleBreaker.h"
#include <iostream>
using namespace std;

/// return number of arcs broken
int CycleBreaker::breakCycles(Athena2DArrayGenome<int>& genome){

// cout << "break cycles" << endl;
// test genome
// 1 4 -1
// -1 -1 -1
// 7 -1 -1
// -1 -1 -1
// 0 -1 -1
// -1 -1 9
// -1 -1 -1
// -1 -1 -1
// -1 -1 -1
// 5 6 -1
// for(int y=0;y<genome.height();y++){
// for(int x=0;x<genome.width(); x++){
// genome.gene(x,y,-1);
// }
// }
// genome.gene(0,0,1);
// genome.gene(1,0,4);
// genome.gene(0,2,7);
// genome.gene(0,4,0);
// genome.gene(2,5,9);
// genome.gene(0,9,5);
// genome.gene(1,9,6);
// cout << "\n";
// genome.write(cout);
// exit(1);




	// set up all the connections in a set
	map<int,node> network;
	int height = genome.height();
	int width = genome.width();
	node blank;
	for(size_t y=0; y<height; y++){
		network[y]=blank;
	}

	for(size_t y=0; y<height; y++){
		for(size_t x=0; x<width; x++){
			// lists all children for one node
			if(genome.gene(x,y) != -1){
				network[genome.gene(x,y)].children.insert(y);
				network[y].parents.insert(genome.gene(x,y));
			}
		}
	}

	map<int,node> origNet = network;

// for(size_t i=0; i<height; i++){
// // 	cout << "node=" << i;
// // 	cout << " parents=";
// // 	for(size_t j=0; j<network[i].parents.size(); j++){
// 	for(set<int>::iterator iter=network[i].parents.begin();
// 		iter != network[i].parents.end(); ++iter){
// 		cout << *iter << " ";
// 	}
// // 	cout << "children=";
// 	for(set<int>::iterator iter=network[i].children.begin();
// 		iter != network[i].children.end(); ++iter){
// 		cout << *iter << " ";
// 	}
// 	cout << endl;
// }
// exit(1);


	vector<int> s1;
	deque<int> s2;
	bool found;
	while(!network.empty()){
// 	cout << "network size=" << network.size() << endl;
		do{
			found = removeSink(network,s2);
		}while(found==true);
		do{
			found = removeSources(network,s1);
		}while(found==true);
		if(!network.empty())
			removeMaxSig(network,s1);
	}

// 	for(size_t i=0; i<s1.size(); i++){
// 		cout << s1[i] << " ";
// 	}
// 	for(size_t i=0; i<s2.size(); i++){
// 		cout << s2[i] << " ";
// 	}
// 	cout<< endl;

	// iterate through vector sequence
	// remove any left-ward edge
	int broken=0;
	set<int> leftSeq;
	for(vector<int>::iterator iter=s1.begin(); iter != s1.end(); ++iter){
		for(set<int>::iterator childIter = origNet[*iter].children.begin();
			childIter != origNet[*iter].children.end(); ++childIter){
			if(leftSeq.find(*childIter) != leftSeq.end()){
// cout << "break " << *iter << " to " << *childIter << endl;
				// matrix[*iter][*childIter]=0;
// 				genome.gene(*iter,*childIter,0);
// 				genome.gene(*childIter, *iter,-1);
				for(int x=0;x<width;x++){
					if(genome.gene(x, *childIter) == *iter){
						genome.gene(x, *childIter,-1);
						break;
					}
				}
				broken++;
			}
		}
		leftSeq.insert(*iter);
	}
	for(deque<int>::iterator iter=s2.begin(); iter != s2.end(); ++iter){
		for(set<int>::iterator childIter = origNet[*iter].children.begin();
			childIter != origNet[*iter].children.end(); ++childIter){
			if(leftSeq.find(*childIter) != leftSeq.end()){
// cout << "break " << *iter << " to " << *childIter << endl;
// 				matrix[*iter][*childIter]=0;
// 				genome.gene(*iter,*childIter,0);
				for(int x=0;x<width;x++){
					if(genome.gene(x, *childIter) == *iter){
						genome.gene(x, *childIter,-1);
						break;
					}
				}
				broken++;
			}
		}
		leftSeq.insert(*iter);
	}
// genome.write(cout);
// exit(1);
	return broken;
}


///
/// remove sink
/// @param
/// @return true when removed one
///
bool CycleBreaker::removeSink(map<int,node>& network, deque<int>& s2){

	for(map<int, node>::iterator iter=network.begin(); iter != network.end(); ++iter){
		if(iter->second.children.empty()){
      s2.push_front(iter->first);
      // remove parents connection to this node
			for(set<int>::iterator parentIter=iter->second.parents.begin();
				parentIter != iter->second.parents.end(); ++parentIter){
				network[*parentIter].children.erase(iter->first);
			}
// cout << "sink remove " << iter->first << endl;
			network.erase(iter);
			return true;
		}
	}
	return false;
}


///
/// remove source
/// @param
/// @return true when removed one
///
bool CycleBreaker::removeSources(map<int,node>& network, vector<int>& s1){

	for(map<int, node>::iterator iter=network.begin(); iter != network.end(); ++iter){
		if(iter->second.parents.empty()){
      s1.push_back(iter->first);
      // remove parents connection to this node
			for(set<int>::iterator childIter=iter->second.children.begin();
				childIter != iter->second.children.end(); ++childIter){
				network[*childIter].parents.erase(iter->first);
			}
			network.erase(iter);
// cout << "source remove " << iter->first << endl;
			return true;
		}
	}
	return false;
}

///
/// Remove node with maximum sigma between input and output connections
/// @param network
/// @param s1 vector to concatenate removed node to
///
void CycleBreaker::removeMaxSig(map<int,node>& network, vector<int>& s1){

	int maxSigma=-1000;
	int maxNode, sigma;
	for(map<int, node>::iterator iter=network.begin(); iter != network.end(); ++iter){
		sigma = iter->second.children.size()-iter->second.parents.size();
		if(sigma > maxSigma){
			maxNode = iter->first;
			maxSigma = sigma;
		}
	}

	// update network
	for(set<int>::iterator childIter=network[maxNode].children.begin();
		childIter != network[maxNode].children.end(); ++childIter){
		network[*childIter].parents.erase(maxNode);
	}
	for(set<int>::iterator parentIter=network[maxNode].parents.begin();
		parentIter != network[maxNode].parents.end(); ++parentIter){
		network[*parentIter].children.erase(maxNode);
	}
	network.erase(maxNode);
	s1.push_back(maxNode);
// cout << "maxSigma=" << maxSigma << " maxNode=" << maxNode << endl;
}

