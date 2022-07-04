#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <queue>
#include <bitset>
#include <sstream>
#include <string>
#include <fstream>


#ifndef GRAPH_H
#define GRAPH_H

#define kBitsize 64

class Node {
	public:
		int ID;
		int T;
		int D;
		int level;
		int parentID;
		int upperbound;
		int exact;
		std::set<int> neighbours;
		std::set<int> children;

		std::bitset<kBitsize> seen;
		std::bitset<kBitsize> visit;
		std::bitset<kBitsize> visitNext;

		Node(int);
		Node(Node*);
		void addNeighbour(int);
		void removeNeighbour(int);

		int getDegree() { return neighbours.size(); }
		int getDistance();

		void setLowerbound(int lb) {
			if (lb > lowerbound) {
				lowerbound = lb;
			}
		}
		int getLowerbound() { return lowerbound; }

	private:
		int lowerbound;
};

Node::Node(int id) {
	ID = id;
	T = 1;
	D = 0;
	level = -1;
	parentID = -1;
	lowerbound = -1;
	upperbound = -1;
	exact = -1;
	seen = 0;
	visit = 0;
	visitNext = 0;
}

Node::Node(Node* node) {
	ID = node->ID;
	T = node->T;
	D = node->D;
	level = node->level;
	parentID = node->parentID;
	lowerbound = node->lowerbound;
	upperbound = node->upperbound;
	exact = node->exact;
	neighbours = node->neighbours;
	children = node->children;
	seen = node->seen;
	visit = node->visit;
	visitNext = node->visitNext;
}

void Node::addNeighbour(int id) {
	neighbours.insert(id);
}

void Node::removeNeighbour(int id) {
	auto it = neighbours.find(id);
	if (it != neighbours.end()) {
		neighbours.erase(it);
	}
}

int Node::getDistance() {
	if (exact > 0) return exact;
	return lowerbound;
}

class BFSLevel {
	public:
		int levelID;
		int levelT;
		int levelD;
		std::set<int> nodeIDs;

		BFSLevel(int);
		void addNode(Node*);
};

BFSLevel::BFSLevel(int ID) {
	levelID = ID;
	levelT = 0;
	levelD = 0;
}

void BFSLevel::addNode(Node* node) {
	levelT += node->T;
	levelD += node->D;
	nodeIDs.insert(node->ID);
}

class BFSTree {
	public:
		std::map<int, BFSLevel*> levels;
		int seedID;

		BFSTree(int);
};

BFSTree::BFSTree(int id) {
	seedID = id;
}

class Graph {
	public:
		int N;
		int M;

		unsigned long edgeVisited = 0;

		std::map<int, std::set<int>> idMap;
		std::map<int, Node*> nodeMap;
		std::set<int> nodeIDs;

		Graph(const char*);
		Graph(Graph*);
		void removeNodeWithID(int ID);
		BFSTree* bfs(int sourceID);
		int getMaxDegreeNodeID(Graph*);
		int getMaxLevelNodeID();
		int getRandomNodeID();

		std::vector<BFSTree*> msBFS(std::vector<int>);
};

Graph::Graph(const char* file) {
	std::string fpath(file);
	std::ifstream fs(fpath);
	std::string line;

	while (std::getline(fs, line)) {
		if (line[0] == '#') continue;

		std::istringstream iss(line);
		int id1, id2;
		if (!(iss >> id1 >> id2)) {
			printf("failed to parse file %s\n", file);
			exit(1);
		}

		if (id1 == id2) {
			printf("found self edge. skipping\n");
			continue;
		}

		Node* n1;
		Node* n2;

		// add id1 -> id2
		auto mapIt = idMap.find(id1);
		if (mapIt == idMap.end()) {
			std::set<int> s;
			s.insert(id2);
			idMap.insert(std::pair<int, std::set<int>>(id1, s));

			n1 = new Node(id1);
		}
		else {
			mapIt->second.insert(id2);

			n1 = nodeMap[id1];
		}

		// add id2 -> id1
		mapIt = idMap.find(id2);
		if (mapIt == idMap.end()) {
			std::set<int> s;
			s.insert(id1);
			idMap.insert(std::pair<int, std::set<int>>(id2, s));

			n2 = new Node(id2);
		}
		else {
			mapIt->second.insert(id1);

			n2 = nodeMap[id2];
		}

		n1->addNeighbour(id2);
		n2->addNeighbour(id1);
		nodeMap.insert(std::pair<int, Node*>(id1, n1));
		nodeMap.insert(std::pair<int, Node*>(id2, n2));

		nodeIDs.insert(id1);
		nodeIDs.insert(id2);

		M++;
	}

	N = nodeIDs.size();
	M /= 2;
}

Graph::Graph(Graph* graph) {
	M = graph->M;
	N = graph->N;
	edgeVisited = graph->edgeVisited;
	idMap = graph->idMap;
	nodeIDs = graph->nodeIDs;
	for (const auto pair : graph->nodeMap) {
		int id = pair.first;
		Node* node = pair.second;
		Node* nodeCopy = new Node(node);
		nodeMap.insert(std::pair<int, Node*>(id, nodeCopy));
	}
}

void Graph::removeNodeWithID(int ID) {
	std::set<int> neighbours = idMap[ID];
	for (int nID : neighbours) {
		idMap[nID].erase(ID);
		nodeMap[nID]->removeNeighbour(ID);
	}
	Node* node = nodeMap[ID];
	delete node;
	idMap.erase(ID);
	nodeMap.erase(ID);
	nodeIDs.erase(ID);
}

BFSTree* Graph::bfs(int sourceID) {
	BFSTree *tree = new BFSTree(sourceID);

	std::queue<int> queue;
	std::set<int> visited;

	queue.push(sourceID);
	visited.insert(sourceID);

	nodeMap[sourceID]->level = 0;

	int rootExactDistance = nodeMap[sourceID]->D;

	while (!queue.empty()) {
		int ID = queue.front();
		queue.pop();

		for (int nID : idMap[ID]) {
			if (visited.find(nID) != visited.end()) {
				continue;
			}
			queue.push(nID);
			visited.insert(nID);

			edgeVisited ++;

			int level = nodeMap[ID]->level + 1;
			Node* nNode = nodeMap[nID];

			rootExactDistance += nNode->T * level + nNode->D;

			nNode->level = level;
			//nNode->parentID = ID;

			if (tree->levels.find(level) != tree->levels.end()) {
				tree->levels[level]->addNode(nNode);
			}
			else {
				BFSLevel *bfsLevel = new BFSLevel(level);
				bfsLevel->addNode(nNode);
				tree->levels.insert(std::pair<int, BFSLevel*>(level, bfsLevel));
			}
		}
	}

	nodeMap[sourceID]->exact = rootExactDistance;

	return tree;
}


int Graph::getMaxDegreeNodeID(Graph* originalGraph) {
	std::vector<int> idAry;
	idAry.reserve(idMap.size());
	for (const auto pair : idMap) {
		idAry.emplace_back(pair.first);
	}
	return *std::max_element(idAry.begin(), idAry.end(), [&originalGraph](const int lhs, const int rhs){
			int d1 = originalGraph->nodeMap[lhs]->getDegree();
			int d2 = originalGraph->nodeMap[rhs]->getDegree();
			return d1 < d2;
			});
}

int Graph::getMaxLevelNodeID() {
	BFSTree* tree = this->bfs(*nodeIDs.begin());
	std::map<int, BFSLevel*> levels = tree->levels;
	return *(levels[levels.size()]->nodeIDs.begin());
}

int Graph::getRandomNodeID() {
	srand(time(NULL));
	size_t randomIdx = rand() % nodeIDs.size();
	auto it = nodeIDs.begin();
	std::advance(it, randomIdx);
	return *it;
}

std::vector<BFSTree*> Graph::msBFS(std::vector<int> sources) {

	std::vector<int> rootExactDistance;
	std::vector<BFSTree*> forest;

	rootExactDistance.reserve(sources.size());
	forest.reserve(sources.size());

	// reset
	for (int i : nodeIDs) {
		nodeMap[i]->visitNext = 0;
		nodeMap[i]->visit = 0;
		nodeMap[i]->seen = 0;
	}

	// initialize
	for (int i = 0; i < sources.size(); i++) {
		Node* sourceNode = nodeMap[sources[i]];
		sourceNode->seen.set(i);
		sourceNode->visit.set(i);
		rootExactDistance.emplace_back(sourceNode->D);
		forest.emplace_back(new BFSTree(sources[i]));
	}

	int level = 0;

	while (true) {
		bool allSkipped = true;
		level++;

		for (int i : nodeIDs) {
			Node* node = nodeMap[i];

			if (node->visit == 0) continue;

			allSkipped = false;

			for (int nid : node->neighbours) {
				Node* nNode = nodeMap[nid];
				std::bitset<kBitsize> D = node->visit & ~nNode->seen;

				if (D == 0) continue;

				nNode->visitNext = nNode->visitNext | D;
				nNode->seen = nNode->seen | D;

				edgeVisited++;

				int distanceIncre = nNode->T * level + nNode->D;

				for (int j = 0; j < sources.size(); j++) {
					if (!D.test(j)) continue;

					rootExactDistance[j] += distanceIncre;
					BFSTree* tree = forest[j];

					if (tree->levels.find(level) != tree->levels.end()) {
						tree->levels[level]->addNode(nNode);
					}
					else {
						BFSLevel *bfsLevel = new BFSLevel(level);
						bfsLevel->addNode(nNode);
						tree->levels.insert(std::pair<int, BFSLevel *>(level, bfsLevel));
					}
				}


			}
		}

		if (allSkipped) break;

		for (int i : nodeIDs) {
			Node* node = nodeMap[i];

			node->visit = node->visitNext;
			node->visitNext = 0;
		}
	}

	for (int i = 0; i < sources.size(); i++) {
		nodeMap[sources[i]]->exact = rootExactDistance[i];
	}

	return forest;
}

std::ostream &operator<<(std::ostream &os, Graph* const graph) { 
	for (const auto pair : graph->idMap) {
		int id = pair.first;
		std::set<int> neighbours = pair.second;
		os << id << " [";
		for (int nID : neighbours) {
			os << " " << nID << " ";
		}
		os << "]" << std::endl;
	}
	return os;
}

std::ostream &operator<<(std::ostream &os, BFSTree* const tree) {
	for (const auto pair : tree->levels) {
		int level = pair.first;
		BFSLevel* levelObj = pair.second;
		os << "level " << level << ": [";
		for (int id : levelObj->nodeIDs) {
			os << " " << id << " ";
		}
		os << "]" << std::endl;
	}
	return os;
}

#endif

