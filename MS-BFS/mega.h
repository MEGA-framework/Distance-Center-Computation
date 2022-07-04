#include "graph.h"
#include "Timer.h"
#include <fstream>

Graph* pruning(Graph* graph) {
	int M = 1;
	int N = graph->idMap.size();

	std::set<int> deleted;
	bool shouldContinue = true;

	Graph* cp = new Graph(graph);

	while (shouldContinue) {
		shouldContinue = false;

		for (const auto pair : graph->nodeMap) {
			int id = pair.first;
			Node* node = pair.second;

			if (deleted.find(id) == deleted.end()) {
				int deg = node->getDegree();
				if (deg == M && N > 1) {
					deleted.insert(id);
					N --;
					shouldContinue = true;

					auto it = node->neighbours.begin();
					int parentID = *it;
					Node* parent = graph->nodeMap[parentID];
					parent->T += node->T;
					parent->D += node->T + node->D;
					node->parentID = parentID;

					Node* parentCP = cp->nodeMap[parentID];
					parentCP->T += node->T;
					parentCP->D += node->T + node->D;
					Node* nodeCP = cp->nodeMap[id];
					nodeCP->parentID = parentID;

					parent->removeNeighbour(id);
				}
				if (N == 2) {
					shouldContinue = false;
				}
			}
		}
	}

	for (int id : deleted) {
		graph->removeNodeWithID(id);
	}

	return cp;
}

void levelBuilding(Graph* graph, BFSTree* tree) {
	for (int i = 1; i <= tree->levels.size(); i++) {
		Node* seedNode = graph->nodeMap[tree->seedID];
		int levelSum = i * seedNode->T + seedNode->D;

		std::vector<int> jAry;
		jAry.reserve(3);

		for (int j = 1; j <= tree->levels.size(); j++) {
			if (abs(i - j) <= 1) {
				jAry.emplace_back(j);
			}
			else {
				BFSLevel* level = tree->levels[j];
				levelSum += level->levelT * abs(i - j) + level->levelD;
			}
		}

		std::vector<int> jMapKeys;
		std::vector<std::pair<int, int>> jMapValues;

		for (int j : jAry) {
			for (int jNodeID : tree->levels[j]->nodeIDs) {
				if (jNodeID == tree->seedID) continue;
				Node* jNode = graph->nodeMap[jNodeID];
				auto tuple = std::pair<int, int>(jNode->T, jNode->D);
				jMapKeys.emplace_back(jNodeID);
				jMapValues.emplace_back(tuple);
			}
		}


		std::set<int> levelNodeIDs = tree->levels[i]->nodeIDs;
		std::vector<int> idVec(levelNodeIDs.begin(), levelNodeIDs.end());
#pragma omp parallel for
		for (int idx = 0; idx < idVec.size(); idx++) {

			int id = idVec[idx];

			Node* node = graph->nodeMap[id];
			std::set<int> nIDs = node->neighbours;

			int T = 0;
			int D = 0;

			for (int nid : nIDs) {
				if (nid != tree->seedID) {
					Node* nNode = graph->nodeMap[nid];
					T += nNode->T;
					D += nNode->D;
				}
			}

			int lb = levelSum + T + D;

			T = 0;
			D = 0;

			for (int ii = 0; ii < jMapKeys.size(); ii++) {
				int jID = jMapKeys[ii];
				std::pair<int, int> tuple = jMapValues[ii];
				if (nIDs.find(jID) == nIDs.end() && jID != id) {
					T += tuple.first;
					D += tuple.second;
				}
			}

			lb += T * 2 + D + node->D;

			node->setLowerbound(lb);
		}
	}
}

void hierachicalClustering(Graph* graph, int seed, bool tightVersion=true) {
	std::cout << "seed: " << seed << std::endl;

	BFSTree* tree = graph->bfs(seed);

	if (tightVersion) {
		return levelBuilding(graph, tree);
	}

	for (int i = 1; i <= tree->levels.size(); i++) {
		Node* seedNode = graph->nodeMap[seed];
		int levelSum = i * seedNode->T + seedNode->D;

		for (int j = 1; j <= tree->levels.size(); j++) {
			BFSLevel* level = tree->levels[j];
			int multiplier = std::max(1, abs(i - j));
			levelSum += multiplier * level->levelT + level->levelD;
		}
		for (int id : tree->levels[i]->nodeIDs) {
			Node* node = graph->nodeMap[id];
			node->setLowerbound(levelSum - node->T);
		}
	}
}

int backTrack(int candidateID, Graph* graph, Graph* originalGraph) {
	Node* node = originalGraph->nodeMap[candidateID];
	std::map<int, int> backTrackMap;
	for (int nID : node->neighbours) {
		if (graph->idMap.find(nID) == graph->idMap.end() && node->parentID != nID) {
			graph->edgeVisited++;
			Node* nNode = originalGraph->nodeMap[nID];
			int n = originalGraph->idMap.size() - nNode->T;
			if (nNode->T > n) {
				backTrackMap.insert(std::pair<int, int>(nID, nNode->T));
			}
		}
	}
	if (backTrackMap.empty()){
		return candidateID;
	}

	int bestID = 0;
	int bestT = 0;
	for (const auto pair : backTrackMap) {
		if (pair.second > bestT) {
			bestID = pair.first;
			bestT = pair.second;
		}
	}

	return backTrack(bestID, graph, originalGraph);
}

int msComputation(Graph* graph, Graph* originalGraph, int sourcesSize) {
	if (sourcesSize > graph->idMap.size()) {
		printf("source size is larger than the size of the prunned graph\n");
		exit(1);
	}

	std::vector<int> queue(graph->nodeIDs.begin(), graph->nodeIDs.end());

	int candidateID;
	int bfsCount = 1;

	while (true) {
		std::partial_sort(queue.begin(), queue.begin() + sourcesSize, queue.end(), [&graph](const int a, const int b) {
				int d1 = graph->nodeMap[a]->getDistance();
				int d2 = graph->nodeMap[b]->getDistance();
				return d1 < d2;
				});

		std::vector<int> sources(queue.begin(), queue.begin() + sourcesSize);


		int headID = sources[0];
		if (graph->nodeMap[headID]->exact > 0) {
			candidateID = headID;
			goto end;
		}

		std::vector<BFSTree*> forest = graph->msBFS(sources);

		bfsCount++;

		for (BFSTree* tree : forest) {
			levelBuilding(graph, tree);
			Node* seedNode = graph->nodeMap[tree->seedID];

			if (seedNode->exact != seedNode->getLowerbound()) continue;

			// seedNode has tight bound. we check if it's the smallest in queue
			auto minIt = std::min_element(queue.begin(), queue.end(), [&graph](const int a, const int b) {
					int d1 = graph->nodeMap[a]->getDistance();
					int d2 = graph->nodeMap[b]->getDistance();
					return d1 < d2;
					});

			if (*minIt == tree->seedID) {
				// seedNode has tight bound and is the smallest in queue --> it is the center
				candidateID = tree->seedID;
				goto end;
			}
		}

		// delete the forest
		for (BFSTree* tree : forest) {
			for (auto pair : tree->levels) {
				BFSLevel* level = pair.second;
				delete level;
			}
			delete tree;
		}

		printf("BFS count: %d\n", bfsCount);
	}
end:

	printf("BFS#: %d\n", bfsCount);

	if (graph->idMap.size() > (originalGraph->idMap.size()/2)) {
		return candidateID;
	}
	printf("starting back track\n");
	candidateID = backTrack(candidateID, graph, originalGraph);

	return candidateID;
}

int computation(Graph* graph, Graph* originalGraph) {
	std::vector<int> queue(graph->nodeIDs.begin(), graph->nodeIDs.end());

	int candidateID;
	int bfsCount = 1;

	while (true) {
		auto minIt = std::min_element(queue.begin(), queue.end(), [&graph](const int a, const int b){
				int d1 = graph->nodeMap[a]->getDistance();
				int d2 = graph->nodeMap[b]->getDistance();
				return d1 < d2;
				});

		int minID = *minIt;
		queue.erase(minIt);

		Node* minNode = graph->nodeMap[minID];

		if (minNode->exact > 0) {
			candidateID = minID;
			break;
		}

		BFSTree* tree = graph->bfs(minID);
		levelBuilding(graph, tree);
		bfsCount++;

		if (minNode->exact == minNode->getLowerbound()) {
			candidateID = minID;
			break;
		}

		queue.emplace_back(minID);

		// delete the tree
		for (auto pair : tree->levels) {
			BFSLevel* level = pair.second;
			delete level;
		}
		delete tree;

		if (bfsCount % 10 == 0) {
			printf("BFS count: %d\n", bfsCount);
		}
	}

	printf("BFS#: %d\n", bfsCount);

	if (graph->idMap.size() > (originalGraph->idMap.size()/2)) {
		return candidateID;
	}
	printf("starting back track\n");
	candidateID = backTrack(candidateID, graph, originalGraph);

	return candidateID;
}
