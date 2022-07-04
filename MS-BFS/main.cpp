#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "mega.h"
#include "Timer.h"

void start(const char* path, bool useMSBFS) {
	Timer timer;

	std::cout << "loading " << path << std::endl;
	std::cout << "MS-BFS: " << useMSBFS << std::endl;

	timer.start("loading");
	Graph* graph = new Graph(path);
	timer.end("loading");

	std::cout << "N: " << graph->N << std::endl;
	std::cout << "M: " << graph->M << std::endl;

	timer.start("pruning");
	Graph* copy = pruning(graph);
	timer.end("pruning");

	//int seed = graph->getMaxLevelNodeID();
	int seed = graph->getMaxDegreeNodeID(copy);

	timer.start("HC");
	hierachicalClustering(graph, seed);
	timer.end("HC");

	timer.start("computation");
	int centerID;
	if (useMSBFS) {
		centerID = msComputation(graph, copy, std::min(kBitsize, (int)graph->idMap.size()));
	}
	else {
		centerID = computation(graph, copy);
	}
	timer.end("computation");


	std::cout << std::endl << "distance center: " << centerID << std::endl;
	std::cout << "edge visited: " << graph->edgeVisited << std::endl;
	std::cout << "speedup: " << (unsigned long)graph->N * (unsigned long)graph->M / (double)(graph->edgeVisited) << std::endl;

	std::ofstream of;
	of.open("result.csv", std::ios_base::app);
	of << path << ",";
	of << useMSBFS << ",";
	of << kBitsize << ",";
	of << centerID << ",";
	of << graph->edgeVisited << ",";
	of << std::endl;
	of.close();
}

int main(int argc, char** argv) {
	//bool useMSBFS = false;
	//if (argc == 1) {
		//printf("please supply the path to the graph data\n");
		//exit(1);
	//}
	//if (argc > 2) {
		//useMSBFS = atoi(argv[2]) > 0;
	//}

	//start(argv[1], useMSBFS);

	std::vector<std::pair<std::string, bool>> inputs = {
		{"data/loc-gowalla_edges.txt", true},
		{"data/loc-gowalla_edges.txt", false},
		{"data/com-youtube.ungraph.txt", true},
		{"data/com-youtube.ungraph.txt", false},
		{"data/as-skitter.txt", true},
		{"data/as-skitter.txt", false}
	};

	for (auto pair : inputs) {
		const char* path = pair.first.c_str();
		bool flag = pair.second;
		start(path, flag);
	}

	return 0;
}
