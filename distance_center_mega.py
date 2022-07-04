class Node:
    def __init__(self, nodeId, dis):
        self.nodeId = nodeId
        self.dis = dis
        self.exact = False
        
    def updateLower(self, lowerBound):
        if lowerBound > self.dis:
            self.dis = lowerBound
        
    def updateDis(self, dis):
        self.dis = dis
        self.exact = True
    
    def isExact(self):
        return self.exact
    
    def returnNodeId(self):
        return self.nodeId
    
    def returnDis(self):
        return self.dis

class PriorityQueue(object): 
    def __init__(self):
        self.queue = [] 
  
    def __str__(self):
        return ' '.join([str(i) for i in self.queue]) 
  
    # for checking if the queue is empty 
    def isEmpty(self):
        return len(self.queue) == [] 
  
    # for inserting an element in the queue 
    def insert(self, data):
        self.queue.append(data) 
  
    # for popping an element based on Priority 
    def extractMin(self):
        try: 
            minS = 0
            for i in range(len(self.queue)): 
                if self.queue[i].returnDis() < self.queue[minS].returnDis(): 
                    minS = i 
            item = self.queue[minS] 
            del self.queue[minS]
            return item 
        except IndexError: 
            print() 
            exit() 
            
def getMinDegreeV(G):
    #node_degree_dict = {}
    minDeg = 99999999
    minDegV = -1
    for NI in G.Nodes():
        if NI.GetDeg() < minDeg:
            minDegV = NI.GetId()
    return minDegV

def calLowerBound(G):
    global nodesList
    seed = G.GetRndNId()
    #print(seed)
    #seed = 107
    
    #find the min degree vertex
    seed = getMinDegreeV(G)
    largestDisFromSeed = snap.GetNodeEcc(G,seed, False)
    NodeNum, NodeVec = G.GetNodesAtHop(seed, largestDisFromSeed, False)
    seed = NodeVec[0]
    
    #BfsTree = snap.GetBfsTree(G, seed, False, False)
    myQueueLower = PriorityQueue()
    myQueueUpper = PriorityQueue()
    C = {}
    j = 1
    Sv = 0
    while True:
        s = snap.TIntV()
        snap.GetNodesAtHop(G, seed, j, s, True)
        if len(s) != 0:
            C[j] = s
            Sv = Sv + len(s) * j
            j = j + 1
        else:
            break
    root = Node(seed, Sv)
    nodesList[seed] = root
    myQueueLower.insert(root)
    '''print 0, " ", seed
    for cluster in C:
        for nodeId in C[cluster]:
            print cluster, " ", nodeId, " ", C[cluster].Len()'''
    maxL = len(C)
    for i in range(1, maxL + 1):
        sumLevel = i
        for j in range(1, maxL + 1):
            if abs(i - j) < 1 or abs(i - j) == 1:
                sumLevel = sumLevel + 2 * C[j].Len()
            else:
                sumLevel = sumLevel + abs(i - j) * C[j].Len()
        for nodeId in C[i]:
            deg = G.GetNI(nodeId).GetDeg()
            lowerBound = sumLevel - deg - 2
            currNode = Node(nodeId, lowerBound)
            nodesList[nodeId] = currNode
            myQueueLower.insert(currNode)
    return myQueueLower

def ExactTotalDis(G, node):
    #print "src ", node.returnNodeId()
    j = 1
    Sv = 0
    while True:
        s = snap.TIntV()
        snap.GetNodesAtHop(G, node.returnNodeId(), j, s, True)
        if len(s) != 0:
            Sv = Sv + len(s) * j
            global edgeVisited
            edgeVisited = edgeVisited + len(s)
            print(len(s))
            j = j + 1
        else:
            break
    return Sv

def UpdateLowerBound(seed, G):
    global nodesList
    C = {}
    j = 1
    Sv = 0
    while True:
        s = snap.TIntV()
        snap.GetNodesAtHop(G, seed, j, s, True)
        if len(s) != 0:
            C[j] = s
            Sv = Sv + len(s) * j
            global edgeVisited
            edgeVisited = edgeVisited + len(s)
            j = j + 1
        else:
            break
    nodesList[seed].updateDis(Sv)
    maxL = len(C)
    for i in range(1, maxL + 1):
        sumLevel = i
        for j in range(1, maxL + 1):
            if abs(i - j) < 1 or abs(i - j) == 1:
                sumLevel = sumLevel + 2 * C[j].Len()
            else:
                sumLevel = sumLevel + abs(i - j) * C[j].Len()
        for nodeId in C[i]:
            deg = G.GetNI(nodeId).GetDeg()
            lowerBound = sumLevel - deg - 2
            nodesList[nodeId].updateLower(lowerBound)

def Pruning(Gp):
    G = snap.ConvertGraph(type(Gp), Gp)
    M = 1
    N = Gp.GetNodes()
    T = {}
    D = {}
    for node in Gp.Nodes():
        nodeId = node.GetId()
        T[nodeId] = 1
        D[nodeId] = 0
        
    continuePruning = True    
    while continuePruning:
        continuePruning = False
        for node in Gp.Nodes():
            nodeId = node.GetId()
            deg = node.GetDeg()
            if deg == M and N > 1:
                parent = node.GetNbrNId(0)
                T[parent] = T[parent] + T[nodeId]
                D[parent] = D[parent] + T[nodeId] + D[nodeId]
                continuePruning = True
                Gp.DelNode(nodeId)
                N = N - 1
    #print(T)
    #print(D)

import pandas as pd
import snap
import time

edgeVisited = 0
nodesList = {}

if __name__ == '__main__':
    graphDF = pd.read_csv('tweets_graph.csv', encoding = "latin-1")

    allNamesList = []     
    for index, row in graphDF.iterrows():
        src = row['Src']
        dst = row['Dst']
        allNamesList.append(src)
        allNamesList.append(dst)

    allNamesList = list(set(allNamesList))

    indexToNameDict = {}
    for index, val in enumerate(allNamesList):
        indexToNameDict[index] = val

    NameToIndexDict = {}
    for index, val in enumerate(allNamesList):
        NameToIndexDict[val] = index

    srcList = []
    complete_G = snap.TUNGraph.New()
    for index, row in graphDF.iterrows():
        src = row['Src']
        dst = row['Dst']
        src = NameToIndexDict[src]
        dst = NameToIndexDict[dst]
        if not complete_G.IsNode(src):
            complete_G.AddNode(src)
        if not complete_G.IsNode(dst):
            complete_G.AddNode(dst)
        complete_G.AddEdge(src, dst)
        srcList.append(src)

    for NI in complete_G.Nodes():
        NID = NI.GetId()
        if complete_G.IsEdge(NID, NID):
            complete_G.DelEdge(NID, NID)
    
    oldV = complete_G.GetNodes()
    oldE = complete_G.GetEdges()

    max_cc_G = complete_G.GetMxScc()
    ccOldV = max_cc_G.GetNodes()
    ccOldE = max_cc_G.GetEdges()
            
    tStart = time.time()
    
    Pruning(max_cc_G)
    newV = max_cc_G.GetNodes()
    newE = max_cc_G.GetEdges()

    queue = calLowerBound(max_cc_G)
    numBFS = 1
    while True:
        minNode = queue.extractMin()
        if minNode.isExact() == True:
            print(minNode.returnNodeId(), " ", minNode.returnDis())
            break
        else:
            lowerBound = minNode.returnDis()
            UpdateLowerBound(minNode.returnNodeId(), max_cc_G)
            numBFS = numBFS + 1
            #print(lowerBound, " ", minNode.returnDis())
            if minNode.returnDis() == lowerBound:
                print(minNode.returnNodeId(), " ", minNode.returnDis())
                break
            else:
                queue.insert(minNode)
    tEnd = time.time()
    speedup = (ccOldV * ccOldE)/(edgeVisited + (ccOldV - newV))
    print(tEnd - tStart)
    print(speedup)