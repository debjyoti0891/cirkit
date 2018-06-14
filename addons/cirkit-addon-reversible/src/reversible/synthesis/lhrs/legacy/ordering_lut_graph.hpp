#ifndef ORDER_LUT_GRAPH
#define ORDER_LUT_GRAPH
#include <algorithm>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <queue>
#include <vector>
#include <random>
#include <stack>
#include <set>
 
namespace cirkit {
namespace legacy {

// Create a struct to hold properties for each vertex

struct VertexProperties
{
    int id;
    int topoLabel;
    int mffcLabel;
    int nodeType;
    std::vector<int> successor_vertices;//list of successors
    std::vector<int> predecessor_vertices;
    std::set<int> output_labels; // list of terminal nodes dependent on this
    int qubit_alloc; //store which qubiit currently allocated. else -1
};

class ComparePrio
{
    public: // tuples are (node id, topolevel, contrib)
        bool operator()(std::tuple<int,int,float> n1,std::tuple<int,int,float> n2) {
            if(std::get<1>(n1) == std::get<1>(n2))
                return std::get<2>(n1) < std::get<2>(n2);
            else
                return std::get<1>(n1) < std::get<1>(n2);
        }
};

class ComparePrioTopo
{
    public: // tuples are (node id, topolevel)
        bool operator()(std::tuple<int,int> n1,std::tuple<int,int> n2) {
                return std::get<1>(n1) < std::get<1>(n2);
        }
};

class ComparePrioOutput
{
    public: // tuples are (node id, topolevel)
        bool operator()(std::tuple<int,int> n1,std::tuple<int,int> n2) {
                return std::get<1>(n1) < std::get<1>(n2);
        }
};


typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexProperties> Graph;

typedef boost::graph_traits < Graph >::in_edge_iterator in_edge_iterator;
typedef boost::graph_traits < Graph >::out_edge_iterator out_edge_iterator;

struct LabelProperties 
{
    int node;
    int numberOfNodes;
    std::vector<int> fanIn;
    std::set<int> dependsOn;
    std::set<int> supports;
};

struct Pending
{
    int nodes;
    std::set<int> dependsOn;
};
typedef struct Pending pendingCone;

struct OrderNode
{
    int topoLabel;
    int node;
    int mffcLabel;
        
    bool operator<(const OrderNode& rhs) const
    {
        // lower topoLabel must be processed later
        return topoLabel < rhs.topoLabel; 
    }
};
typedef struct OrderNode prioNode;

class LutGraph
{
    private:
        // LutGraph
        Graph lutGraph;
        int actual_po_count;
         // map from cone root ---label, number of nodes 
        // and fan_in_vector
        std::map<int, struct LabelProperties> labelMap;

        // vector that acts as a stack to store current label
        std::vector<int> labelStack;
        int currentLabel;


        //TODO :  vector  of steps
        // 0,node = compute node 
        // 1,node = uncompute node
        std::vector<std::pair<int,int>> steps;


        std::vector<int> pi, po; 

        // 0,node,qubit no., output_label = compute node 
        // 1,node, qubit no.,output_label = uncompute node
        std::vector<std::tuple<int,int,int,int> > pebble_steps;//
        
        // set of nodes currently computed
        std::set<int> computedSet;

        // set of output labels that have been computed 
        std::set<int> computedOutput;

        // set of nodes that have been computed once
        std::set<int> computedOnce;
        
        //priority queue for pebbling tuples are (node id, topolevel, contrib)
        std::priority_queue<std::tuple<int,int,float>,std::vector<std::tuple<int,int,float>>,ComparePrio> pebble_prio;

        //priority queue for level order sort (node id, topolevel)
        std::priority_queue<std::tuple<int,int>,std::vector<std::tuple<int,int>>,ComparePrioTopo> topo_prio;

        std::priority_queue<std::tuple<int,int>,std::vector<std::tuple<int,int>>,ComparePrioOutput> output_prio;

        std::set<int> dependency_set; //needed in compute_during_uncomp()

        int qubits_available;//indicates total no. of qubits at beginnging
        int outputs_comp_yet;

        int qubits_used_yet;//qubits reqd till now
        int current_qubits;

        std::stack<int> qubits_list;// add and remove qubits from here. 1 at index 
        std::vector<bool> qubits_bool;

        int total_steps;
        int current_steps; // no of steps reqd for current output



        // stores the order of computed labels
        std::vector<int> labelOrdering;
        
        // map to store computation order of individual
        // nodes in a cone
        std::map<int, std::vector<int>> labelCompute;

        std::map<int, std::vector<int>> output_stats; // map output node to 
        //(no.of nodes feeding it, level, qubits reqd to compute it, steps reqd)

        int perfCounter;        
        void labelBFS();
        void topoOrdering(int vertex);
        void fanin(int child); 
        void visitChild(int child);
        void determineSupport(int output);
        void labelAndInit();
        std::set<int> dependencies(int root);
        void countPending(int label, pendingCone& p);
        void computeCone(int coneRoot);
        void compute(int root);

        bool uncomputeCone(int coneRoot, bool output, bool force);
        void uncompute(int root, bool output);

        
        void genSolution(std::vector<int> perm, std::vector<int> &sol, int k, double alpha);
        int evaluateCost(std::vector<int> poOrder);
        void evalRemoveNodes(int nodes, int &allocatedNodes, int &freeNodes, bool output);
        void evalAddNodes(int nodes, int &allocatedNodes, int &freeNodes);

        void preprocess_pebble();
        void pebble_compute(int output_label);
        void pebble_uncompute(int output_label);
        float contrib(int n, int label);
        int get_compute_qubit(int top);
        bool eligible(int node, int output_label);
        bool succDep(int node, int output_label);
        bool succUncomp(int node, int output_label);
        void compute_during_uncomp(int node,int parent,int output_label);
        void find_dep_set(int node);
        void reorder_prio(int output_label);
    public:
        LutGraph(); //default constructor
        LutGraph(std::vector<std::pair<int,int>> edges); // create graph from list of edges
        LutGraph(std::vector<std::pair<int,int>> edges, std::vector<int> output);
        LutGraph(std::vector<std::pair<int,int>> edges, std::vector<int> output, int po_count);
        LutGraph(Graph g); // A boost graph
        //~LutGraph();
        void determineMFFC(); 
        void pebble();
        unsigned int scheduleGraph(int order, bool crit );
        void transitiveMFFC();
        inline std::vector<std::pair<int,int>> getSteps(){ return steps; }
        inline std::vector<std::tuple<int,int,int,int> > getPebbleSteps(){ return pebble_steps; }
        inline int getQubitCount(){ return qubits_used_yet; }
        inline int getQubitsFree(){ return qubits_available - outputs_comp_yet; }
};
}} // end of namespace
#endif
