#include "ordering_lut_graph.hpp"
#undef ORDERING_DEBUG 
#include <fstream> 
namespace cirkit {
namespace legacy {
LutGraph::LutGraph(Graph g)
{ // A boost graph
	lutGraph = g;
	currentLabel = 0;
	std::cout << "Lut Graph size : " << num_vertices(lutGraph)<<"\n";
}

LutGraph::LutGraph(std::vector<std::pair<int,int>> edges)
{
	// the vertex must be from 0 to n.
	for(auto edge: edges)
	{
		boost::add_edge((int)edge.first, (int)edge.second, lutGraph);
	}
	currentLabel = 0;

	std::cout << "Lut Graph size : " << num_vertices(lutGraph) << "\n";
}

LutGraph::LutGraph(std::vector<std::pair<int,int>> edges, std::vector<int> output)
{
    // the vertex must be from 0 to n.
    for(auto edge: edges)
    {
        boost::add_edge((int)edge.first, (int)edge.second, lutGraph);
    }
    currentLabel = 0;
    #ifdef ORDERING_DEBUG
    std::cout<<"Orig out:" ;
    #endif
    for(auto out : output)
    {
        po.push_back(out) ;
#ifdef ORDERING__DEBUG 
        std::cout << out << " ";
    #endif
    }
    std::cout<< "\n";
    std::cout << "Lut Graph size : " << num_vertices(lutGraph) << "\n";
    std::cout << "Number of outputs: " << output.size() << "\n";
    
}

LutGraph::LutGraph(std::vector<std::pair<int,int>> edges, std::vector<int> output, int po_count)
{
    // the vertex must be from 0 to n.
    for(auto edge: edges)
    {
        boost::add_edge((int)edge.first, (int)edge.second, lutGraph);
    }
    currentLabel = 0;
    #ifdef ORDERING_DEBUG
    std::cout<<"Orig out:" ;
    #endif
    for(auto out : output)
    {
        po.push_back(out) ;
#ifdef ORDERING__DEBUG 
        std::cout << out << " ";
    #endif
    }
    actual_po_count = po_count;  
    std::cout<< "\n";
    std::cout << "Lut Graph size : " << num_vertices(lutGraph) << "\n";
    std::cout << "Number of outputs: " << output.size() << "\n";
    

}
LutGraph::LutGraph()
{
	
}
/*
LutGraph::~LutGraph()
{
    pebble_steps.clear();
    pebble_steps.shrink_to_fit();
} */
// void LutGraph::topoOrdering(int v)
// {
//     //std::cout << "ordering : " << v <<":" << lutGraph[v].topoLabel << "\n";
//     std::pair<in_edge_iterator, in_edge_iterator> inEdges = in_edges(v, lutGraph);

//     for(in_edge_iterator iter = inEdges.first; iter != inEdges.second; ++iter)
//     {
		
//         int s = source(*iter, lutGraph);
//         if(lutGraph[s].topoLabel  < lutGraph[v].topoLabel+1)
//         {
//             lutGraph[s].topoLabel = lutGraph[v].topoLabel+1;
//             topoOrdering(s);
//         }
//         else
//             continue;
//     }
// }

void LutGraph::topoOrdering(int v)
{
	//std::cout << "ordering : " << v <<":" << lutGraph[v].topoLabel << "\n";
	std::pair<out_edge_iterator, out_edge_iterator> outEdges = out_edges(v, lutGraph);

	for(out_edge_iterator iter = outEdges.first; iter != outEdges.second; ++iter)
	{
		
		int t = target(*iter, lutGraph);
		if(lutGraph[t].topoLabel  < lutGraph[v].topoLabel+1)
		{
			lutGraph[t].topoLabel = lutGraph[v].topoLabel+1;
			topoOrdering(t);
		}
		else
			continue;
	}
}

void LutGraph::transitiveMFFC()                                                                                                                               
{
	std::set<int> added;
	std::vector<int> toProcess;                                                
	for(auto x : labelMap)                                                     
	{
		toProcess = std::vector<int>();
		int label = x.first;  
		for(auto d : labelMap[label].dependsOn)                                
		{
			toProcess.push_back(d);                                                  
		}

		while(!toProcess.empty())                                              
		{
			int l = toProcess[0];                                              
			toProcess.erase(toProcess.begin());                                                
			for(auto d : labelMap[l].dependsOn)                                
			{
				labelMap[label].dependsOn.insert(d);                                     
				if(added.find(d) == added.end())                               
				{ 
					toProcess.push_back(d);                                          
				}
			}                                                                  
			added.insert(l);                                                   

		std::cout << label <<  "depend : " << l << "\n";
		}
	} 
	for(auto const &ent : labelMap)
	{
		fanin(ent.second.node); // determine fanins
		#ifdef ORDERING_DEBUG

		std::cout << ent.first << ": " << ent.second.node << "," << ent.second.numberOfNodes << "d<-[ ";

		for(auto const &ele : ent.second.dependsOn)
		 std::cout << ele << " ";
		std::cout<< "]    s->[";

		
		for(auto const &ele : ent.second.supports)
		{
			std::cout<< ele << " " ;
		}
	   std::cout << "] [";
		for(auto const &ele : ent.second.fanIn)
			std::cout << ele << " ";

		std::cout << "]\n";
		#endif

	}
} 
void LutGraph::labelAndInit()
{
	// po
	Graph::vertex_iterator vIt, vItEnd;
	boost::tie(vIt,vItEnd) = vertices(lutGraph);
	int idCount = 0;
	for(; vIt != vItEnd; vIt++)
	{
		lutGraph[*vIt].topoLabel = -1;
		lutGraph[*vIt].mffcLabel = 0;
		//  std::cout << "vertex " << *vIt << "\n";
		// TODO : check if pi == po posssible
		if(boost::out_degree(*vIt,lutGraph) == 0)
		{    
			po.push_back(*vIt);
			lutGraph[*vIt].nodeType = 1;
			lutGraph[*vIt].topoLabel = 0;
		}
		else if(boost::in_degree(*vIt,lutGraph) == 0)
		{    
			pi.push_back(*vIt);
			lutGraph[*vIt].nodeType = 2;
		}
		else
			lutGraph[*vIt].nodeType = 3;
		lutGraph[*vIt].id = idCount++;
		#ifdef ORDERING_DEBUG
		std::cout<< lutGraph[*vIt].mffcLabel << " " << lutGraph[*vIt].topoLabel<< "\n"; 
		#endif
	}    

	// determine topo label
	for(auto output: po)      
		topoOrdering(output);
	
	// determine mffcLabel
	labelBFS();
	/*
	for(auto output: po)
	{ 
		perfCounter = 0; 
		std::cout << "Processing Cone PO: " << output << ":l" << lutGraph[output].mffcLabel << "\n";


		if(lutGraph[output].mffcLabel == 0)// new output
			visitChild(output);        
		else
			std::cout << "Labelled output." << output << " TODO\n";
	}
	*/
	boost::tie(vIt,vItEnd) = vertices(lutGraph);
   // update the node count per label 
	for(; vIt != vItEnd; vIt++)
	{
		#ifdef ORDERING_DEBUG

		std::cout << "vertex " << *vIt << " Label : " << lutGraph[*vIt].mffcLabel<< "TopoLabel:" << lutGraph[*vIt].topoLabel << "\n";
		#endif

		int currLabel = lutGraph[*vIt].mffcLabel;
		labelMap[currLabel].numberOfNodes = labelMap[currLabel].numberOfNodes+1;
	}
	// update the support 
	std::cout << "Updating support \n";

	for(auto const output: po)
	{
		determineSupport(output);
	}
	
	std::cout << "Determine fanin \n";
	for(auto const &ent : labelMap)
	{
		fanin(ent.second.node); // determine fanins
		#ifdef ORDERING_DEBUG

		std::cout << ent.first << ": " << ent.second.node << "," << ent.second.numberOfNodes << "d<-[ ";

		for(auto const &ele : ent.second.dependsOn)
		 std::cout << ele << " ";
		std::cout<< "]    s->[";

		
		for(auto const &ele : ent.second.supports)
		{
			std::cout<< ele << " " ;
		}
	   std::cout << "] [";
		for(auto const &ele : ent.second.fanIn)
			std::cout << ele << " ";

		std::cout << "]\n";
		#endif

	}
	
	

}    

void LutGraph::fanin(int child)
{
	std::vector<int> fanIn;
	std::vector<int> toProcess;
	std::set<int> queued;
	int label = lutGraph[child].mffcLabel;
	queued.insert(child);
	toProcess.push_back(child);
	bool isFanIn;
	while(!toProcess.empty())
	{
		int node = toProcess.front();
		toProcess.erase(toProcess.begin());
		std::pair<in_edge_iterator, in_edge_iterator> inEdges = in_edges(node, lutGraph);
		isFanIn = true;
		for(in_edge_iterator iter = inEdges.first; iter != inEdges.second; ++iter)
		{
			// std::cout << *iter << " ";
			int s = source(*iter, lutGraph);
			if(lutGraph[s].mffcLabel == label)
			{   
				if(queued.find(s) == queued.end())
				{
					toProcess.push_back(s);
					queued.insert(s);
				}
				isFanIn = false;
			}
		}   
		if(isFanIn)
			fanIn.push_back(node);
	}
	labelMap[label].fanIn = fanIn;   
	/*
	int fanIn = true;
	if(boost::in_degree(child,lutGraph) == 0 && fanIn)
	{
		labelMap[lutGraph[child].mffcLabel].fanIn.push_back(child);
		return;
	}
	std::pair<in_edge_iterator, in_edge_iterator> inEdges = in_edges(child, lutGraph);

	for(in_edge_iterator iter = inEdges.first; iter != inEdges.second; ++iter)
	{
		// std::cout << *iter << " ";
		int s = source(*iter, lutGraph);
		if(lutGraph[child].mffcLabel == lutGraph[s].mffcLabel)
		{   fanin(s);
			fanIn = false;
		}
	}   
	if(fanIn)
		labelMap[lutGraph[child].mffcLabel].fanIn.push_back(child);
	*/
}
void LutGraph::determineSupport(int output)
{
	int outlabel = lutGraph[output].mffcLabel;
	std::set<int> depends;
	for(auto label: labelMap[outlabel].dependsOn)
	{
		labelMap[label].supports.insert(outlabel);
		for(auto dlabel : labelMap[label].dependsOn)
		{
			depends.insert(dlabel);
		}

	}
	while(!depends.empty())
	{
		int label = *depends.begin();
		depends.erase(depends.begin());
		labelMap[label].supports.insert(outlabel);
		for(auto dlabel : labelMap[label].dependsOn)
		{
			depends.insert(dlabel);
		}


	}

}

void LutGraph::labelBFS()
{
	std::set<int> queued;
	std::priority_queue<prioNode> q;
	int currentLabel =  0;
	for(auto output : po)
	{   
			currentLabel++;
			LabelProperties newProp = {output, 0, std::vector<int>(), std::set<int>()};
			labelMap[currentLabel] = newProp; 
			queued.insert(output);
		lutGraph[output].mffcLabel = currentLabel;
			// add the predecessors
	  std::pair<in_edge_iterator, in_edge_iterator> inEdges = boost::in_edges(output, lutGraph);
		for(; inEdges.first != inEdges.second; ++inEdges.first)
		{
			#ifdef ORDERING_DEBUG
			std::cout << "Processing child : " << *inEdges.first << "\n";
			#endif
			int vertex = source(*inEdges.first,lutGraph);
			if(queued.find(vertex) == queued.end())
			{
				queued.insert(vertex);
				prioNode dNode = {-1*lutGraph[vertex].topoLabel,vertex,0 };
				q.push(dNode);               
			}    
		 }

	} // completed adding the predecessors of the output
	while(!q.empty())
	{
	   prioNode pnode = q.top();
	   q.pop();
	   #ifdef ORDERING_DEBUG
	   std::cout << "Labelling : " << pnode.node << "\n";
		#endif
	   // determine the label based on successors
	   std::pair<out_edge_iterator, out_edge_iterator> outEdges = boost::out_edges( pnode.node, lutGraph);
	   std::set<int> succLabelSet;
	   int existingLabel;
	   for(; outEdges.first != outEdges.second; ++outEdges.first)
	   {
		   int vertex = target(*outEdges.first, lutGraph);
		   existingLabel = (lutGraph[vertex].mffcLabel);
		   succLabelSet.insert(existingLabel);
	   }
	   
	   // check if more than one label ---> conflict
	   if(succLabelSet.size() == 1)
		   lutGraph[pnode.node].mffcLabel = existingLabel;
	   else
	   {
		   currentLabel++;
		   lutGraph[pnode.node].mffcLabel = currentLabel;
		   for(int label : succLabelSet)
		   {
			   labelMap[label].dependsOn.insert(currentLabel);
		   }
		   // add a new label to LabelMap
		   LabelProperties newProp = {pnode.node, 0, std::vector<int>(), std::set<int>()};
		   labelMap[currentLabel] = newProp; 
	   }
	   // add the predecessors
	  std::pair<in_edge_iterator, in_edge_iterator> inEdges = boost::in_edges( pnode.node, lutGraph);
		for(; inEdges.first != inEdges.second; ++inEdges.first)
		{
			#ifdef ORDERING_DEBUG
			std::cout << "Processing child : " << *inEdges.first << "\n";
			#endif
			int vertex = source(*inEdges.first,lutGraph);
			if(queued.find(vertex) == queued.end())
			{
				queued.insert(vertex);
				prioNode dNode = {-1*lutGraph[vertex].topoLabel,vertex,0 };
				q.push(dNode);               
			}    
		 }
	}   // end of while 
}

void LutGraph::visitChild(int child)
{
	perfCounter ++;
	if(perfCounter % 500 == 0)
		std::cout << "visited " << perfCounter << "nodes \n";
	int popLabel = false;
	int childLabel = lutGraph[child].mffcLabel;
	if(childLabel == 0) // no label
	{
		if(labelStack.empty()) 
		{
			currentLabel = currentLabel + 1;
			labelStack.push_back(currentLabel);
			popLabel = true;
			LabelProperties newProp = {child, 0, std::vector<int>(), std::set<int>()};
			labelMap[currentLabel] = newProp; 
		}
		lutGraph[child].mffcLabel = labelStack.back();
	
	}
	else // there is a label
	{
		// check for conflict
		bool isDependent = false;
		int activeLabel = labelStack.back();
		if(childLabel != activeLabel) // conflict
		{
			if(labelMap.find(childLabel) != labelMap.end())
			{
				isDependent = labelMap[childLabel].dependsOn.find(activeLabel)!=labelMap[childLabel].dependsOn.end();
			}
			if(isDependent) // no new label is needed
			{
				lutGraph[child].mffcLabel = activeLabel;
			
			}
			else if(labelMap[childLabel].node == child)
			{
				labelMap[activeLabel].dependsOn.insert(childLabel);
				labelStack.push_back(childLabel);
				popLabel = true;
			}
			else
			{
				#ifdef ORDERING_DEBUG
			   
				std::cout<< "conflict @" << child << " " << childLabel << " x " << activeLabel <<"\n";
				#endif
			   
				// add a new label to stack
				currentLabel = currentLabel+1;
				labelStack.push_back(currentLabel);
				labelMap[activeLabel].dependsOn.insert(currentLabel);
				labelMap[childLabel].dependsOn.insert(currentLabel);
				lutGraph[child].mffcLabel = currentLabel;
				popLabel = true;
				LabelProperties newProp = {child, 0, std::vector<int>(), std::set<int>()};
				labelMap[currentLabel] = newProp; 
			}
		}
	 }

	//. proces
	std::pair<in_edge_iterator, in_edge_iterator> inEdges = boost::in_edges(child, lutGraph);
	#ifdef ORDERING_DEBUG
		
	std::cout << "In edges: " << child << std::endl;
	#endif
	for(; inEdges.first != inEdges.second; ++inEdges.first)
	{
		#ifdef ORDERING_DEBUG
		std::cout << "Processing child : " << *inEdges.first << "\n";
		#endif
		int vertex = source(*inEdges.first,lutGraph);
		visitChild(vertex);
	}
	
	// remove label from stack if any added
	if(popLabel)
		labelStack.pop_back();


}

std::set<int> LutGraph::dependencies(int root)
{
	std::set<int> dependsOn;
	std::set<int> queued;
	std::vector<int> process;

	for(auto label: labelMap[lutGraph[root].mffcLabel].dependsOn)
	{

		if(computedSet.find(labelMap[label].node) == computedSet.end())
		{
			process.push_back(label);
			dependsOn.insert(label);
			queued.insert(label);
		}

	}

	while(!process.empty())
	{
		int label = process.front();
		process.erase(process.begin());
		
		for(auto clabel : labelMap[label].dependsOn)
		{
			if(queued.find(clabel) == queued.end() && 
					computedSet.find(labelMap[clabel].node) == computedSet.end())
			{
				process.push_back(clabel);
				dependsOn.insert(clabel);
				queued.insert(clabel);
			}
		}

	}
	/*
	for(auto label: labelMap[lutGraph[root].mffcLabel].dependsOn)
	{
		if(computedSet.find(labelMap[label].node) == computedSet.end())
		{
			dependsOn.insert(label);
			std::set<int> ddependsOn = dependencies(labelMap[label].node);
			dependsOn.insert(ddependsOn.begin(), ddependsOn.end());
				
		}
	 }*/
	#ifdef ORDERING_DEBUG
	std::cout << "dependencies complete\n";
	#endif
	return dependsOn;
}

void LutGraph::computeCone(int coneRoot)
{
	 perfCounter++;
	if(perfCounter % 10 == 0)
		std::cout << "perf compute: " << perfCounter << "\n";
	
	// check if the root has already been computed
	if(computedSet.find(coneRoot) != computedSet.end())
		return;
	std::priority_queue<prioNode> q;
	std::set <int> dependOn = dependencies(coneRoot);

	// order dependencies in a priority queue
	for(int dlabel: dependOn)            
	{
		prioNode dNode = {lutGraph[labelMap[dlabel].node].topoLabel, labelMap[dlabel].node, dlabel };
		if(computedSet.find(labelMap[dlabel].node) == computedSet.end())
			q.push(dNode);
	   // std::cout << "q size:" << q.size() << "\n"; 
	}
			

	// compute dependencies
	while(!q.empty())
	{
		prioNode dNode = q.top();
		#ifdef ORDERING_DEBUG
		std::cout << "depComp" << dNode.node << " : level" << lutGraph[dNode.node].topoLabel << "\n";
		#endif
		q.pop();
		if(computedSet.find(dNode.node) == computedSet.end())
		   computeCone(dNode.node); 
	}        

	compute(coneRoot);
	#ifdef ORDERING_DEBUG

	std::cout << "computed cone :n" << coneRoot << ":l" << lutGraph[coneRoot].mffcLabel << "\n";
	#endif

}
void LutGraph::compute(int root)
{
   // order fanins 
	std::priority_queue<prioNode> q;
	for(auto const node : labelMap[lutGraph[root].mffcLabel].fanIn)
	{
		prioNode fnode = {lutGraph[node].topoLabel, node, 0};
		q.push(fnode);
	}    
	std::vector<int> computeOrder;
	while(!q.empty())
	{
		prioNode fnode = q.top();
		q.pop();
		if(computedSet.find(fnode.node) != computedSet.end())
		{
			// an already computed node
			continue;
		}
		#ifdef ORDERING_DEBUG

		std::cout << "lut" << fnode.node << " computed\n";
		#endif

		// add the step here globally
		steps.push_back(std::make_pair(0, fnode.node));
	
		// add to the computed set
		computedSet.insert(fnode.node);

		// add step to the global compute order root -> compute
		computeOrder.push_back(fnode.node);
	   
		// add step to the mffcLabel steps
		labelCompute[lutGraph[fnode.node].mffcLabel].push_back(fnode.node);
		// add the successors of current node to q with same label!
		std::pair<out_edge_iterator, out_edge_iterator> outEdges = boost::out_edges(fnode.node, lutGraph);
				 
		for(; outEdges.first != outEdges.second; ++outEdges.first)
		{
				int vertex = target(*outEdges.first,lutGraph);
				if(lutGraph[vertex].mffcLabel != lutGraph[fnode.node].mffcLabel)
					continue;

				prioNode dnode = {lutGraph[vertex].topoLabel,
						vertex, lutGraph[vertex].mffcLabel};
				q.push(dnode);
		}
			
	}
	// add current label to labelOrdering 
	labelOrdering.push_back(lutGraph[root].mffcLabel);
	if(std::find(po.begin(), po.end(),root) != po.end())
	{
		// This is a primary output
		computedOutput.insert(lutGraph[root].mffcLabel);
	}
}

bool LutGraph::uncomputeCone(int coneRoot, bool output, bool force= false)
{
	#ifdef ORDERING_DEBUG

	std::cout << "Attempting uncompute n" << coneRoot << ":l" << lutGraph[coneRoot].mffcLabel << "\n";
	#endif

	if(std::find(labelOrdering.begin(),labelOrdering.end(),lutGraph[coneRoot].mffcLabel) == labelOrdering.end())
	{
		#ifdef ORDERING_DEBUG  
		std::cout << coneRoot << " not computed yet!";
		#endif
		return false;
	}
	//  check if successors have already been uncomputed
	//  else abort
	bool predUncomputed = true;
	std::pair<out_edge_iterator, out_edge_iterator> outEdges = boost::out_edges(coneRoot, lutGraph);
		 
	for(; outEdges.first != outEdges.second; ++outEdges.first)
	{
		
		int vertex = target(*outEdges.first,lutGraph);
		// check if successor is present in computed set and is 
		// not an output
		if(computedSet.find(vertex) != computedSet.end() && 
		  std::find(po.begin(), po.end(), vertex) == po.end())
	   {

		   #ifdef ORDERING_DEBUG
		   std::cout << "n" << coneRoot << " successor n" << vertex << "has not been uncomputed\n";
		   #endif
		   predUncomputed = false;
		   break;
	   } 
	
	}
	if(!predUncomputed)
		return false;

	if(! force)
	{ // determine if cone needed in computation of other outputs
		 bool supportCone = false; 
		for(int supportedCone : labelMap[lutGraph[coneRoot].mffcLabel].supports)
		{
			// one of the output cone the current cone supports 
			// has not been computed
			if(computedOutput.find(supportedCone) == computedOutput.end())
			{
				#ifdef ORDERING_DEBUG

				std::cout << "Supporting l" << supportedCone << "<- n" << coneRoot <<"\n";
				#endif

				supportCone = true;
				break;
			}
		}
		if(supportCone)
			return false;
		
	}
	uncompute(coneRoot, output);
	#ifdef ORDERING_DEBUG
	std::cout << "uncomputed cone :n" << coneRoot << ":l" << lutGraph[coneRoot].mffcLabel << "\n";
	#endif
	// order the dependencies wrt to inverse of labelOrdering
	std::priority_queue<prioNode> q;

	// order dependencies in a priority queue
	for(auto const dlabel: labelMap[lutGraph[coneRoot].mffcLabel].dependsOn)            
	{
		std::vector<int>::iterator pos;
		pos = std::find(labelOrdering.begin(),labelOrdering.end(), dlabel);
		if(pos == labelOrdering.end())
		{
			#ifdef ORDERING_DEBUG

			std::cout << " Already uncomputed l" << dlabel << "\n";
			#endif

			continue;

		}
		int index = pos - labelOrdering.begin();
		prioNode dNode = {index, labelMap[dlabel].node, dlabel};
	   
		q.push(dNode);
	}
			
	int curr = -1;
	// uncompute dependency cones
	while(!q.empty())
	{
		#ifdef ORDERING_DEBUG
		std::priority_queue<prioNode> q1;
		q1 = q;
		std::cout << "q:";
		while(!q1.empty())
		{
			std::cout << q1.top().node << " ";
			q1.pop();
		}
		#endif
		prioNode dNode = q.top();
		q.pop();
		if(dNode.node == curr) // retrying same uncompute again -- skip
			continue;
		curr = dNode.node;
		// the node has already been uncomputed!
		if(computedSet.find(dNode.node) == computedSet.end() ||
				std::find(labelOrdering.begin(), labelOrdering.end(), lutGraph[dNode.node].mffcLabel) == labelOrdering.end())
			continue;
		// not output, do not force uncompute
		bool unc_succ = uncomputeCone(dNode.node, false, false); 
		if(!unc_succ)
			continue;
		
		 // order dependencies in a priority queue
		for(auto const dlabel: labelMap[lutGraph[dNode.node].mffcLabel].dependsOn)            
		{
			std::vector<int>::iterator pos;
			pos  = std::find(labelOrdering.begin(), labelOrdering.end(), dlabel);
			if(pos == labelOrdering.end())
			{ 
				#ifdef ORDERING_DEBUG

				std::cout << " Already uncomputed l" << dlabel << "\n";
				#endif

				continue; 
			}
			int index = pos - labelOrdering.begin();
			prioNode dNode = {index*-1, labelMap[dNode.node].node, dlabel};
		   
			q.push(dNode);
		}  

		#ifdef ORDERING_DEBUG
		 q1 = q;
		std::cout << "\nq1:";
		while(!q1.empty())
		{
			std::cout << q1.top().node << " ";
			q1.pop();
		}
		#endif
	}        
	return true;
}


void LutGraph::uncompute(int root, bool output)
{
	int mffcLabel = lutGraph[root].mffcLabel;
	std::vector<int> computeOrder = labelCompute[mffcLabel];
	if(output)
		computeOrder.pop_back(); // no need to uncompute output
	while(!computeOrder.empty())
	{

		int node =  computeOrder.back();
		computeOrder.pop_back();
		steps.push_back(std::make_pair(1, node));
		#ifdef ORDERING_DEBUG
		std::cout << "lut" << node << "uncomputed for root" << root <<"\n";
		#endif
		// remove node from computedSet
		computedSet.erase(node);
	} 

	// remove label from labelOrdering vector
	std::vector<int>::iterator it;
	it = std::find(labelOrdering.begin(), labelOrdering.end(), mffcLabel);
	if(it != labelOrdering.end())
	{
		labelOrdering.erase(it);
	}
}

void LutGraph::determineMFFC()
{
	labelAndInit();
}
void LutGraph::countPending(int label, pendingCone& p)
{
	if(computedSet.find(labelMap[label].node) != computedSet.end())
		return ;
	else
	{
		p.nodes = p.nodes + labelMap[label].numberOfNodes;
		for(int dlabel : labelMap[label].dependsOn)
		{
			if(p.dependsOn.find(dlabel) == p.dependsOn.end())
			{
				p.dependsOn.insert(dlabel);
				countPending(dlabel,p);
			}
		}
	}
}

void LutGraph::evalAddNodes(int nodes, int &allocatedNodes, int &freeNodes)
{
	if(freeNodes > nodes)
	{
		freeNodes = freeNodes - nodes;    
	}
	else
		freeNodes = 0;
	allocatedNodes = allocatedNodes + nodes; 

}
void LutGraph::evalRemoveNodes(int nodes, int &allocatedNodes, int &freeNodes, bool output)
{
	int adjust = output? 1 : 0;
	freeNodes = freeNodes + nodes - adjust;    
	assert(allocatedNodes >= nodes && " Invalid removal of nodes ");     
	allocatedNodes = allocatedNodes - nodes + adjust; 
}

int LutGraph::evaluateCost(std::vector<int> poOrder)
{
	std::set<int> computed;
	int allocatedNodes, freeNodes;
	allocatedNodes = freeNodes = 0;
	std::queue<int> pendingQ; 
	std::set<int> addedSet;

	int poLabel, succLabel;
	for(int po : poOrder)
	{
		std::cout << po << " outputs" ;
		std::cout << "\n";
		poLabel = lutGraph[po].mffcLabel;
		evalAddNodes(labelMap[poLabel].numberOfNodes, allocatedNodes, freeNodes);
		computed.insert(poLabel);
		// add successor
		for(int succ : labelMap[poLabel].dependsOn)
		{
			if(computed.find(succ) == computed.end())
			{
					addedSet.insert(succ);
					pendingQ.push(succ);
			}        
		}

		while(!pendingQ.empty())
		{
			succLabel = pendingQ.front();
			pendingQ.pop();
			if(computed.find(succLabel) == computed.end())
			{
				evalAddNodes(labelMap[succLabel].numberOfNodes, allocatedNodes, freeNodes);
				computed.insert(succLabel);
			}
			for(int succ : labelMap[succLabel].dependsOn)
			{
				if(computed.find(succ) == computed.end() && addedSet.find(succ) == addedSet.end())
				{
					addedSet.insert(succ);
					pendingQ.push(succ);
				}
			}// end of for

		}//end of while
		//*
		#ifdef COST_DEBUG
		std::cout << "computed:";
		for(int node:computed)
			std::cout<< node << " " ;
		std::cout << "\n";
		#endif
		//*/

		// uncompute   
		evalRemoveNodes(labelMap[poLabel].numberOfNodes, allocatedNodes, freeNodes, true); //output uncompute
		pendingQ = std::queue<int>();
		addedSet = std::set<int>(); 
		// add successor
		for(int succ : labelMap[poLabel].dependsOn)
		{
			if(computed.find(succ) != computed.end())
			{
					addedSet.insert(succ);
					pendingQ.push(succ);
			}        
		}

		while(!pendingQ.empty())
		{
			//std::cout<< "attempt : " << succLabel << " "; 
			succLabel = pendingQ.front();
			pendingQ.pop();
			if(computed.find(succLabel) != computed.end())
			{
				// check if supports any uncomputed output
				bool support = false;
				for(int suppPo : labelMap[succLabel].supports)
					if(computed.find(suppPo) == computed.end())
					{
						//std::cout << suppPo << " -< \n";
						support = true;
						break;
					}
				if(!support)
				{
					evalRemoveNodes(labelMap[succLabel].numberOfNodes, allocatedNodes, freeNodes, false);
					//std::cout<<"u" << succLabel << " ";
					computed.erase(succLabel);
				}
				for(int succ : labelMap[succLabel].dependsOn)
				{
					if(computed.find(succ) != computed.end() && addedSet.find(succ) == addedSet.end())
					{
						addedSet.insert(succ);
						pendingQ.push(succ);
					}
				}// end of for

			}
	  
		}//end of while
		/*std::cout << "\nuncomputed:";
		for(node:computed)
			std::cout<< node << " " ;
		std::cout << "\n";
		*/


	} // end of for
	std::cout << "Allocated: " << allocatedNodes << " Free:" << freeNodes << " Total:" 
			<< allocatedNodes+freeNodes << "\n";
	// for(std::map<int,struct LabelProperties>::iterator it = labelMap.begin(); it != labelMap.end(); ++it) {
	//     std::cout << it->first << "\n";
	//}
	return allocatedNodes+freeNodes;
}


void LutGraph::genSolution(std::vector<int> solution, std::vector<int> &sol, int k, double alpha)
{
	int oldCost = evaluateCost(solution); 
	std::vector<int> bestSolution = solution;
	int bestCost = oldCost, newCost;

	double T = 1.0; 
	double T_min = 0.00001;
	int i;
	// set up random number related variables
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 eng(rd()); // seed the generator
	std::uniform_int_distribution<> distr(0, solution.size()-1); // define the range
	
	std::mt19937_64 rng;
	std::uniform_real_distribution<double> unif(0, 1);
	while(T > T_min && i < k)
	{	
		// generate a new permutation with two positions swapped
		int p1, p2, count = 0;
		do
		{
			p1 = distr(eng);
			p2 = distr(eng);	    
		}while(p1 == p2 && count == 100);
		std::vector<int> newSolution;
		newSolution = solution;
	   int move = newSolution[p1];
		newSolution[p1] = newSolution[p2];
		newSolution[p2] = move; 
		std::cout << p1 << "," << p2 << " : ";
		for(int out: newSolution)
			std::cout << out << " ";
		std::cout << "\n";
	 
		newCost = evaluateCost(newSolution);
		double ap = exp((newCost - oldCost)/T); // acceptance_probability
		if(ap > unif(rng))
		{
			solution = newSolution;
			oldCost = newCost;
			std::cout << "iteration "<< i << " : " << newCost << "\n";
		}
		if(bestCost > newCost)
		{
			bestSolution = newSolution;
			bestCost = newCost;
		}
		i++;
		T = T*alpha;
	} // End of while
	sol = bestSolution;
}

unsigned int LutGraph::scheduleGraph(int orderType, bool crit = false)
{

	// determine the minimum or maximum
	std::set<int> pendingOutputSet;
	std::vector<int> poComputeOrder; 
	for(auto const output : po)
	{
		pendingOutputSet.insert(output);
	}
	bool order = orderType == 0;// 0 = ascending 1 = descending
	int compVal = order ? 10000000 : -1;   
	while(!pendingOutputSet.empty())
	{

		pendingCone minCone = {compVal, std::set<int>()};
		int minOut = -1;
		for(auto const output : pendingOutputSet)
		{
			pendingCone cone = {0, std::set<int>()};
			countPending(lutGraph[output].mffcLabel, cone);
			if((order && cone.nodes < minCone.nodes)
				|| (!order && cone.nodes > minCone.nodes))
			{
				minCone = cone;
				minOut = output;

			}

		} 
		
		pendingOutputSet.erase(minOut);

	//for(auto const minout: po)
	//{   
		// compute the output
		std::cout << "Computing PO n" << minOut << ": l"<< lutGraph[minOut].mffcLabel <<" #nodes = " <<  minCone.nodes << "\n";
		perfCounter = 0;
		computeCone(minOut);
		poComputeOrder.push_back(minOut);
		#ifdef ORDERING_DEBUG

		std::cout << " Computed :";
		for(auto out : computedOutput)
			std::cout <<"l" << out << " ";
		std::cout << "\n";
		#endif

		// Uncompute the output cone
		uncomputeCone(minOut,true);


	}


	#ifdef ORDERING_DEBUG
	for(auto comp:steps) 
	{
		std::cout << comp.first << ":" << comp.second << " ";
	}
	std::cout << "\n";
	std::cout << "\n Final computed set:";
	for(auto el : computedSet)
		std::cout << el << " ";
	std::cout << "\n";
	#endif
	unsigned int reqQbits = evaluateCost(poComputeOrder);
	//std::vector<int> optiPerm;
	//genSolution(poComputeOrder, optiPerm, 100, 0.99);
	return reqQbits;
}
//} // end of namespace

/**************************************************************************************************



       					EVERYTHING BELOW THIS AND topoOrdering() ABOVE




***************************************************************************************************/
void LutGraph::preprocess_pebble()
{
	// finding leaf and target nodes
	Graph::vertex_iterator vIt, vItEnd;
	boost::tie(vIt,vItEnd) = vertices(lutGraph);
	int idCount = 0;
	for(; vIt != vItEnd; vIt++)
	{
		lutGraph[*vIt].id = idCount++; // setting ID for vertices
		lutGraph[*vIt].topoLabel = -1;
		lutGraph[*vIt].qubit_alloc = -1;

		if(boost::in_degree(*vIt,lutGraph) == 0)
		{    
			pi.push_back(*vIt);
			lutGraph[*vIt].topoLabel = 0;
		}
	}
	#ifdef ORDERING_DEBUG
	for(int h=0;h<po.size();h++)
		std::cout<<"OUT : "<<po[h]<<"\n";

	for(int h=0;h<pi.size();h++)
		std::cout<<"IN : "<<pi[h]<<"\n";
	#endif
	/*******************/



	// order vertices topologically --- changed topoOrdering above.
	for(auto input: pi)      
		topoOrdering(input);

	#ifdef ORDERING_DEBUG
	boost::tie(vIt,vItEnd) = vertices(lutGraph);
	for(; vIt != vItEnd; vIt++)
		std::cout<<"vertex  :"<<*vIt<< "    ID :  " <<lutGraph[*vIt].id <<"   TOPO:  "<<lutGraph[*vIt].topoLabel<<"\n";
	#endif
	/*******************/



	// set up vectors for successor/predecessor for each vertices  
	boost::tie(vIt,vItEnd) = vertices(lutGraph);
	for(; vIt != vItEnd; vIt++)
	{
		std::pair<out_edge_iterator, out_edge_iterator> outEdges = out_edges(*vIt, lutGraph);        
		for(out_edge_iterator iter = outEdges.first; iter != outEdges.second; ++iter)
		{
			
			int t = target(*iter, lutGraph);
			lutGraph[*vIt].successor_vertices.push_back(t);
		}

		std::pair<in_edge_iterator, in_edge_iterator> inEdges = in_edges(*vIt, lutGraph);        
		for(in_edge_iterator iter = inEdges.first; iter != inEdges.second; ++iter)
		{
			
			int s = source(*iter, lutGraph);
			lutGraph[*vIt].predecessor_vertices.push_back(s);
		}
	}

	#ifdef ORDERING_DEBUG
	std::cout<<"successor \n";
	boost::tie(vIt,vItEnd) = vertices(lutGraph);
	for(; vIt != vItEnd; vIt++)
	{
		std::cout<<lutGraph[*vIt].id<<": ";
		for(int h=0;h<lutGraph[*vIt].successor_vertices.size();h++)
		{
			std::cout<<lutGraph[*vIt].successor_vertices[h]<<" ";
		}
		std::cout<<"\n";
	}
	
	std::cout<<"predecessor \n";
	boost::tie(vIt,vItEnd) = vertices(lutGraph);
	for(; vIt != vItEnd; vIt++)
	{
		std::cout<<lutGraph[*vIt].id<<": ";
		for(int h=0;h<lutGraph[*vIt].predecessor_vertices.size();h++)
		{
			std::cout<<lutGraph[*vIt].predecessor_vertices[h]<<" ";
		}
		std::cout<<"\n";
	}
	#endif
	/********************/



	//bfs from each terminal node to find labels of node supporting output nodes
	//priority queue for level order sort (node id, topolevel)
    std::priority_queue<std::tuple<int,int>, std::vector<std::tuple<int,int>>,ComparePrioTopo> bfs_prio;
    std::set<int> added; 
    int dependents, v;
#ifdef ORDERING_DEBUG 
        std::cout << "Number of vertices : " << num_vertices(lutGraph) << "\n";
#endif
	for(auto output:po)
	{
        for(int i=0;i<4;i++)
            output_stats[output].push_back(0);
		//output_stats[output].resize(4);
        //		bfs_queue.push(output);
        added = std::set<int>();
        bfs_prio.push(std::make_tuple(output, lutGraph[output].topoLabel));
				
        dependents = 0; 

        int n = num_vertices(lutGraph) ;
        int count = 0;
        while(!bfs_prio.empty())
        {
            count ++ ;
            v = std::get<0>(bfs_prio.top());
            if(count % 1000 == 0)
            std::cout << "v" << v << " " << n << " " << count << " \n" ;
            bfs_prio.pop();

            if(added.find(v) != added.end())
                continue; 
            added.insert(v); 
			lutGraph[v].output_labels.insert(output);
            dependents++;
            std::pair<in_edge_iterator, in_edge_iterator> inEdges = in_edges(v, lutGraph);        
			for(in_edge_iterator iter = inEdges.first; iter != inEdges.second; ++iter)
			{                
				int s = source(*iter, lutGraph);
				bfs_prio.push(std::make_tuple(s, lutGraph[s].topoLabel));
			}
        }
       	output_stats[output][0] = dependents; 
		output_stats[output][1] = lutGraph[output].topoLabel;

		lutGraph[output].output_labels.insert(output);
        /*	
		//std::set<int> no_of_dependent;
		lutGraph[output].output_labels.insert(output);
		while(bfs_queue.size()!=0)
		{
			int v = bfs_queue.front();
			//std::cout<<v<<"\n";
			bfs_queue.pop();
			std::pair<in_edge_iterator, in_edge_iterator> inEdges = in_edges(v, lutGraph);        
			for(in_edge_iterator iter = inEdges.first; iter != inEdges.second; ++iter)
			{                
				int s = source(*iter, lutGraph);
				lutGraph[s].output_labels.insert(output);
				bfs_queue.push(s);
				no_of_dependent.insert(s);
			}
		}
        std::cout <<output << ":"<< no_of_dependent.size() << " " << dependents << " << Dep\n"; 
		output_stats[output][0] = no_of_dependent.size();
:w
output_stats[output][1] = lutGraph[output].topoLabel;
        */
	}

	#ifdef ORDERING_DEBUG
	std::cout<<" Label test \n";
	boost::tie(vIt,vItEnd) = vertices(lutGraph);
	for(; vIt != vItEnd; vIt++)
	{
		std::cout<<lutGraph[*vIt].id<<": ";
		std::set<int>::iterator it;
		for (it=lutGraph[*vIt].output_labels.begin(); it!=lutGraph[*vIt].output_labels.end(); ++it)
			std::cout << ' ' << *it;
		std::cout<<"\n";
	}
	#endif
	/************************/
}

float LutGraph::contrib(int n,int output_label)
{
	float contrib = 0;        
	for(uint i=0;i<lutGraph[n].successor_vertices.size();i++)
	{   
		int succ = lutGraph[n].successor_vertices[i];
		std::set<int> vec = lutGraph[succ].output_labels;
		int count = 0;
		if (std::find(vec.begin(), vec.end(), output_label) != vec.end())
		{
			for(uint j=0;j<lutGraph[succ].predecessor_vertices.size();j++)
			{   
				int pred = lutGraph[succ].predecessor_vertices[j];     
				if(lutGraph[pred].qubit_alloc == -1)        
					count++;
			}
			contrib += (1.0/count);
		}
	}
	return contrib;
}

bool LutGraph::eligible(int node, int output_label)
{
	bool ret = false;    
	for(uint j=0;j<lutGraph[node].predecessor_vertices.size();j++)
	{   
		int pred = lutGraph[node].predecessor_vertices[j];     
		if(lutGraph[pred].qubit_alloc != -1)        
			ret = true;
		else
		{
			ret = false;
			return ret;
		}
	}
	if (std::find(pi.begin(), pi.end(), node) != pi.end()) // it is a leaf node
		return true;
	return ret;
}
bool LutGraph::succDep(int node, int output_label) // doing as described in text
{
	bool ret = false;    
	for(uint j=0;j<lutGraph[node].successor_vertices.size();j++)
	{   
		int succ = lutGraph[node].successor_vertices[j];
		std::set<int> vec = lutGraph[succ].output_labels;
		if (std::find(vec.begin(), vec.end(), output_label) != vec.end())     
		{
			#ifdef ORDERING_DEBUG
			if (node == 250)
			{
				std::cout<<"DEBUG4\n";
				std::cout<<succ<<" "<<( computedOnce.find(succ) != computedOnce.end() )<<" "<<output_label<<"\n";
			}
			#endif

			if(computedOnce.find(succ) != computedOnce.end())        
				ret = true;
			else
			{
				ret = false;
				return ret;
			}	
		}
	}
	return ret;
}
bool LutGraph::succUncomp(int node, int output_label) // doing as described in text
{
	bool ret = false;    
	for(uint j=0;j<lutGraph[node].successor_vertices.size();j++)
	{   
		int succ = lutGraph[node].successor_vertices[j];
		std::set<int> vec = lutGraph[succ].output_labels;
		if (std::find(vec.begin(), vec.end(), output_label) != vec.end())    
		{
			if(lutGraph[succ].qubit_alloc == -1)        
				ret = true;
			else
			{
				ret = false;
				return ret;
			}
		}
	}
	return ret;
}


int LutGraph::get_compute_qubit(int top)
{
	if(qubits_list.size() > 0) // if free qubits already there
	{
		int r = qubits_list.top();
		qubits_list.pop();
		if(qubits_bool[r] == true)
		{
			qubits_used_yet++;
			qubits_bool[r] = false;
		}
		current_qubits ++;
		return r;
	}
    std::set<int> added_node; 
	//do topo sort to find the node to uncompute according to algo
	topo_prio.push(std::make_tuple(top,lutGraph[top].topoLabel));
	int found = -1;
	int uncomp = -1;//final node selected to uncompute

	while(topo_prio.size()>0)
	{
		int n = std::get<0>(topo_prio.top());
		//std::cout<<n <<" "<<eligible(n)<<" "<<succDep(n)<<" "<<succUncomp(n)<<"\n";
		topo_prio.pop();
        added_node.insert(n);
		if(found == 1 && lutGraph[n].topoLabel < lutGraph[uncomp].topoLabel)//to stop going to nxt level once found
			break;

		//std::cout<<"IN "<<n<<"\n";
		if(lutGraph[n].qubit_alloc !=-1 && eligible(n,top))
		{
			#ifdef ORDERING_DEBUG
			if (n == 250)
			{
				std::cout<<" SUCCDEP "<<succDep(n,top)<<"\n";

			}
			#endif

			if(succDep(n,top))
			{
				found = 1;uncomp = n;
				if(succUncomp(n,top))
					break;
			}
		}
		for(uint i=0;i<lutGraph[n].predecessor_vertices.size();i++)
        {
            if(added_node.find(lutGraph[n].predecessor_vertices[i]) == added_node.end())
            {
				topo_prio.push(std::make_tuple(lutGraph[n].predecessor_vertices[i],lutGraph[lutGraph[n].predecessor_vertices[i]].topoLabel));
                added_node.insert(lutGraph[n].predecessor_vertices[i]);
            }
        
        }
	}
	if(found == 1)
	{
		//uncompute
		int q = lutGraph[uncomp].qubit_alloc;
		lutGraph[uncomp].qubit_alloc = -1;
		computedSet.erase(uncomp);

		pebble_steps.push_back(std::make_tuple(1,uncomp,q,top));
        #ifdef ORDERING_DEBUG 
		std::cout<<"uncompute  node: "<<uncomp<<"  qubit freed: "<<q<<"\n";
        #endif
		current_steps++;
		return q;
	}
	else
	{
		qubits_available++;
		current_qubits++;
		qubits_used_yet++;
		qubits_bool.push_back(false);
		std::cout<<"ADDED EXTRA QUBIT, TOTAL QUBITS NOW : "<<qubits_available<<"\n";
		return qubits_available;
	}

}

void LutGraph::reorder_prio(int output_label)
{
	std::vector<int> v;
	while(pebble_prio.size()>0)
	{
		v.push_back(std::get<0>(pebble_prio.top()));
		pebble_prio.pop();
	}
	for(uint i =0;i<v.size();i++)
	{
		pebble_prio.push(std::make_tuple(v[i], lutGraph[v[i]].topoLabel, contrib(v[i],output_label)));
	}
}

void LutGraph::pebble_compute(int output_label)
{
	//create the  priority Q for the 
	computedOnce = std::set<int>();
	for(auto input:pi)
	{
		std::set<int> vec = lutGraph[input].output_labels;
		if (std::find(vec.begin(), vec.end(), output_label) != vec.end())
		{
 		pebble_prio.push(std::make_tuple(input, lutGraph[input].topoLabel, contrib(input,output_label)));
		}	
	}
	/*if(output_label == 9)
	{
		while (!pebble_prio.empty())
		{
			std::cout << ' ' << std::get<0>(pebble_prio.top()) <<" "<<std::get<2>(pebble_prio.top());
			pebble_prio.pop();
		}
		return;
	}*/

	#ifdef ORDERING_DEBUG
	int counterp = 0;

	std::queue<int> bfs_temp;
	std::set<int> dep_set;
	bfs_temp.push(output_label);
	while(!bfs_temp.empty())
	{
		int n_temp = bfs_temp.front();
		dep_set.insert(n_temp);
		bfs_temp.pop();
		std::cout<<n_temp<<" < ";
		for(uint j=0;j<lutGraph[n_temp].predecessor_vertices.size();j++)
		{   
			int pred = lutGraph[n_temp].predecessor_vertices[j];
			std::cout<<pred<<" ";
			// if (dep_set.find(pred) != dep_set.end())
			bfs_temp.push(pred);
		}
		std::cout<<"\n";
	}
	#endif
	

	while(computedSet.count(output_label) != 1)
	{
		int node;

		#ifdef ORDERING_DEBUG
		std::cout<<"Curr Allooc \n";
		for (auto nn : dep_set)
		{		
			std::cout<<"["<<nn<<" : "<<lutGraph[nn].qubit_alloc<< " : "<< eligible(nn,output_label)<<" ]"<<" ";
		}
		std::cout<<"\n";
		#endif

		if (!pebble_prio.empty())
		{
			node = std::get<0>(pebble_prio.top());
			pebble_prio.pop();
		}
		else
		{
			std::cout<<"no of nodes = "<<output_stats[output_label][0]<<"\n";
			std::cout<<"empty!!!\n";
			exit(0);
		}
		
		#ifdef ORDERING_DEBUG
		counterp ++;
		std::cout<<"no of nodes = "<<output_stats[output_label][0]<<"\n";
		std::cout<<"ccpos "<<counterp<<"\n";
		#endif
		// while(topo_prio.size() > 0)
		// 	topo_prio.pop();
		topo_prio = std::priority_queue<std::tuple<int,int>,std::vector<std::tuple<int,int>>,ComparePrioTopo> ();
		int q = get_compute_qubit(output_label); // get the qubit tht will be used to compute

		//compute
		computedSet.insert(node);
		computedOnce.insert(node);
		lutGraph[node].qubit_alloc = q;

		pebble_steps.push_back( std::make_tuple(0,node,q,output_label));
        #ifdef ORDERING_DEBUG 
		std::cout<<"compute  node: "<<node<<"  qubit: "<<q<<"\n";
        #endif
		current_steps ++;

		// reorder_prio(output_label);//need to reorder the priority Q after one node is computed


		//push successsor only if all predecessors computed
		for(uint i=0;i<lutGraph[node].successor_vertices.size();i++)
		{	
			int succ = lutGraph[node].successor_vertices[i];
			std::set<int> vec = lutGraph[succ].output_labels;
			if (std::find(vec.begin(), vec.end(), output_label) != vec.end())// checking if curretn ouput supported by this node
			{
				bool toAdd = false;    
				for(uint j=0;j<lutGraph[succ].predecessor_vertices.size();j++)
				{   
					int pred = lutGraph[succ].predecessor_vertices[j];     
					if(lutGraph[pred].qubit_alloc != -1) //check if all pred of succ are computed
						toAdd = true;
					else
					{
						toAdd = false;
						break;
					}
				}
				if(toAdd)
				{
					pebble_prio.push(std::make_tuple(succ, lutGraph[succ].topoLabel, contrib(succ,output_label)));
				}
			}
		}
		reorder_prio(output_label);//need to reorder the priority Q after one node is computed
	}
	computedOutput.insert(output_label);
}

void LutGraph::find_dep_set(int node)
{
	for(uint i = 0; i< lutGraph[node].predecessor_vertices.size() ; i++) //my computed preds
		if(lutGraph[lutGraph[node].predecessor_vertices[i]].qubit_alloc != -1)
			dependency_set.insert(lutGraph[node].predecessor_vertices[i]);
		else
			find_dep_set(lutGraph[node].predecessor_vertices[i]); //calling again on uncomputed preds

}

void LutGraph::compute_during_uncomp(int node, int parent, int output_label)//parent is the original node whch I wanted to uncompute. root
{
	if(!eligible(node,output_label)) // one or preds are uncomputed. compute them first
	{
		std::vector<int> pred_not_computed;
		for(uint j=0;j<lutGraph[node].predecessor_vertices.size();j++)
		{
			if(lutGraph[lutGraph[node].predecessor_vertices[j]].qubit_alloc == -1)
				pred_not_computed.push_back(lutGraph[node].predecessor_vertices[j]);
		}
		for(uint i =0; i < pred_not_computed.size(); i++)
		{
			int node = pred_not_computed[i];
			compute_during_uncomp(node,parent,output_label);
		}
	} 

	if(eligible(node,output_label))// all preds computed ?
	{
		//find a qubit and compute

		// if qubits are free in list 
		if(qubits_list.size()>0)
		{
			//compute
			int q = qubits_list.top();
			qubits_list.pop();
			computedSet.insert(node);
			lutGraph[node].qubit_alloc = q;

			pebble_steps.push_back(std::make_tuple(0,node,q,output_label));
            #ifdef ORDERING_DEBUG 
			std::cout<<"compute  node: "<<node<<"  qubit: "<<q<<"\n";
            #endif
			current_steps++;
		}
		else
		{
			// need to find a suitable node to uncompute and get free qubit for computing this one

			//build dependency_set
			dependency_set.clear();
			// do not pick a qubit from a) computed children of original parent b) my computed children recursively
			for(uint i = 0; i< lutGraph[parent].predecessor_vertices.size() ; i++)
				if(lutGraph[lutGraph[parent].predecessor_vertices[i]].qubit_alloc != -1)
					dependency_set.insert(lutGraph[parent].predecessor_vertices[i]);

			find_dep_set(node);
			
			/*std::cout<<"For node "<<node<<"\n";
			for (std::set<int>::iterator it=dependency_set.begin(); it!=dependency_set.end(); ++it)
			    std::cout << ' ' << *it;
			std::cout << '\n';*/

			bool found = false;
			int uncomp = -1;// node to uncompute

			//build prio Q for topo sort from root to find free qbits
			while(topo_prio.size() > 0)
				topo_prio.pop();// clearing any junk values

			topo_prio.push(std::make_tuple(output_label,lutGraph[output_label].topoLabel));//push output(root)

			while(topo_prio.size()>0)
			{
				int curr = std::get<0>(topo_prio.top());
				topo_prio.pop();

				if(curr != output_label && lutGraph[curr].qubit_alloc != -1 && dependency_set.count(curr)!=1 && eligible(curr,output_label))
				{
					found = 1;uncomp = curr;
					break;
				}
				
				for(uint i=0;i<lutGraph[curr].predecessor_vertices.size();i++)
					if(lutGraph[lutGraph[curr].predecessor_vertices[i]].qubit_alloc != -1)
						topo_prio.push(std::make_tuple(lutGraph[curr].predecessor_vertices[i],lutGraph[lutGraph[curr].predecessor_vertices[i]].topoLabel));
			}

			if(found)
			{
				//uncompute found qubit and compute node
				int q = lutGraph[uncomp].qubit_alloc;
				lutGraph[uncomp].qubit_alloc = -1;
				computedSet.erase(uncomp);

				pebble_steps.push_back(std::make_tuple(1,uncomp,q,output_label));
  #ifdef ORDERING_DEBUG 
				std::cout<<"uncompute  node: "<<uncomp<<"  qubit freed: "<<q<<"\n";
    #endif
				current_steps++;

				////
				computedSet.insert(node);
				lutGraph[node].qubit_alloc = q;

				pebble_steps.push_back( std::make_tuple(0,node,q,output_label));
#ifdef ORDERING_DEBUG 
                std::cout<<"compute  node: "<<node<<"  qubit: "<<q<<"\n";
#endif
				current_steps++;
			}
			else
			{
				//exrta qubit reqd
				qubits_available++;
				qubits_bool.push_back(false);
				current_qubits++;
				qubits_used_yet++;
				std::cout<<"ADDED EXTRA QUBIT, TOTAL QUBITSc NOW : "<<qubits_available<<"\n";
				int q = qubits_available;

				computedSet.insert(node);
				lutGraph[node].qubit_alloc = q;

				pebble_steps.push_back( std::make_tuple(0,node,q,output_label));
#ifdef ORDERING_DEBUG
				std::cout<<"compute  node: "<<node<<"  qubit: "<<q<<"\n";
    #endif
				current_steps++;
			}
		}
	}

}

void LutGraph::pebble_uncompute(int output_label)
{
	//creating FCO
    std::vector<int> first_compute_order; // initialise to blank everytime for every output label
	for(int i=pebble_steps.size()-1;i>=0;i--)
	{
		if(std::get<0>(pebble_steps[i]) == 0 && std::get<3>(pebble_steps[i]) == output_label)
			first_compute_order.push_back(std::get<1>(pebble_steps[i]));
	}

	for(uint i = 1;i<first_compute_order.size();i++) // leave 1st node as it is target - index from 1
	{
		int n = first_compute_order[i];
		if(lutGraph[n].qubit_alloc == -1)// Current node already uncomputed
			continue;

		//if any successor in computed list, do not uncompute
		bool check = true;
		
		/* not checcking if one of successors still computed. Case should not arise normally by algo.
		// but checking sometimes causes problems if node has a successor which was final o/p of other tree(has diff o/p label) 

		for(uint j=0;j<lutGraph[n].successor_vertices.size();j++)
		{
			if(lutGraph[n].successor_vertices[j] != output_label && lutGraph[lutGraph[n].successor_vertices[j]].qubit_alloc != -1)
			{
				check = false;
				std::cout<<"ERROR!!! One of successors still computed for "<<n<<"\n";//can this case arise?
				break;
			}

		}
		*/
		

		if (check and eligible(n,output_label)) //all successros uncomputed and all pred computed
		{
			int q = lutGraph[n].qubit_alloc;
			lutGraph[n].qubit_alloc = -1;
			computedSet.erase(n);

			pebble_steps.push_back( std::make_tuple(1,n,q,output_label));
#ifdef ORDERING_DEBUG
			std::cout<<"uncompute  node: "<<n<<"  qubit freed: "<<q<<"\n";
    #endif
			current_steps ++;
			qubits_list.push(q);
		}
		
		else // some pred not computed
		{
			//find pred not computed if any
			std::vector<int> pred_not_computed;
			for(uint j=0;j<lutGraph[n].predecessor_vertices.size();j++)
			{
				if(lutGraph[lutGraph[n].predecessor_vertices[j]].qubit_alloc == -1)
					pred_not_computed.push_back(lutGraph[n].predecessor_vertices[j]);
			}

			// compute pred one by one
			for(uint i =0; i < pred_not_computed.size(); i++)
			{
				int node = pred_not_computed[i];
				compute_during_uncomp(node,n,output_label);
			}

			// now uncompute the node we were about to
			int q = lutGraph[n].qubit_alloc;
			lutGraph[n].qubit_alloc = -1;
			computedSet.erase(n);

			pebble_steps.push_back( std::make_tuple(1,n,q,output_label));
#ifdef ORDERING_DEBUG
			std::cout<<"uncompute  node: "<<n<<"  qubit freed: "<<q<<"\n";
    #endif
			current_steps ++;
			qubits_list.push(q);
		}
		
	}
}
void LutGraph::pebble()
{
	std::cout<<"STARTING PREPROCESSING\n";
	preprocess_pebble();
	std::cout<<"FINISHING PREPROCESSING\n";
	for(auto output:po)
	{
		output_prio.push(std::make_tuple(output,output_stats[output][0]) );
	}

    int output_count = 0;
    int max_cone =  output_stats[std::get<0>(output_prio.top())][0] ;
1;	
/*
    while(!output_prio.empty())
    {
        output_count++;
        int cone_size = output_stats[std::get<0>(output_prio.top())][0]+output_count ;
        max_cone = max_cone < cone_size ? cone_size : max_cone;
        output_prio.pop();   
    } */

	qubits_available = ceil(0.8*(max_cone + actual_po_count) );
	// qubits_available = 2;
    int start_qubits = qubits_available;
	total_steps = 0;
/*	for(auto output:po)
	{
        #ifdef ORDERING_DEBUG 
		std::cout <<"output "<< output <<" has " << output_stats[output][0];
        std::cout <<" nodes \n";
        #endif
		output_prio.push(std::make_tuple(output,output_stats[output][0]) );
	}*/

	std::cout<<"SIZ "<<po.size()<<" "<<output_prio.size()<<"\n";

    for(int i=0; i<qubits_available;i++)
        qubits_bool.push_back(false);
    //qubits_bool.resize(qubits_available);

	outputs_comp_yet = 0;
	qubits_used_yet = 0; // stores total qubits used

	int qubits_till_now = 0; //stores max qubits reqd to compute a cone
	qubits_list = std::stack <int> ();

	for(int i = qubits_available ;i >= 1 ;i--) // set up list of qubits always fill from 1
    {
        qubits_list.push(i); // update again while uncomputing
        qubits_bool[i] = true;
    }

//qubits_used_yet may differ from qubits_till_now as some qubits lost in computing output, they are not computed back.
	std::vector<int> temp;


	int max_cone_label = std::get<0>(output_prio.top());
    std::ofstream debugStats ("debugStats.txt",std::ios::out|std::ios::app);
    debugStats << " \n--------------------------\n";
    int po_processed = 0;
	while(!output_prio.empty())
	{   
        po_processed++;
		int output_label = std::get<0>(output_prio.top());
		temp.push_back(output_label);

		output_prio.pop();

		current_steps = 0;
		current_qubits = 0; // qubits reqd to compute current cone

        debugStats << output_label << ","<< output_stats[output_label][0] << "," << qubits_used_yet << "," ;
		std::cout<<"\nCOMPUTING output: "<<output_label<<"(" << output_stats[output_label][0] <<") " ;
        std::cout << "[" << po_processed << "/" << po.size() << "]\n"; 
		pebble_compute(output_label);
		std::cout<<"\nFINISHED COMPUTING : "<<output_label<<"\n";

		qubits_till_now = current_qubits > qubits_till_now ? current_qubits : qubits_till_now;

		std::cout<<"TOTAL QUBITS USED YET : "<<qubits_used_yet<<"\n\n";
debugStats << "C," << current_qubits << "," << qubits_used_yet << "," ;     

		std::cout<<"\nSTARTING UNCOMPUTATION OF "<<output_label<<"\n\n";
		pebble_uncompute(output_label);
		std::cout<<"\nFINISHED UNCOMPUTING : "<<output_label<<"\n\n";
debugStats << "U," << current_qubits << "," << qubits_used_yet << "\n" ;     
		qubits_till_now = current_qubits > qubits_till_now ? current_qubits : qubits_till_now;

		std::cout<<"TOTAL QUBITS USED YET : "<<qubits_used_yet<<"\n\n";
		// std::cout<<"TOTAL QUBITS TILL NOW : "<<qubits_till_now<<" "<<current_qubits<<"\n\n";

		outputs_comp_yet++;
		total_steps += current_steps;

		output_stats[output_label][3] = current_steps;
		output_stats[output_label][2] = current_qubits;	
	}
	std::cout<<"TOTAL STEPS : "<<total_steps<<"\n\n";
    /* iterate second time to reduce steps by using actual number of qubits 
        this may lead to more qubits being used */
    /*
    int initial_step, initial_qubits;
    initial_step = total_steps;
    initial_qubits = qubits_used_yet; 
	Graph::vertex_iterator vIt, vItEnd;
	boost::tie(vIt,vItEnd) = vertices(lutGraph);
	int idCount = 0;
	for(; vIt != vItEnd; vIt++)
	{
	    lutGraph[*vIt].qubit_alloc = -1;
    }
    // reset internal variables 
    computedSet = std::set<int>();
    computedOutput = std::set<int>();
    computedOnce = std::set<int>();
    qubits_available = qubits_used_yet;
	total_steps = 0;

	qubits_bool.resize(qubits_available);

	outputs_comp_yet = 0;
	qubits_used_yet = 0; // stores total qubits used

	qubits_till_now = 0; //stores max qubits reqd to compute a cone
	qubits_list = std::stack <int> ();

	for(int i = qubits_available ;i >= 1 ;i--) // set up list of qubits always fill from 1
    {
        qubits_list.push(i); // update again while uncomputing
        qubits_bool[i] = true;
    } 

    debugStats << " \n--------------------------\n";
    po_processed = 0;
    pebble_steps = std::vector<std::tuple<int,int,int,int>>();
    std::cout << "Recomputing with " << qubits_available << "\n";
    for(auto output_label:temp)
    {
        po_processed++;

		current_steps = 0;
		current_qubits = 0; // qubits reqd to compute current cone

        debugStats << output_label << ","<< output_stats[output_label][0] << "," << qubits_used_yet << "," ;
		std::cout<<"\nCOMPUTING output: "<<output_label<<"(" << output_stats[output_label][0] <<") " ;
        std::cout << "[" << po_processed << "/" << po.size() << "]\n"; 
		pebble_compute(output_label);
		std::cout<<"\nFINISHED COMPUTING : "<<output_label<<"\n";

		qubits_till_now = current_qubits > qubits_till_now ? current_qubits : qubits_till_now;

		std::cout<<"TOTAL QUBITS USED YET : "<<qubits_used_yet<<"\n\n";
debugStats << "C," << current_qubits << "," << qubits_used_yet << "," ;     

		std::cout<<"\nSTARTING UNCOMPUTATION OF "<<output_label<<"\n\n";
		pebble_uncompute(output_label);
		std::cout<<"\nFINISHED UNCOMPUTING : "<<output_label<<"\n\n";
debugStats << "U," << current_qubits << "," << qubits_used_yet << "\n" ;     
		qubits_till_now = current_qubits > qubits_till_now ? current_qubits : qubits_till_now;

		std::cout<<"TOTAL QUBITS USED YET : "<<qubits_used_yet<<"\n\n";
		// std::cout<<"TOTAL QUBITS TILL NOW : "<<qubits_till_now<<" "<<current_qubits<<"\n\n";

		outputs_comp_yet++;
		total_steps += current_steps;

		output_stats[output_label][3] = current_steps;
		output_stats[output_label][2] = current_qubits;	
	}
    debugStats << "Available Q:"<< start_qubits << " " << max_cone  ;
    debugStats << " Initial Q:" << initial_qubits ;
    debugStats << " S:" << initial_step << "\n";
    debugStats << "Final Q:" << qubits_used_yet ;
    debugStats << "S_recomp:" << total_steps << "\n"; 
	

    std::cout<<"STATS "<<temp.size()<<"\n";
	// for(uint i=0; i < temp.size(); i++)
		// std::cout<<temp[i]<<"\n";
	for(uint i=0; i < temp.size(); i++)
	{
		// std::cout<<"ABC\n";
		std::vector<int> v = output_stats[temp[i]];
		std::cout<< temp[i]<< " size: "<<v[0]<<" levels: "<<v[1]<<" qubits: "<<v[2]<<" steps: "<<v[3]<<" \n";
	}
	std::cout<<" s : "<<output_stats[max_cone_label][0]<<" pi : "<<pi.size()<<" po : "<<po.size()<<"\n"; 
    bool redundancy_remove = true; 
    if(redundancy_remove)
    {
    int step_count = pebble_steps.size();
    int c_node, c_qubit, c_type;
    int o_node, o_qubit, o_type;
    std::tuple<int,int,int,int> c;
    std::tuple<int,int,int,int> o;
    for(int i = 0; i < step_count-1; i ++ )
    {
        o = pebble_steps[i];
        c = pebble_steps[i+1]; 
;
#ifdef ORDERING_DEBUG 
	 	std::cout<< "i" << i << " " << std::get<0>(pebble_steps[i])<<" "<<std::get<1>(pebble_steps[i])<<" "<<std::get<2>(pebble_steps[i])<<" "<<std::get<3>(pebble_steps[i])<<"\n";
	       
	 	std::cout<<std::get<0>(pebble_steps[i+1])<<" "<<std::get<1>(pebble_steps[i+1])<<" "<<std::get<2>(pebble_steps[i+1])<<" "<<std::get<3>(pebble_steps[i+1])<<"\n";
#endif
       
        while(!(std::get<1>(o) != std::get<1>(c) || std::get<2>(o) != std::get<2>(c)))
        {
           if((std::get<0>(o) == 1 && std::get<0>(c) == 0) ||
                  (std::get<0>(o) == 0 && std::get<0>(c) == 1))
           {
                pebble_steps.erase(pebble_steps.begin()+i+1);
                pebble_steps.erase(pebble_steps.begin()+i);
                i = i -1;
                if(i < 0)
                    break;
           } 
           else 
               break;
           if(i+1 >= pebble_steps.size() -1 )
               break;
           o = pebble_steps[i];
           c = pebble_steps[i+1];
#ifdef ORDERING_DEBUG 
	 	std::cout<< "wi" << i << " " << std::get<0>(pebble_steps[i])<<" "<<std::get<1>(pebble_steps[i])<<" "<<std::get<2>(pebble_steps[i])<<" "<<std::get<3>(pebble_steps[i])<<"\n";
	       
	 	std::cout<<std::get<0>(pebble_steps[i+1])<<" "<<std::get<1>(pebble_steps[i+1])<<" "<<std::get<2>(pebble_steps[i+1])<<" "<<std::get<3>(pebble_steps[i+1])<<"\n";
#endif 
        }
        if(i+1 >= pebble_steps.size() -1 )
               break;
          
    }
    std::cout << "Original step count: " << step_count << "\n";
    std::cout << "Optimized step count:" << pebble_steps.size() << "\n";
   debugStats << " S_opt:" << pebble_steps.size() << "\n"; 
    }
    */
    debugStats.close();
    while(!pebble_prio.empty())
    {
        pebble_prio.pop();
    }
}

}} // end of namespace legacy cirkit 
