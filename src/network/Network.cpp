// Network.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#include "stdio.h"
#include "Network.h"
#include "Transistor.h"
#include "InputOutput.h"

namespace CASENA
{
	Network* Network::network = 0;
	int Network::analysis_type = 0;
	Network* Network::GetNetwork()
	{
		if(network == 0)		network = new Network;
		return network;
	}
	void Network::ReleaseNetwork()
	{
		if(network != 0)
		{
			delete network;
			network = 0;
		}
	}
	Network::Network(){Initialize();}
	Network::Network(const Network& network){*this = network;}
	Network::~Network(){Reset();}
	Network& Network::operator=(const Network& network)
	{
		// do nothing, networks cannot be equated
		analysis_type = network.analysis_type;
		return *this;
	}
	void Network::Reset()
	{
		nodes.ForEachDelete();
		components.ForEachDelete();
		Initialize();
	}
	bool Network::Read(const char* filename)
	{
		Reset();
		FILE* file = fopen(filename,"r");
		if(file == 0)		return false;
		char line[2048];
		// get maximum node id and component counts
		unsigned int max_node_id = 0;
		unsigned int node_id = 0;
		unsigned int branch_count = 0;
		while(!feof(file))
		{
			if(!EZ::IO::ReadLine(line,2048,file))		continue;
			if(strlen(line) < 3)		continue;
			if(Branch::IsBranchComponent(line))
			{
				node_id = Branch::ReadMaxNodeID(line);
				branch_count++;
			}
			else if(Transistor::IsTransistor(line))
			{
				node_id = Transistor::ReadMaxNodeID(line);
				branch_count += Transistor::BranchCount(line);
			}
			if(node_id > max_node_id)		max_node_id = node_id;
		}
		// create node array, node IDs are zero based
		nodes.Allocate(max_node_id + 1);
		for(unsigned int i = 0 ; i <= max_node_id ; i++)
		{
			nodes[i] = new Node;
			nodes[i]->ID(i);
		}
		unsigned int id = max_node_id + 1;
		first_branch_id = id;
		// read components
		fseek(file,0,SEEK_SET);
		Component* component = 0;
		unsigned int start_id = 0;
		branches.Allocate(branch_count);
		while(!feof(file))
		{
			if(!EZ::IO::ReadLine(line,2048,file))		continue;
			if(strlen(line) < 3)		continue;
			component = 0;
			start_id = id;
			if(Branch::IsBranchComponent(line))
			{
				component = Branch::Create(line);
				if(!component->Read(line))
				{
					delete component;
					continue;
				}
				id = component->ClaimIDs(id);
				branches[start_id] = (Branch*)component;
				components.PushBack(component);
			}
			else if(Transistor::IsTransistor(line))
			{
				component = Transistor::Create(line);
				if(!component->Read(line))
				{
					delete component;
					continue;
				}
				id = component->ClaimIDs(id);
				// skip as many branches in the branches array as 
				// claimed by the transistor
				for(unsigned int i = start_id ; i < id ; i++)
				{
					branches[i] = 0;
				}
				components.PushBack(component);
			}
		}
		fclose(file);
		return true;
	}
	bool Network::SteadyState(){return (analysis_type == 1);}
	bool Network::Transient(){return (analysis_type == 2);}
	Node* Network::GetNode(const unsigned int& id) const{return nodes[id];}
	Branch* Network::GetBranch(const unsigned int& id) const{return branches[id - first_branch_id];}
	void Network::Print() const
	{
		printf("network:\n");
		for(EZ::ListItem<Component*>* item = components.Start() ; item != 0 ; item = item->Next())
		{
			item->Data()->Print();
		}
	}
	void Network::Initialize()
	{
		analysis_type = 0;
		nodes.Reset();
		components.Reset();
		first_branch_id = 0;
		branches.Reset();
	}
}

int main(int argc,char** argv)
{
	if(argc < 2)
	{
		printf("error: missing input file name\n");
		printf("usage: casena input_file_name\n");
		return 1;
	}
	CASENA::Network* network = CASENA::Network::GetNetwork();
	if(!network->Read(argv[1]))
	{
		printf("error: failed to read input file %s\n",argv[1]);
	}
	network->Print();
	CASENA::Network::ReleaseNetwork();
	return 0;
}

