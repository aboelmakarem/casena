// Network.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#include "stdio.h"
#include "Network.h"

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
	void Network::Reset(){Initialize();}
	bool Network::SteadyState(){return (analysis_type == 1);}
	bool Network::Transient(){return (analysis_type == 2);}
	void Network::Initialize()
	{
		analysis_type = 0;
		components.Reset();
	}
}

int main(int argc,char** argv)
{
	
	return 0;
}

