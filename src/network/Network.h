// Network.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/14/2021

#ifndef NETWORK_H_
#define NETWORK_H_

#include "Node.h"
#include "Branch.h"

namespace CASENA
{
	class Network
	{
	public:
		static Network* GetNetwork();
		static void ReleaseNetwork();
		~Network();
		void Reset();
		static bool SteadyState();
		static bool Transient();
		
	private:
		static Network* network;
		Network();
		Network(const Network& network);
		Network& operator=(const Network& network);
		void Initialize();
		static int analysis_type;
		EZ::List<Component*> components;
	};
}

#endif

