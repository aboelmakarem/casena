// Network.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/14/2021

#ifndef NETWORK_H_
#define NETWORK_H_

#include "Node.h"
#include "Branch.h"
#include "Array.h"
#include "List.h"

namespace CASENA
{
	class Network
	{
	public:
		static Network* GetNetwork();
		static void ReleaseNetwork();
		~Network();
		void Reset();
		bool Read(const char* filename);
		static bool SteadyState();
		static bool Transient();
		Node* GetNode(const unsigned int& id) const;
		Branch* GetBranch(const unsigned int& id) const;
		void Print() const;
		void Run() const;
		
	private:
		static Network* network;
		Network();
		Network(const Network& network);
		Network& operator=(const Network& network);
		void Initialize();
		static int analysis_type;
		EZ::Array<Node*> nodes;
		EZ::List<Component*> components;
		unsigned int first_branch_id;
		EZ::Array<Branch*> branches;
	};
}

#endif

