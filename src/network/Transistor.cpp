// Transistor.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/28/2021

#include "Transistor.h"
#include "Network.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"

namespace CASENA
{
	Transistor::Transistor(){Initialize();}
	Transistor::Transistor(const Transistor& transistor) : Component(transistor){*this = transistor;}
	Transistor::~Transistor(){Reset();}
	Transistor& Transistor::operator=(const Transistor& transistor)
	{
		Component::operator=(transistor);
		return *this;
	}
	void Transistor::Reset()
	{
		Initialize();
		Component::Reset();
	}
	bool Transistor::IsTransistor(const char* line)
	{
		if(strncmp(line,"mosfet",6) == 0)	return true;
		if(strncmp(line,"bjt",3) == 0)		return true;
		return false;
	}
	unsigned int Transistor::ReadMaxNodeID(const char* line)
	{
		EZ::String line_string(line);
		EZ::List<EZ::String*> tokens;
		line_string.Tokenize(tokens," \t");
		// read component type
		EZ::ListItem<EZ::String*>* token = tokens.Start();
		delete token->Data();
		tokens.PopFront();
		// read component name
		token = tokens.Start();
		delete token->Data();
		tokens.PopFront();
		// read node IDs
		token = tokens.Start();
		unsigned int node1_id = atoi((*(token->Data()))());
		delete token->Data();
		tokens.PopFront();
		token = tokens.Start();
		unsigned int node2_id = atoi((*(token->Data()))());
		delete token->Data();
		tokens.PopFront();
		token = tokens.Start();
		unsigned int node3_id = atoi((*(token->Data()))());
		delete token->Data();
		tokens.PopFront();
		unsigned int node4_id = 0;
		if(strncmp(line,"mosfet",6) == 0)
		{
			// This is a MOSFET, read the 4th node ID
			token = tokens.Start();
			node4_id = atoi((*(token->Data()))());
			delete token->Data();
			tokens.PopFront();
		}
		for(EZ::ListItem<EZ::String*>* item = tokens.Start() ; item != 0 ; item = item->Next())
		{
			delete item->Data();
		}
		tokens.Reset();
		unsigned int max_node_id = (node1_id > node2_id ? node1_id:node2_id);
		max_node_id = (max_node_id > node3_id ? max_node_id:node3_id);
		max_node_id = (max_node_id > node4_id ? max_node_id:node4_id);
		return max_node_id;
	}
	unsigned int Transistor::BranchCount(const char* line)
	{
		if(strncmp(line,"mosfet",6) == 0)	return 4;
		if(strncmp(line,"bjt",3) == 0)		return 3;
		return 0;
	}
	Component* Transistor::Create(const char* line)
	{
		if(strncmp(line,"mosfet",6) == 0)	return new MOSFET;
		if(strncmp(line,"bjt",3) == 0)		return new MOSFET;
		return 0;
	}
	bool Transistor::Read(const char* line)
	{
		EZ::String line_string(line);
		EZ::List<EZ::String*> tokens;
		line_string.Tokenize(tokens," \t");
		// read component type
		EZ::ListItem<EZ::String*>* token = tokens.Start();
		delete token->Data();
		tokens.PopFront();
		// read component name
		token = tokens.Start();
		Name(*(token->Data()));
		delete token->Data();
		tokens.PopFront();
		// read the rest of the string and pass it to properties reader
		bool read_properties = ReadProperties(&tokens);
		if(!read_properties)
		{
			printf("error: incorrect definition for transistor %s\n",Name()());
		}
		for(EZ::ListItem<EZ::String*>* item = tokens.Start() ; item != 0 ; item = item->Next())
		{
			delete item->Data();
		}
		tokens.Reset();
		return read_properties;
	}
	void Transistor::Initialize(){}

	MOSFET::MOSFET(){Initialize();}
	MOSFET::MOSFET(const MOSFET& transistor) : Transistor(transistor){*this = transistor;}
	MOSFET::~MOSFET(){Reset();}
	MOSFET& MOSFET::operator=(const MOSFET& transistor)
	{
		Transistor::operator=(transistor);
		gate_node_id = transistor.gate_node_id;
		source_node_id = transistor.source_node_id;
		drain_node_id = transistor.drain_node_id;
		body_node_id = transistor.body_node_id;
		np_type = transistor.np_type;
		threshold_voltage = transistor.threshold_voltage;
		early_voltage = transistor.early_voltage;
		mu = transistor.mu;
		cox = transistor.cox;
		w = transistor.w;
		l = transistor.l;
		gamma = transistor.gamma;
		phi = transistor.phi;
		memcpy(currents,transistor.currents,HistoryCount*4*sizeof(double));
		return *this;
	}
	void MOSFET::Reset()
	{
		Initialize();
		Transistor::Reset();
	}
	unsigned int MOSFET::ClaimIDs(const unsigned int& start_id)
	{
		// a MOSFET owns 4 branches with the following IDs in order
		// gate branch: id
		// drain branch: id + 1
		// source branch: id + 2
		// body branch: id + 3
		ID(start_id);
		return start_id + 4;
	}
	void MOSFET::Equation(EZ::Math::Matrix& f) const
	{
		double i_gate = 0.0;
		double i_drain = 0.0;
		double i_source = 0.0;
		double i_body = 0.0;
		unsigned int id = ID();
		// convention: gate and drain current are flowing into the transistor 
		// and source and body currents are flowing out of the transistor
		if(Network::SteadyState())
		{
			SteadyStateCurrents(i_gate,i_drain,i_source,i_body);
			f(gate_node_id,0,f(gate_node_id,0) + i_gate);
			f(drain_node_id,0,f(drain_node_id,0) + i_drain);
			f(source_node_id,0,f(source_node_id,0) - i_source);
			f(body_node_id,0,f(body_node_id,0) - i_body);
			f(id,0,currents[0] - i_gate);
			f(id + 1,0,currents[HistoryCount] - i_drain);
			f(id + 2,0,currents[2*HistoryCount] - i_source);
			f(id + 3,0,currents[3*HistoryCount] - i_body);
		}
		else
		{
			TransientCurrents(i_gate,i_drain,i_source,i_body);
		}
	}
	void MOSFET::Gradients(EZ::Math::Matrix& A) const
	{
		unsigned int id = ID();
		if(Network::SteadyState())
		{
			A(gate_node_id,id,1.0);
			A(drain_node_id,id + 1,1.0);
			A(source_node_id,id + 2,-1.0);
			A(body_node_id,id + 3,-1.0);
			A(id,id,1.0);
			A(id + 1,id + 1,1.0);
			A(id + 2,id + 2,1.0);
			A(id + 3,id + 3,1.0);
		}
		else
		{
			
		}
	}
	void MOSFET::Update(const EZ::Math::Matrix& x,const unsigned int& id_offset)
	{
		
	}
	void MOSFET::Print() const
	{
		printf("mosfet %s\n",Name()());
	}
	void MOSFET::Initialize()
	{
		gate_node_id = 0;
		source_node_id = 0;
		drain_node_id = 0;
		body_node_id = 0;
		np_type = 1;
		threshold_voltage = 0.0;
		early_voltage = 0.0;
		mu = 0.0;
		cox = 0.0;
		w = 0.0;
		l = 0.0;
		gamma = 0.0;
		phi = 0.0;
		memset(currents,0,HistoryCount*4*sizeof(double));
	}
	bool MOSFET::ReadProperties(const EZ::List<EZ::String*>* line_tokens)
	{
		if(line_tokens->Size() != 13)		return false;
		EZ::ListItem<EZ::String*>* token = line_tokens->Start();
		token = token->Next();
		np_type = atoi((*(token->Data()))());
		gate_node_id = atoi((*(token->Data()))());
		token = token->Next();
		source_node_id = atoi((*(token->Data()))());
		token = token->Next();
		drain_node_id = atoi((*(token->Data()))());
		token = token->Next();
		body_node_id = atoi((*(token->Data()))());
		token = token->Next();
		threshold_voltage = atof((*(token->Data()))());
		token = token->Next();
		early_voltage = atof((*(token->Data()))());
		token = token->Next();
		mu = atof((*(token->Data()))());
		token = token->Next();
		cox = atof((*(token->Data()))());
		token = token->Next();
		w = atof((*(token->Data()))());
		token = token->Next();
		l = atof((*(token->Data()))());
		token = token->Next();
		gamma = atof((*(token->Data()))());
		token = token->Next();
		phi = atof((*(token->Data()))());
		return true;
	}
	void MOSFET::SteadyStateCurrents(double& i_gate,double& i_drain,double& i_source,double& i_body) const
	{
		// at steady state, the gate and body currents are zero
		i_gate = 0.0;
		i_body = 0.0;
		Network* network = Network::GetNetwork();
		// node voltages
		double vg = network->GetNode(gate_node_id)->Voltage();
		double vd = network->GetNode(drain_node_id)->Voltage();
		double vs = network->GetNode(source_node_id)->Voltage();
		double vb = network->GetNode(body_node_id)->Voltage();
		double vt = threshold_voltage + gamma*(sqrt(2.0*phi + fabs(vs - vb)) - sqrt(2.0*phi));
		double vgs = vg - vs;
		double vds = vd - vs;
		if(np_type == 1)
		{
			// NMOS transistor
			double vov = vgs - vt;
			if(vov < 0.0)
			{
				// cut-off mode operation
				i_drain = 0.0;
			}
			else
			{
				if(vds >= vov)
				{
					// saturation mode operation
					i_drain = 0.5*mu*cox*w/l*vov*vov*(1.0 + vds/early_voltage);
				}
				else
				{
					// tridoe mode operation
					i_drain = mu*cox*w/l*(vov*vds - 0.5*vds*vds);
				}
			}
		}
		else if(np_type == 2)
		{
			// PMOS transistor
			
		}
		// for steady state, the drain and source currents are the same 
		// but with opposite directions
		i_source = i_drain;
	}
	void MOSFET::TransientCurrents(double& i_gate,double& i_drain,double& i_source,double& i_body) const
	{
		// fix this
		i_gate = 0.0;
		i_body = 0.0;
		i_drain = 0.0;
		i_source = i_drain;
	}
}


