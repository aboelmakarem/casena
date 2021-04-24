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
		if(strncmp(line,"mosfet",6) == 0)	return MOSFET::CurrentCount();
		if(strncmp(line,"bjt",3) == 0)		return BJT::CurrentCount();
		return 0;
	}
	Component* Transistor::Create(const char* line)
	{
		if(strncmp(line,"mosfet",6) == 0)	return new MOSFET;
		if(strncmp(line,"bjt",3) == 0)		return new BJT;
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
		printf("passing transistor RP\n");
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
		memcpy(currents,transistor.currents,HistoryCount*MOSFET_CURRENT_COUNT*sizeof(double));
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
		return (start_id + MOSFET_CURRENT_COUNT);
	}
	void MOSFET::Equation(EZ::Math::Matrix& f) const
	{
		printf("eqn mosfet\n");
		unsigned int id = ID();
		// convention: all currents are flowing into the transistor (gate, 
		// drain, source, body)
		if(Network::SteadyState())
		{
			// for steady state operation, gate and body currents 
			// are zero and source current is equal to drain current 
			// but in an opposite direction
			double i_drain = SteadyStateDrainCurrent();
			printf("got id = %e\n",i_drain);
			f(drain_node_id,0,f(drain_node_id,0) + i_drain);
			f(source_node_id,0,f(source_node_id,0) - i_drain);
			f(id,0,0.0);
			f(id + 1,0,currents[HistoryCount] - i_drain);
			// source current is the same as drain current but with an 
			// opposite sign
			f(id + 2,0,currents[2*HistoryCount] + i_drain);
			f(id + 3,0,0.0);
		}
		else
		{
			//TransientCurrents(i_gate,i_drain,i_source,i_body);
		}
		printf("done eqn mosfet\n");fflush(0);
	}
	void MOSFET::Gradients(EZ::Math::Matrix& A) const
	{
		printf("grad mosfet\n");
		unsigned int id = ID();
		if(Network::SteadyState())
		{
			printf("set grad @ SS : %u,%u,%u,%u\n",
			drain_node_id,source_node_id,gate_node_id,body_node_id);
			A(drain_node_id,id + 1,1.0);
			A(source_node_id,id + 2,1.0);
			// gate branch equation: ig = 0
			A(id,id,1.0);
			// drain branch equation: is = -id or is + id = 0			
			A(id + 1,id + 1,1.0);
			double vg_grad = 0.0;
			double vd_grad = 0.0;
			double vs_grad = 0.0;
			double vb_grad = 0.0;
			SteadyStateDrainCurrentGradients(vg_grad,vd_grad,vs_grad,vb_grad);
			A(id + 1,gate_node_id,-vg_grad);
			A(id + 1,drain_node_id,-vd_grad);
			A(id + 1,source_node_id,-vs_grad);
			A(id + 1,body_node_id,-vb_grad);
			// source branch equation: is = -id or is + id = 0
			A(id + 2,id + 1,1.0);
			A(id + 2,id + 2,1.0);
			// body branch equation: ib = 0
			A(id + 3,id + 3,1.0);
		}
		else
		{
			
		}
	}
	void MOSFET::Update(const EZ::Math::Matrix& x,const unsigned int& id_offset)
	{
		printf("updating mosfet\n");
		unsigned int id = ID() - id_offset;
		unsigned int target_index = 0;
		unsigned int source_index = 0;
		for(unsigned int i = HistoryCount - 1 ; i > 0 ; i--)
		{
			target_index = i*MOSFET_CURRENT_COUNT;
			source_index = (i - 1)*MOSFET_CURRENT_COUNT;
			memcpy(&currents[target_index],&currents[source_index],MOSFET_CURRENT_COUNT*sizeof(double));
		}
		currents[0] = x(id,0);
		currents[1] = x(id + 1,0);
		currents[2] = x(id + 2,0);
		currents[3] = x(id + 3,0);
		printf("done updating\n");fflush(0);
	}
	void MOSFET::Print() const
	{
		if(np_type == 1)		printf("nmos %s\n",Name()());
		else if(np_type == 2)	printf("pmos %s\n",Name()());
		printf("\tnodes (G/D/S/B): %u,%u,%u,%u\n",gate_node_id,drain_node_id,source_node_id,body_node_id);
		printf("\tthreshold voltage = %e, Early voltage = %e\n",threshold_voltage,early_voltage);
		printf("\tW = %e, L = %e, Cox = %e, mu = %e\n",w,l,cox,mu);
		printf("\tgamma = %e, phi = %e\n",gamma,phi);
	}
	unsigned int MOSFET::CurrentCount(){return MOSFET_CURRENT_COUNT;}
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
		memset(currents,0,HistoryCount*MOSFET_CURRENT_COUNT*sizeof(double));
	}
	bool MOSFET::ReadProperties(const EZ::List<EZ::String*>* line_tokens)
	{
		if(line_tokens->Size() != 13)		return false;
		EZ::ListItem<EZ::String*>* token = line_tokens->Start();
		np_type = atoi((*(token->Data()))());
		token = token->Next();
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
	double MOSFET::SteadyStateDrainCurrent() const
	{
		printf("SSDC 00\n");fflush(0);
		double i_drain = 0.0;
		Network* network = Network::GetNetwork();
		// node voltages
		double vg = network->GetNode(gate_node_id)->Voltage();
		double vd = network->GetNode(drain_node_id)->Voltage();
		double vs = network->GetNode(source_node_id)->Voltage();
		double vb = network->GetNode(body_node_id)->Voltage();
		double vt = threshold_voltage + gamma*(sqrt(2.0*phi + fabs(vs - vb)) - sqrt(2.0*phi));
		double vgs = vg - vs;
		double vds = vd - vs;
		printf("SSDC 01\n");fflush(0);
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
					// triode mode operation
					i_drain = mu*cox*w/l*(vov*vds - 0.5*vds*vds);
				}
			}
		}
		else if(np_type == 2)
		{
			// PMOS transistor
			
		}
		printf("SSDC 05\n");fflush(0);
		return i_drain;
	}
	void MOSFET::TransientCurrents(double& i_gate,double& i_drain,double& i_source,double& i_body) const
	{
		// fix this
		i_gate = 0.0;
		i_body = 0.0;
		i_drain = 0.0;
		i_source = i_drain;
	}
	void MOSFET::SteadyStateDrainCurrentGradients(double& vg_grad,double& vd_grad,double& vs_grad,double& vb_grad) const
	{
		vg_grad = 0.0;
		vd_grad = 0.0;
		vs_grad = 0.0;
		vb_grad = 0.0;
		Network* network = Network::GetNetwork();
		// node voltages
		double vg = network->GetNode(gate_node_id)->Voltage();
		double vd = network->GetNode(drain_node_id)->Voltage();
		double vs = network->GetNode(source_node_id)->Voltage();
		double vb = network->GetNode(body_node_id)->Voltage();
		double vt = threshold_voltage + gamma*(sqrt(2.0*phi + fabs(vs - vb)) - sqrt(2.0*phi));
		double sqrt_term = sqrt(2.0*phi + fabs(vs - vb));
		double dvtdvs = 0.5*gamma/sqrt_term;
		if(vs < vb) 	dvtdvs = -dvtdvs;
		double dvtdvb = -dvtdvs;
		double vgs = vg - vs;
		double vds = vd - vs;
		if(np_type == 1)
		{
			// NMOS transistor
			double vov = vgs - vt;
			// for cut-off mode operation, all gradients are zero
			if(vov < 0.0)					return;
			else
			{
				if(vds >= vov)
				{
					// saturation mode operation
					double factor = mu*cox*w/l;
					vg_grad = factor*vov*(1.0 + vds/early_voltage);
					vd_grad = 0.5*factor*vov*vov*(1.0/early_voltage);
					vs_grad = -factor*vov*((1.0 + vds/early_voltage)*(1.0 + dvtdvs) + 0.5*vov*(1.0/early_voltage));
					vb_grad = -factor*vov*(1.0 + vds/early_voltage)*dvtdvb;
				}
				else
				{
					// triode mode operation
					mu*cox*w/l*(vov*vds - 0.5*vds*vds);
					double factor = mu*cox*w/l;
					vg_grad = factor*vds;
					vd_grad = factor*(vov - vds);
					vs_grad = -factor*(vov + vds*dvtdvs);
					vb_grad = -factor*vds*dvtdvb;
				}
			}
		}
		else if(np_type == 2)
		{
			// PMOS transistor
			
		}
	}
	
	BJT::BJT(){Initialize();}
	BJT::BJT(const BJT& transistor) : Transistor(transistor){*this = transistor;}
	BJT::~BJT(){Reset();}
	BJT& BJT::operator=(const BJT& transistor)
	{
		Transistor::operator=(transistor);
		base_node_id = transistor.base_node_id;
		emitter_node_id = transistor.emitter_node_id;
		collector_node_id = transistor.collector_node_id;
		np_type = transistor.np_type;
		memcpy(currents,transistor.currents,HistoryCount*BJT_CURRENT_COUNT*sizeof(double));
		return *this;
	}
	void BJT::Reset()
	{
		Initialize();
		Transistor::Reset();
	}
	unsigned int BJT::ClaimIDs(const unsigned int& start_id)
	{
		// a BJT owns 3 branches with the following IDs in order
		// base branch: id
		// emitter branch: id + 1
		// collector branch: id + 2
		ID(start_id);
		return (start_id + BJT_CURRENT_COUNT);
	}
	void BJT::Equation(EZ::Math::Matrix& f) const
	{
		
	}
	void BJT::Gradients(EZ::Math::Matrix& A) const
	{
		
	}
	void BJT::Update(const EZ::Math::Matrix& x,const unsigned int& id_offset)
	{
		
	}
	void BJT::Print() const
	{
		
	}
	unsigned int BJT::CurrentCount(){return BJT_CURRENT_COUNT;}
	void BJT::Initialize()
	{
		np_type = 0;
		base_node_id = 0;
		emitter_node_id = 0;
		collector_node_id = 0;
		memset(currents,0,HistoryCount*BJT_CURRENT_COUNT*sizeof(double));
	}
	bool BJT::ReadProperties(const EZ::List<EZ::String*>* line_tokens)
	{
		
	}
}


