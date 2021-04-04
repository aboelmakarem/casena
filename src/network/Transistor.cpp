// Transistor.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/28/2021

#include "Transistor.h"

namespace CASENA
{
	Transistor::Transistor(){Initialize();}
	Transistor::Transistor(const Transistor& transistor) : Component(transistor){*this = transistor;}
	Transistor::~Transistor(){Reset();}
	Transistor& Transistor::operator=(const Transistor& transistor)
	{
		Component::operator=(transistor);
		name = transistor.name;
		return *this;
	}
	void Transistor::Reset()
	{
		Initialize();
		Component::Reset();
	}
	bool Transistor::IsTransistor(const char* line)
	{
		
	}
	unsigned int Transistor::ReadMaxNodeID(const char* line)
	{
		
	}
	bool Transistor::Read(const char* line)
	{
		
	}
	void Transistor::Initialize()
	{
		name.Reset();
	}	

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
		return *this;
	}
	void MOSFET::Reset()
	{
		Initialize();
		Transistor::Reset();
	}
	unsigned int MOSFET::ClaimIDs(const unsigned int& start_id)
	{
		// a MOSFET owns 6 branches with the following IDs in order
		// gate-source branch: id
		// gate-drain branch: id + 1
		// gate-body branch: id + 2
		// source-body branch: id + 3
		// drain-body branch: id + 4
		// drain-source branch: id + 5
		ID(start_id);
		return start_id + 6;
	}
	void MOSFET::Equation(EZ::Math::Matrix& f) const
	{
		
	}
	void MOSFET::Gradients(EZ::Math::Matrix& A) const
	{
		
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
	}
	bool MOSFET::ReadProperties(const EZ::List<EZ::String*>* line_tokens)
	{
		
	}
}


