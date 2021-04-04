// Transistor.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/28/2021

#ifndef TRANSISTOR_H_
#define TRANSISTOR_H_

#include "Component.h"
#include "List.h"
#include "String.h"

namespace CASENA
{
	class Transistor : public Component
	{
	public:
		Transistor();
		Transistor(const Transistor& transistor);
		virtual ~Transistor();
		virtual Transistor& operator=(const Transistor& transistor);
		void Reset();
		static bool IsTransistor(const char* line);
		static unsigned int ReadMaxNodeID(const char* line);
		bool Read(const char* line);
		static Component* Create(const char* line);
		
	private:
		void Initialize();
		EZ::String name;
		
	protected:
		virtual bool ReadProperties(const EZ::List<EZ::String*>* line_tokens) = 0;
	};
	
	class MOSFET : public Transistor
	{
	public:
		MOSFET();
		MOSFET(const MOSFET& transistor);
		~MOSFET();
		MOSFET& operator=(const MOSFET& transistor);
		void Reset();
		unsigned int ClaimIDs(const unsigned int& start_id);
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
		bool ReadProperties(const EZ::List<EZ::String*>* line_tokens);
		// np_type: 1 for NMOS and 2 for PMOS, NMOS is default
		int np_type;
		unsigned int gate_node_id;
		unsigned int source_node_id;
		unsigned int drain_node_id;
		unsigned int body_node_id;
		double threshold_voltage;
		double early_voltage;
		double mu;
		double cox;
		double w;
		double l;
	};
}

#endif

