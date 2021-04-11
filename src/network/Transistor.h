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
		static unsigned int BranchCount(const char* line);
		static Component* Create(const char* line);
		bool Read(const char* line);
		
	private:
		void Initialize();
		
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
		void Print() const;
		
	private:
		void Initialize();
		bool ReadProperties(const EZ::List<EZ::String*>* line_tokens);
		void SteadyStateCurrents(double& i_gate,double& i_drain,double& i_source,double& i_body) const;
		void TransientCurrents(double& i_gate,double& i_drain,double& i_source,double& i_body) const;
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
		double gamma;
		double phi;
		double currents[4*HistoryCount];
	};
	
	class BJT : public Transistor
	{
	public:
		BJT();
		BJT(const BJT& transistor);
		~BJT();
		BJT& operator=(const BJT& transistor);
		void Reset();
		unsigned int ClaimIDs(const unsigned int& start_id);
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		void Print() const;
		
	private:
		void Initialize();
		bool ReadProperties(const EZ::List<EZ::String*>* line_tokens);
		// np_type: 1 for npn BJT and 2 for pnp BJT, npn BJT is default
		int np_type;
		unsigned int base_node_id;
		unsigned int emitter_node_id;
		unsigned int collector_node_id;
		double currents[3*HistoryCount];
	};
}

#endif

