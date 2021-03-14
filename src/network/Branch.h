// Branch.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#ifndef BRANCH_H_
#define BRANCH_H_

#include "Node.h"

namespace CASENA
{
	class Branch : public Component
	{
	public:
		Branch();
		Branch(const Branch& branch);
		virtual ~Branch();
		virtual Branch& operator=(const Branch& branch);
		void Reset();
		double Current() const;
		double PreviousCurrent() const;
		void Current(const double& value);
		void StartNode(const Node* node);
		const Node* StartNode() const;
		void EndNode(const Node* node);
		const Node* EndNode() const;
		
	private:
		void Initialize();
		double currents[HistoryCount];
		const Node* start;
		const Node* end;
	};
	
	class IndSource : public Branch
	{
	public:
		IndSource();
		IndSource(const IndSource& source);
		virtual ~IndSource();
		virtual IndSource& operator=(const IndSource& source);
		void Reset();
		void Form(const double& value);
		double Form() const;
		
	private:
		void Initialize();
		double form;
	};
	
	class IndVolSource : public IndSource
	{
	public:
		IndVolSource();
		IndVolSource(const IndVolSource& source);
		~IndVolSource();
		IndVolSource& operator=(const IndVolSource& source);
		void Reset();
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
	};
	
	class IndCurrSource : public IndSource
	{
	public:
		IndCurrSource();
		IndCurrSource(const IndCurrSource& resistor);
		~IndCurrSource();
		IndCurrSource& operator=(const IndCurrSource& resistor);
		void Reset();
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
	};
	
	class VCSource : public Branch
	{
	public:
		VCSource();
		VCSource(const VCSource& source);
		virtual ~VCSource();
		virtual VCSource& operator=(const VCSource& source);
		void Reset();
		void Coefficient(const double& value);
		double Coefficient() const;
		void SourceStartNode(const Node* node);
		const Node* SourceStartNode() const;
		void SourceEndNode(const Node* node);
		const Node* SourceEndNode() const;
		double SourceVoltage() const;
		
	private:
		void Initialize();
		double coefficient;
		const Node* source_start;
		const Node* source_end;
	};
	
	class VCVolSource : public VCSource
	{
	public:
		VCVolSource();
		VCVolSource(const VCVolSource& source);
		~VCVolSource();
		VCVolSource& operator=(const VCVolSource& source);
		void Reset();
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
	};
	
	class VCCurrSource : public VCSource
	{
	public:
		VCCurrSource();
		VCCurrSource(const VCCurrSource& source);
		~VCCurrSource();
		VCCurrSource& operator=(const VCCurrSource& source);
		void Reset();
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
	};
	
	class CCSource : public Branch
	{
	public:
		CCSource();
		CCSource(const CCSource& source);
		virtual ~CCSource();
		virtual CCSource& operator=(const CCSource& source);
		void Reset();
		void Coefficient(const double& value);
		double Coefficient() const;
		void SourceBranch(const Branch* branch);
		const Branch* SourceBranch() const;
		double SourceCurrent() const;
		
	private:
		void Initialize();
		double coefficient;
		const Branch* source_branch;
	};
	
	class CCVolSource : public CCSource
	{
	public:
		CCVolSource();
		CCVolSource(const CCVolSource& source);
		~CCVolSource();
		CCVolSource& operator=(const CCVolSource& source);
		void Reset();
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
	};
	
	class CCCurrSource : public CCSource
	{
	public:
		CCCurrSource();
		CCCurrSource(const CCCurrSource& source);
		~CCCurrSource();
		CCCurrSource& operator=(const CCCurrSource& source);
		void Reset();
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
	};
	
	class Resistor : public Branch
	{
	public:
		Resistor();
		Resistor(const Resistor& resistor);
		~Resistor();
		Resistor& operator=(const Resistor& resistor);
		void Reset();
		void Resistance(const double& value);
		double Resistance() const;
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
		double resistance;
	};
	
	class Capacitor : public Branch
	{
	public:
		Capacitor();
		Capacitor(const Capacitor& capacitor);
		~Capacitor();
		Capacitor& operator=(const Capacitor& capacitor);
		void Reset();
		void Capacitance(const double& value);
		double Capacitance() const;
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
		double capacitance;
	};
	
	class Inductor : public Branch
	{
	public:
		Inductor();
		Inductor(const Inductor& inductor);
		~Inductor();
		Inductor& operator=(const Inductor& inductor);
		void Reset();
		void Inductance(const double& value);
		double Inductance() const;
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
		double inductance;
	};
}

#endif

