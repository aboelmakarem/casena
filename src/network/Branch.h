// Branch.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#ifndef BRANCH_H_
#define BRANCH_H_

#include "Component.h"

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
		static bool IsBranchComponent(const char* line);
		static unsigned int ReadMaxNodeID(const char* line);
		bool Read(const char* line);
		static Component* Create(const char* line);
		double Current() const;
		double PreviousCurrent() const;
		void Current(const double& value);
		void StartNodeID(const unsigned int& target_id);
		unsigned int StartNodeID() const;
		void EndNodeID(const unsigned int& target_id);
		unsigned int EndNodeID() const;
		void Equation(EZ::Math::Matrix& f) const;
		void Gradients(EZ::Math::Matrix& A) const;
		
	private:
		void Initialize();
		double currents[HistoryCount];
		
	protected:
		virtual bool ReadProperties(const EZ::List<EZ::String*>* line_tokens);
		virtual void BranchEquation(EZ::Math::Matrix& f) const = 0;
		virtual void BranchGradients(EZ::Math::Matrix& A) const = 0;
		unsigned int start_node_id;
		unsigned int end_node_id;
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
		bool ReadProperties(const EZ::List<EZ::String*>* line_tokens);
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
		void Print() const;
		
	private:
		void Initialize();
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
	};
	
	class IndCurrSource : public IndSource
	{
	public:
		IndCurrSource();
		IndCurrSource(const IndCurrSource& resistor);
		~IndCurrSource();
		IndCurrSource& operator=(const IndCurrSource& resistor);
		void Reset();
		void Print() const;
		
	private:
		void Initialize();
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
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
		void SourceStartNodeID(const unsigned int& id);
		unsigned int SourceStartNodeID() const;
		void SourceEndNodeID(const unsigned int& id);
		unsigned int SourceEndNodeID() const;
		double SourceVoltage() const;
		
	private:
		void Initialize();
		bool ReadProperties(const EZ::List<EZ::String*>* line_tokens);
		double coefficient;
		
	protected:
		unsigned int source_start_node_id;
		unsigned int source_end_node_id;
	};
	
	class VCVolSource : public VCSource
	{
	public:
		VCVolSource();
		VCVolSource(const VCVolSource& source);
		~VCVolSource();
		VCVolSource& operator=(const VCVolSource& source);
		void Reset();
		void Print() const;
		
	private:
		void Initialize();
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
	};
	
	class VCCurrSource : public VCSource
	{
	public:
		VCCurrSource();
		VCCurrSource(const VCCurrSource& source);
		~VCCurrSource();
		VCCurrSource& operator=(const VCCurrSource& source);
		void Reset();
		void Print() const;
		
	private:
		void Initialize();
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
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
		void SourceBranchID(const unsigned int& id);
		unsigned int SourceBranchID() const;
		double SourceCurrent() const;
		
	private:
		void Initialize();
		bool ReadProperties(const EZ::List<EZ::String*>* line_tokens);
		double coefficient;
		
	protected:
		unsigned int source_branch_id;
	};
	
	class CCVolSource : public CCSource
	{
	public:
		CCVolSource();
		CCVolSource(const CCVolSource& source);
		~CCVolSource();
		CCVolSource& operator=(const CCVolSource& source);
		void Reset();
		void Print() const;
		
	private:
		void Initialize();
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
	};
	
	class CCCurrSource : public CCSource
	{
	public:
		CCCurrSource();
		CCCurrSource(const CCCurrSource& source);
		~CCCurrSource();
		CCCurrSource& operator=(const CCCurrSource& source);
		void Reset();
		void Print() const;
		
	private:
		void Initialize();
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
	};
	
	class ShortBranch : public Branch
	{
	public:
		ShortBranch();
		ShortBranch(const ShortBranch& branch);
		~ShortBranch();
		ShortBranch& operator=(const ShortBranch& branch);
		void Reset();
		void Print() const;
		
	private:
		void Initialize();
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
	};
	
	class OpenBranch : public Branch
	{
	public:
		OpenBranch();
		OpenBranch(const OpenBranch& branch);
		~OpenBranch();
		OpenBranch& operator=(const OpenBranch& branch);
		void Reset();
		void Print() const;
		
	private:
		void Initialize();
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
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
		void Print() const;
		
	private:
		void Initialize();
		bool ReadProperties(const EZ::List<EZ::String*>* line_tokens);
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
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
		void Print() const;
		
	private:
		void Initialize();
		bool ReadProperties(const EZ::List<EZ::String*>* line_tokens);
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
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
		void Print() const;
		
	private:
		void Initialize();
		bool ReadProperties(const EZ::List<EZ::String*>* line_tokens);
		void BranchEquation(EZ::Math::Matrix& f) const;
		void BranchGradients(EZ::Math::Matrix& A) const;
		double inductance;
	};
}

#endif

