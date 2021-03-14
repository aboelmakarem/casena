// Branch.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#include "Branch.h"
#include "Network.h"

namespace CASENA
{
	Branch::Branch(){Initialize();}
	Branch::Branch(const Branch& branch) : Component(branch) {*this = branch;}
	Branch::~Branch(){Reset();}
	Branch& Branch::operator=(const Branch& branch)
	{
		Component::operator=(branch);
		for(unsigned int i = 0 ; i < HistoryCount ; i++)
		{
			currents[i] = branch.currents[i];
		}
		start = branch.start;
		end = branch.end;
		return *this;
	}
	void Branch::Reset()
	{
		Initialize();
		Component::Reset();
	}
	double Branch::Current() const{return currents[0];}
	double Branch::PreviousCurrent() const{return currents[1];}
	void Branch::Current(const double& value)
	{
		for(unsigned int i = HistoryCount ; i > 0 ; i--)
		{
			currents[i] = currents[i - 1];
		}
		currents[0] = value;
	}
	void Branch::StartNode(const Node* node){start = node;}
	const Node* Branch::StartNode() const{return start;}
	void Branch::EndNode(const Node* node){end = node;}
	const Node* Branch::EndNode() const{return end;}
	void Branch::Initialize()
	{
		for(unsigned int i = 0 ; i < HistoryCount ; i++)
		{
			currents[i] = 0.0;
		}
		start = 0;
		end = 0;
	}
	
	IndSource::IndSource(){Initialize();}
	IndSource::IndSource(const IndSource& source) : Branch(source){*this = source;}
	IndSource::~IndSource(){Reset();}
	IndSource& IndSource::operator=(const IndSource& source)
	{
		Branch::operator=(source);
		form = source.form;
		return *this;
	}
	void IndSource::Reset()
	{
		Initialize();
		Branch::Reset();
	}
	void IndSource::Form(const double& value){form = value;}
	double IndSource::Form() const{return form;}
	void IndSource::Initialize(){form = 0.0;}
	
	IndVolSource::IndVolSource(){Initialize();}
	IndVolSource::IndVolSource(const IndVolSource& source) : IndSource(source){*this = source;}
	IndVolSource::~IndVolSource(){Reset();}
	IndVolSource& IndVolSource::operator=(const IndVolSource& source)
	{
		IndSource::operator=(source);
		return *this;
	}
	void IndVolSource::Reset()
	{
		Initialize();
		IndSource::Reset();
	}
	void IndVolSource::Equation(EZ::Math::Matrix& f) const
	{
		// voltage source equation is v_end - v_start - v_source = 0
		f(ID(),0,EndNode()->Voltage() - StartNode()->Voltage() - Form());
	}
	void IndVolSource::Gradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		A(branch_id,StartNode()->ID(),-1.0);
		A(branch_id,EndNode()->ID(),1.0);
	}
	void IndVolSource::Initialize(){}
	
	IndCurrSource::IndCurrSource(){Initialize();}
	IndCurrSource::IndCurrSource(const IndCurrSource& source) : IndSource(source){*this = source;}
	IndCurrSource::~IndCurrSource(){Reset();}
	IndCurrSource& IndCurrSource::operator=(const IndCurrSource& source)
	{
		IndSource::operator=(source);
		return *this;
	}
	void IndCurrSource::Reset()
	{
		Initialize();
		IndSource::Reset();
	}
	void IndCurrSource::Equation(EZ::Math::Matrix& f) const
	{
		// current source equation is i_branch - i_source = 0
		f(ID(),0,Current() - Form());
	}
	void IndCurrSource::Gradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		A(branch_id,branch_id,1.0);
	}
	void IndCurrSource::Initialize(){}
	
	VCSource::VCSource(){Initialize();}
	VCSource::VCSource(const VCSource& source) : Branch(source){*this = source;}
	VCSource::~VCSource(){Reset();}
	VCSource& VCSource::operator=(const VCSource& source)
	{
		Branch::operator=(source);
		coefficient = source.coefficient;
		source_start = source.source_start;
		source_end = source.source_end;
		return *this;
	}
	void VCSource::Reset()
	{
		Initialize();
		Branch::Reset();
	}
	void VCSource::Coefficient(const double& value){coefficient = value;}
	double VCSource::Coefficient() const{return coefficient;}
	void VCSource::SourceStartNode(const Node* node){source_start = node;}
	const Node* VCSource::SourceStartNode() const{return source_start;}
	void VCSource::SourceEndNode(const Node* node){source_end = node;}
	const Node* VCSource::SourceEndNode() const{return source_end;}
	double VCSource::SourceVoltage() const{return (source_end->Voltage() - source_start->Voltage());}
	void VCSource::Initialize()
	{
		coefficient = 0.0;
		source_start = 0;
		source_end = 0;
	}
	
	VCVolSource::VCVolSource(){Initialize();}
	VCVolSource::VCVolSource(const VCVolSource& source) : VCSource (source) {*this = source;}
	VCVolSource::~VCVolSource(){Reset();}
	VCVolSource& VCVolSource::operator=(const VCVolSource& source)
	{
		VCSource::operator=(source);
		return *this;
	}
	void VCVolSource::Reset(){VCSource::Reset();}
	void VCVolSource::Equation(EZ::Math::Matrix& f) const
	{
		// voltage controlled voltage source equation is
		// v_end - v_start - coeff*(v_source_end - v_source_start) = 0
		f(ID(),0,EndNode()->Voltage() - StartNode()->Voltage() - Coefficient()*SourceVoltage());
	}
	void VCVolSource::Gradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		double source_coefficient = Coefficient();
		A(branch_id,StartNode()->ID(),-1.0);
		A(branch_id,EndNode()->ID(),1.0);
		A(branch_id,SourceStartNode()->ID(),source_coefficient);
		A(branch_id,SourceEndNode()->ID(),-source_coefficient);
	}
	void VCVolSource::Initialize(){}
	
	VCCurrSource::VCCurrSource(){Initialize();}
	VCCurrSource::VCCurrSource(const VCCurrSource& source) : VCSource (source) {*this = source;}
	VCCurrSource::~VCCurrSource(){Reset();}
	VCCurrSource& VCCurrSource::operator=(const VCCurrSource& source)
	{
		VCSource::operator=(source);
		return *this;
	}
	void VCCurrSource::Reset(){VCSource::Reset();}
	void VCCurrSource::Equation(EZ::Math::Matrix& f) const
	{
		// voltage controlled current source equation is
		// i_branch - coeff*(v_source_end - v_source_start) = 0
		f(ID(),0,Current() - Coefficient()*SourceVoltage());
	}
	void VCCurrSource::Gradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		double source_coefficient = Coefficient();
		A(branch_id,branch_id,1.0);
		A(branch_id,SourceStartNode()->ID(),source_coefficient);
		A(branch_id,SourceEndNode()->ID(),-source_coefficient);
	}
	void VCCurrSource::Initialize(){}
	
	CCSource::CCSource(){Initialize();}
	CCSource::CCSource(const CCSource& source) : Branch(source){*this = source;}
	CCSource::~CCSource(){Reset();}
	CCSource& CCSource::operator=(const CCSource& source)
	{
		Branch::operator=(source);
		coefficient = source.coefficient;
		source_branch = source.source_branch;
		return *this;
	}
	void CCSource::Reset()
	{
		Initialize();
		Branch::Reset();
	}
	void CCSource::Coefficient(const double& value){coefficient = value;}
	double CCSource::Coefficient() const{return coefficient;}
	void CCSource::SourceBranch(const Branch* branch){source_branch = branch;}
	const Branch* CCSource::SourceBranch() const{return source_branch;}
	double CCSource::SourceCurrent() const{return source_branch->Current();}
	void CCSource::Initialize()
	{
		coefficient = 0.0;
		source_branch = 0;
	}
	
	CCVolSource::CCVolSource(){Initialize();}
	CCVolSource::CCVolSource(const CCVolSource& source): CCSource(source){*this = source;}
	CCVolSource::~CCVolSource(){Reset();}
	CCVolSource& CCVolSource::operator=(const CCVolSource& source)
	{
		CCSource::operator=(source);
		return *this;
	}
	void CCVolSource::Reset(){CCSource::Reset();}
	void CCVolSource::Equation(EZ::Math::Matrix& f) const
	{
		// current controlled voltage source equation is
		// v_end - v_start - coeff*i_source = 0
		f(ID(),0,EndNode()->Voltage() - StartNode()->Voltage() - Coefficient()*SourceCurrent());
	}
	void CCVolSource::Gradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		double source_coefficient = Coefficient();
		A(branch_id,StartNode()->ID(),-1.0);
		A(branch_id,EndNode()->ID(),1.0);
		A(branch_id,SourceBranch()->ID(),-source_coefficient);
	}
	void CCVolSource::Initialize(){}
	
	CCCurrSource::CCCurrSource(){Initialize();}
	CCCurrSource::CCCurrSource(const CCCurrSource& source): CCSource(source){*this = source;}
	CCCurrSource::~CCCurrSource(){Reset();}
	CCCurrSource& CCCurrSource::operator=(const CCCurrSource& source)
	{
		CCSource::operator=(source);
		return *this;
	}
	void CCCurrSource::Reset(){CCSource::Reset();}
	void CCCurrSource::Equation(EZ::Math::Matrix& f) const
	{
		// current controlled current source equation is
		// i_branch - coeff*i_source = 0
		f(ID(),0,Current() - Coefficient()*SourceCurrent());
	}
	void CCCurrSource::Gradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		A(branch_id,branch_id,1.0);
		A(branch_id,SourceBranch()->ID(),-Coefficient());
	}
	void CCCurrSource::Initialize(){}
	
	Resistor::Resistor(){Initialize();}
	Resistor::Resistor(const Resistor& resistor) : Branch(resistor){*this = resistor;}
	Resistor::~Resistor(){Reset();}
	Resistor& Resistor::operator=(const Resistor& resistor)
	{
		Branch::operator=(resistor);
		resistance = resistor.resistance;
		return *this;
	}
	void Resistor::Reset()
	{
		Initialize();
		Branch::Reset();
	}
	void Resistor::Resistance(const double& value){resistance = value;}
	double Resistor::Resistance() const{return resistance;}
	void Resistor::Equation(EZ::Math::Matrix& f) const
	{
		// resistor equation is
		// v_end - v_start - resistance*i_branch = 0
		f(ID(),0,EndNode()->Voltage() - StartNode()->Voltage() - resistance*Current());
	}
	void Resistor::Gradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		A(branch_id,StartNode()->ID(),-1.0);
		A(branch_id,EndNode()->ID(),1.0);
		A(branch_id,branch_id,-resistance);
	}
	void Resistor::Initialize(){resistance = 0.0;}
	
	Capacitor::Capacitor(){Initialize();}
	Capacitor::Capacitor(const Capacitor& capacitor) : Branch(capacitor){*this = capacitor;}
	Capacitor::~Capacitor(){Reset();}
	Capacitor& Capacitor::operator=(const Capacitor& capacitor)
	{
		Branch::operator=(capacitor);
		capacitance = capacitor.capacitance;
		return *this;
	}
	void Capacitor::Reset()
	{
		Initialize();
		Branch::Reset();
	}
	void Capacitor::Capacitance(const double& value){capacitance = value;}
	double Capacitor::Capacitance() const{return capacitance;}
	void Capacitor::Equation(EZ::Math::Matrix& f) const
	{
		if(Network::SteadyState())
		{
			// for steady state analysis, capacitors are open circuits
			// capacitor equation is i_branch = 0
			f(ID(),0,Current());
		}
		else
		{
			
		}
	}
	void Capacitor::Gradients(EZ::Math::Matrix& A) const
	{
		if(Network::SteadyState())
		{
			unsigned int branch_id = ID();
			A(branch_id,branch_id,1.0);
		}
		else
		{
			
		}
	}
	void Capacitor::Initialize(){capacitance = 0.0;}
	
	Inductor::Inductor(){Initialize();}
	Inductor::Inductor(const Inductor& inductor) : Branch(inductor){*this = inductor;}
	Inductor::~Inductor(){Reset();}
	Inductor& Inductor::operator=(const Inductor& inductor)
	{
		Branch::operator=(inductor);
		inductance = inductor.inductance;
		return *this;
	}
	void Inductor::Reset()
	{
		Initialize();
		Branch::Reset();
	}
	void Inductor::Inductance(const double& value){inductance = value;}
	double Inductor::Inductance() const{return inductance;}
	void Inductor::Equation(EZ::Math::Matrix& f) const
	{
		if(Network::SteadyState())
		{
			// for steady state analysis, inductors are short circuits
			// inductor equation is v_end - v_start = 0
			f(ID(),0,EndNode()->Voltage() - StartNode()->Voltage());
		}
		else
		{
			
		}
	}
	void Inductor::Gradients(EZ::Math::Matrix& A) const
	{
		if(Network::SteadyState())
		{
			unsigned int branch_id = ID();
			A(branch_id,StartNode()->ID(),-1.0);
			A(branch_id,EndNode()->ID(),1.0);
		}
		else
		{
			
		}
	}
	void Inductor::Initialize(){inductance = 0.0;}
}



