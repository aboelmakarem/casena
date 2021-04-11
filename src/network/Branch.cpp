// Branch.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/13/2021

#include "Branch.h"
#include "Node.h"
#include "Network.h"
#include "string.h"
#include "stdlib.h"

namespace CASENA
{
	Branch::Branch(){Initialize();}
	Branch::Branch(const Branch& branch) : Component(branch) {*this = branch;}
	Branch::~Branch(){Reset();}
	Branch& Branch::operator=(const Branch& branch)
	{
		Component::operator=(branch);
		memcpy(currents,branch.currents,HistoryCount*sizeof(double));
		start_node_id = branch.start_node_id;
		end_node_id = branch.end_node_id;
		return *this;
	}
	void Branch::Reset()
	{
		Initialize();
		Component::Reset();
	}
	bool Branch::IsBranchComponent(const char* line)
	{
		if(strncmp(line,"ivs",3) == 0)		return true;
		if(strncmp(line,"ics",3) == 0)		return true;
		if(strncmp(line,"vcvs",4) == 0)		return true;
		if(strncmp(line,"vccs",4) == 0)		return true;
		if(strncmp(line,"ccvs",4) == 0)		return true;
		if(strncmp(line,"cccs",4) == 0)		return true;
		if(strncmp(line,"short",5) == 0)	return true;
		if(strncmp(line,"open",4) == 0)		return true;
		if(strncmp(line,"res",3) == 0)		return true;
		if(strncmp(line,"cap",3) == 0)		return true;
		if(strncmp(line,"ind",3) == 0)		return true;
		return false;
	}
	unsigned int Branch::ReadMaxNodeID(const char* line)
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
		for(EZ::ListItem<EZ::String*>* item = tokens.Start() ; item != 0 ; item = item->Next())
		{
			delete item->Data();
		}
		tokens.Reset();
		return (node1_id > node2_id ? node1_id:node2_id);
	}
	Component* Branch::Create(const char* line)
	{
		if(strncmp(line,"ivs",3) == 0)		return new IndVolSource;
		if(strncmp(line,"ics",3) == 0)		return new IndCurrSource;
		if(strncmp(line,"vcvs",4) == 0)		return new VCVolSource;
		if(strncmp(line,"vccs",4) == 0)		return new VCCurrSource;
		if(strncmp(line,"ccvs",4) == 0)		return new CCVolSource;
		if(strncmp(line,"cccs",4) == 0)		return new CCCurrSource;
		if(strncmp(line,"short",5) == 0)	return new ShortBranch;
		if(strncmp(line,"open",4) == 0)		return new OpenBranch;
		if(strncmp(line,"res",3) == 0)		return new Resistor;
		if(strncmp(line,"cap",3) == 0)		return new Capacitor;
		if(strncmp(line,"ind",3) == 0)		return new Inductor;
		return 0;
	}
	bool Branch::Read(const char* line)
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
		// read node IDs
		token = tokens.Start();
		start_node_id = atoi((*(token->Data()))());
		delete token->Data();
		tokens.PopFront();
		token = tokens.Start();
		end_node_id = atoi((*(token->Data()))());
		delete token->Data();
		tokens.PopFront();
		// read the rest of the string and pass it to properties reader
		bool read_properties = ReadProperties(&tokens);
		if(!read_properties)
		{
			printf("error: incorrect definition for branch %s\n",Name()());
		}
		for(EZ::ListItem<EZ::String*>* item = tokens.Start() ; item != 0 ; item = item->Next())
		{
			delete item->Data();
		}
		tokens.Reset();
		return read_properties;
	}
	double Branch::Current() const{return currents[0];}
	double Branch::PreviousCurrent() const{return currents[1];}
	void Branch::Current(const double& value)
	{
		for(unsigned int i = HistoryCount - 1 ; i > 0 ; i--)
		{
			currents[i] = currents[i - 1];
		}
		currents[0] = value;
	}
	void Branch::StartNodeID(const unsigned int& target_id){start_node_id = target_id;}
	unsigned int Branch::StartNodeID() const{return start_node_id;}
	void Branch::EndNodeID(const unsigned int& target_id){end_node_id = target_id;}
	unsigned int Branch::EndNodeID() const{return end_node_id;}
	void Branch::Equation(EZ::Math::Matrix& f) const
	{
		// write the end nodes current equations
		f(start_node_id,0,f(start_node_id,0) + currents[0]);
		f(end_node_id,0,f(end_node_id,0) - currents[0]);
		BranchEquation(f);
	}
	void Branch::Gradients(EZ::Math::Matrix& A) const
	{
		// write the end nodes current equations gradients
		unsigned int branch_id = ID();
		A(start_node_id,branch_id,1.0);
		A(end_node_id,branch_id,-1.0);
		BranchGradients(A);
	}
	void Branch::Initialize()
	{
		memset(currents,0,HistoryCount*sizeof(double));
		start_node_id = 0;
		end_node_id = 0;
	}
	bool Branch::ReadProperties(const EZ::List<EZ::String*>* line_tokens){return (line_tokens->Size() == 0);}
		
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
	bool IndSource::ReadProperties(const EZ::List<EZ::String*>* line_tokens)
	{
		if(line_tokens->Size() != 1)		return false;
		form = atof((*(line_tokens->Front()))());
		return true;
	}
	
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
	void IndVolSource::Initialize(){}
	void IndVolSource::BranchEquation(EZ::Math::Matrix& f) const
	{
		// voltage source equation is v_end - v_start - v_source = 0
		Network* network = Network::GetNetwork();
		f(ID(),0,network->GetNode(end_node_id)->Voltage() - network->GetNode(start_node_id)->Voltage() - Form());
	}
	void IndVolSource::BranchGradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		A(branch_id,start_node_id,-1.0);
		A(branch_id,end_node_id,1.0);
	}
	void IndVolSource::Print() const
	{
		printf("voltage source %s: %u --> %u : voltage = %e\n",Name()(),start_node_id,end_node_id,Form());
	}
	
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
	void IndCurrSource::Initialize(){}
	void IndCurrSource::BranchEquation(EZ::Math::Matrix& f) const
	{
		// current source equation is i_branch - i_source = 0
		f(ID(),0,Current() - Form());
	}
	void IndCurrSource::BranchGradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		A(branch_id,branch_id,1.0);
	}
	void IndCurrSource::Print() const
	{
		printf("current source %s: %u --> %u : current = %e\n",Name()(),start_node_id,end_node_id,Form());
	}
	
	VCSource::VCSource(){Initialize();}
	VCSource::VCSource(const VCSource& source) : Branch(source){*this = source;}
	VCSource::~VCSource(){Reset();}
	VCSource& VCSource::operator=(const VCSource& source)
	{
		Branch::operator=(source);
		coefficient = source.coefficient;
		source_start_node_id = source.source_start_node_id;
		source_end_node_id = source.source_end_node_id;
		return *this;
	}
	void VCSource::Reset()
	{
		Initialize();
		Branch::Reset();
	}
	void VCSource::Coefficient(const double& value){coefficient = value;}
	double VCSource::Coefficient() const{return coefficient;}
	void VCSource::SourceStartNodeID(const unsigned int& id){source_start_node_id = id;}
	unsigned int VCSource::SourceStartNodeID() const{return source_start_node_id;}
	void VCSource::SourceEndNodeID(const unsigned int& id){source_end_node_id = id;}
	unsigned int VCSource::SourceEndNodeID() const{return source_end_node_id;}
	double VCSource::SourceVoltage() const
	{
		Network* network = Network::GetNetwork();
		double source_end_voltage = network->GetNode(source_end_node_id)->Voltage();
		double source_start_voltage = network->GetNode(source_start_node_id)->Voltage();
		return (source_end_voltage - source_start_voltage);
	}
	void VCSource::Initialize()
	{
		coefficient = 0.0;
		source_start_node_id = 0;
		source_end_node_id = 0;
	}
	bool VCSource::ReadProperties(const EZ::List<EZ::String*>* line_tokens)
	{
		if(line_tokens->Size() != 3)		return false;
		EZ::ListItem<EZ::String*>* token = line_tokens->Start();
		source_start_node_id = atoi((*(token->Data()))());
		token = token->Next();
		source_end_node_id = atoi((*(token->Data()))());
		token = token->Next();
		coefficient = atof((*(token->Data()))());
		return true;
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
	void VCVolSource::Initialize(){}
	void VCVolSource::BranchEquation(EZ::Math::Matrix& f) const
	{
		// voltage controlled voltage source equation is
		// v_end - v_start - coeff*(v_source_end - v_source_start) = 0
		Network* network = Network::GetNetwork();
		f(ID(),0,network->GetNode(end_node_id)->Voltage() - network->GetNode(start_node_id)->Voltage() - Coefficient()*SourceVoltage());
	}
	void VCVolSource::BranchGradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		double source_coefficient = Coefficient();
		A(branch_id,start_node_id,-1.0);
		A(branch_id,end_node_id,1.0);
		A(branch_id,source_start_node_id,source_coefficient);
		A(branch_id,source_end_node_id,-source_coefficient);
	}
	void VCVolSource::Print() const
	{
		printf("voltage controlled voltage source %s: %u --> %u : source voltage %u --> %u : coefficient = %e\n",
		Name()(),start_node_id,end_node_id,source_start_node_id,source_end_node_id,
		Coefficient());
	}
	
	VCCurrSource::VCCurrSource(){Initialize();}
	VCCurrSource::VCCurrSource(const VCCurrSource& source) : VCSource (source) {*this = source;}
	VCCurrSource::~VCCurrSource(){Reset();}
	VCCurrSource& VCCurrSource::operator=(const VCCurrSource& source)
	{
		VCSource::operator=(source);
		return *this;
	}
	void VCCurrSource::Reset(){VCSource::Reset();}
	void VCCurrSource::Initialize(){}
	void VCCurrSource::BranchEquation(EZ::Math::Matrix& f) const
	{
		// voltage controlled current source equation is
		// i_branch - coeff*(v_source_end - v_source_start) = 0
		f(ID(),0,Current() - Coefficient()*SourceVoltage());
	}
	void VCCurrSource::BranchGradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		double source_coefficient = Coefficient();
		A(branch_id,branch_id,1.0);
		A(branch_id,source_start_node_id,source_coefficient);
		A(branch_id,source_end_node_id,-source_coefficient);
	}
	void VCCurrSource::Print() const
	{
		printf("voltage controlled current source %s: %u --> %u : source voltage %u --> %u : coefficient = %e\n",
		Name()(),start_node_id,end_node_id,source_start_node_id,source_end_node_id,
		Coefficient());
	}
	
	CCSource::CCSource(){Initialize();}
	CCSource::CCSource(const CCSource& source) : Branch(source){*this = source;}
	CCSource::~CCSource(){Reset();}
	CCSource& CCSource::operator=(const CCSource& source)
	{
		Branch::operator=(source);
		coefficient = source.coefficient;
		source_branch_id = source.source_branch_id;
		return *this;
	}
	void CCSource::Reset()
	{
		Initialize();
		Branch::Reset();
	}
	void CCSource::Coefficient(const double& value){coefficient = value;}
	double CCSource::Coefficient() const{return coefficient;}
	void CCSource::SourceBranchID(const unsigned int& id){source_branch_id = id;}
	unsigned int CCSource::SourceBranchID() const{return source_branch_id;}
	double CCSource::SourceCurrent() const{return Network::GetNetwork()->GetBranch(source_branch_id)->Current();}
	void CCSource::Initialize()
	{
		coefficient = 0.0;
		source_branch_id = 0;
	}
	bool CCSource::ReadProperties(const EZ::List<EZ::String*>* line_tokens)
	{
		if(line_tokens->Size() != 2)		return false;
		EZ::ListItem<EZ::String*>* token = line_tokens->Start();
		source_branch_id = atoi((*(token->Data()))());
		token = token->Next();
		coefficient = atof((*(token->Data()))());
		return true;
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
	void CCVolSource::Initialize(){}
	void CCVolSource::BranchEquation(EZ::Math::Matrix& f) const
	{
		// current controlled voltage source equation is
		// v_end - v_start - coeff*i_source = 0
		Network* network = Network::GetNetwork();
		f(ID(),0,network->GetNode(end_node_id)->Voltage() - network->GetNode(start_node_id)->Voltage() - Coefficient()*SourceCurrent());
	}
	void CCVolSource::BranchGradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		double source_coefficient = Coefficient();
		A(branch_id,start_node_id,-1.0);
		A(branch_id,end_node_id,1.0);
		A(branch_id,source_branch_id,-source_coefficient);
	}
	void CCVolSource::Print() const
	{
		printf("current controlled voltage source %s: %u --> %u\n",Name()(),start_node_id,end_node_id);
	}
	
	CCCurrSource::CCCurrSource(){Initialize();}
	CCCurrSource::CCCurrSource(const CCCurrSource& source): CCSource(source){*this = source;}
	CCCurrSource::~CCCurrSource(){Reset();}
	CCCurrSource& CCCurrSource::operator=(const CCCurrSource& source)
	{
		CCSource::operator=(source);
		return *this;
	}
	void CCCurrSource::Reset(){CCSource::Reset();}
	void CCCurrSource::Initialize(){}
	void CCCurrSource::BranchEquation(EZ::Math::Matrix& f) const
	{
		// current controlled current source equation is
		// i_branch - coeff*i_source = 0
		f(ID(),0,Current() - Coefficient()*SourceCurrent());
	}
	void CCCurrSource::BranchGradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		A(branch_id,branch_id,1.0);
		A(branch_id,source_branch_id,-Coefficient());
	}
	void CCCurrSource::Print() const
	{
		printf("current controlled current source %s: %u --> %u\n",Name()(),start_node_id,end_node_id);
	}
	
	ShortBranch::ShortBranch(){Initialize();}
	ShortBranch::ShortBranch(const ShortBranch& branch) : Branch(branch){*this = branch;}
	ShortBranch::~ShortBranch(){Reset();}
	ShortBranch& ShortBranch::operator=(const ShortBranch& branch)
	{
		Branch::operator=(branch);
		return *this;
	}
	void ShortBranch::Reset()
	{
		Initialize();
		Branch::Reset();
	}
	void ShortBranch::Initialize(){}
	void ShortBranch::BranchEquation(EZ::Math::Matrix& f) const
	{
		// short circuit equation is v_end - v_start = 0
		Network* network = Network::GetNetwork();
		f(ID(),0,network->GetNode(end_node_id)->Voltage() - network->GetNode(start_node_id)->Voltage());
	}
	void ShortBranch::BranchGradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		A(branch_id,start_node_id,-1.0);
		A(branch_id,end_node_id,1.0);
	}
	void ShortBranch::Print() const
	{
		printf("short branch %s: %u --> %u\n",Name()(),start_node_id,end_node_id);
	}
	
	OpenBranch::OpenBranch(){Initialize();}
	OpenBranch::OpenBranch(const OpenBranch& branch) : Branch(branch){*this = branch;}
	OpenBranch::~OpenBranch(){Reset();}
	OpenBranch& OpenBranch::operator=(const OpenBranch& branch)
	{
		Branch::operator=(branch);
		return *this;
	}
	void OpenBranch::Reset()
	{
		Initialize();
		Branch::Reset();
	}
	void OpenBranch::Initialize(){}
	void OpenBranch::BranchEquation(EZ::Math::Matrix& f) const
	{
		// open circuit equation is i_branch = 0
		f(ID(),0,Current());
	}
	void OpenBranch::BranchGradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		A(branch_id,branch_id,1.0);
	}
	void OpenBranch::Print() const
	{
		printf("open branch %s: %u --> %u\n",Name()(),start_node_id,end_node_id);
	}
	
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
	void Resistor::Initialize(){resistance = 0.0;}
	bool Resistor::ReadProperties(const EZ::List<EZ::String*>* line_tokens)
	{
		if(line_tokens->Size() != 1)		return false;
		resistance = atof((*(line_tokens->Front()))());
		return true;
	}
	void Resistor::BranchEquation(EZ::Math::Matrix& f) const
	{
		// resistor equation is
		// v_end - v_start - resistance*i_branch = 0
		Network* network = Network::GetNetwork();
		f(ID(),0,network->GetNode(end_node_id)->Voltage() - network->GetNode(start_node_id)->Voltage() - resistance*Current());
	}
	void Resistor::BranchGradients(EZ::Math::Matrix& A) const
	{
		unsigned int branch_id = ID();
		A(branch_id,start_node_id,-1.0);
		A(branch_id,end_node_id,1.0);
		A(branch_id,branch_id,-resistance);
	}
	void Resistor::Print() const
	{
		printf("resistor %s: %u --> %u : resistance = %e\n",Name()(),start_node_id,end_node_id,resistance);
	}
	
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
	void Capacitor::Initialize(){capacitance = 0.0;}
	bool Capacitor::ReadProperties(const EZ::List<EZ::String*>* line_tokens)
	{
		if(line_tokens->Size() != 1)		return false;
		capacitance = atof((*(line_tokens->Front()))());
		return true;
	}
	void Capacitor::BranchEquation(EZ::Math::Matrix& f) const
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
	void Capacitor::BranchGradients(EZ::Math::Matrix& A) const
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
	void Capacitor::Print() const
	{
		printf("capacitor %s: %u --> %u : capacitance = %e\n",Name()(),start_node_id,end_node_id,capacitance);
	}
	
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
	void Inductor::Initialize(){inductance = 0.0;}
	bool Inductor::ReadProperties(const EZ::List<EZ::String*>* line_tokens)
	{
		if(line_tokens->Size() != 1)		return false;
		inductance = atof((*(line_tokens->Front()))());
		return true;
	}
	void Inductor::BranchEquation(EZ::Math::Matrix& f) const
	{
		if(Network::SteadyState())
		{
			// for steady state analysis, inductors are short circuits
			// inductor equation is v_end - v_start = 0
			Network* network = Network::GetNetwork();
			f(ID(),0,network->GetNode(end_node_id)->Voltage() - network->GetNode(start_node_id)->Voltage());
		}
		else
		{
			
		}
	}
	void Inductor::BranchGradients(EZ::Math::Matrix& A) const
	{
		if(Network::SteadyState())
		{
			unsigned int branch_id = ID();
			A(branch_id,start_node_id,-1.0);
			A(branch_id,end_node_id,1.0);
		}
		else
		{
			
		}
	}
	void Inductor::Print() const
	{
		printf("inductor %s: %u --> %u : inductance = %e\n",Name()(),start_node_id,end_node_id,inductance);
	}
}



