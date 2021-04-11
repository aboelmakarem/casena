// Node.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 03/14/2021

#include "Component.h"

namespace CASENA
{
	Component::Component(){Initialize();}
	Component::Component(const Component& component){*this = component;}
	Component::~Component(){Reset();}
	Component& Component::operator=(const Component& component)
	{
		id = component.id;
		name = component.name;
		return *this;
	}
	void Component::Reset(){Initialize();}
	bool Component::Read(const char* line)
	{
		// a dummy check to avoid unused parameter compilation warning
		if(line == 0)		return true;
		return true;
	}
	void Component::ID(const unsigned int& value){id = value;}
	unsigned int Component::ID() const{return id;}
	void Component::Name(const EZ::String& target_name){name = target_name;}
	const EZ::String& Component::Name() const{return name;}
	unsigned int Component::ClaimIDs(const unsigned int& start_id)
	{
		// claim only 1 ID
		id = start_id;
		return id + 1;
	}
	void Component::Initialize()
	{
		id = 0;
		name.Reset();
	}
}

