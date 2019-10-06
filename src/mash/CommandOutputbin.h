//output the binfile in text formatting.

#ifndef INCLUDED_CommandOutputbin
#define INCLUDED_CommandOutputbin

#include "Command.h"

namespace mash{

class CommandOutputbin : public Command
{
public:
	CommandOutputbin();
	int run() const; // override?
	
private:
	

};

} // namespace mash



#endif
