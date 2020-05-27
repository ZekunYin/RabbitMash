#ifndef INCLUDED_COMMANDDUMPDIST
#define INCLUDED_COMMANDDUMPDIST

#include "Command.h"

namespace mash{

class CommandDumpdist: public Command
{
public:
	CommandDumpdist();
	int run() const; // override?
	
private:
	

};

} // namespace mash



#endif
