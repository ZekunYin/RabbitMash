//output the binfile in text formatting.

#ifndef INCLUDED_CommandDumpTri
#define INCLUDED_CommandDumpTri

#include "Command.h"

namespace mash{

class CommandDumptri : public Command
{
public:
	CommandDumptri();
	int run() const; // override?
	
private:
	

};

} // namespace mash



#endif
