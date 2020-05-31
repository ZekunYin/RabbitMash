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
int64_t getFileSize( const char * fileName);

} // namespace mash



#endif
