// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandBounds.h"
#include "CommandList.h"
#include "CommandSketch.h"
#include "CommandFind.h"
#include "CommandDistance.h"
#include "CommandScreen.h"
#include "CommandTriangle.h"
#include "CommandContain.h"
#include "CommandInfo.h"
#include "CommandPaste.h"

#include "CommandDumptri.h"
#include "CommandDumpdist.h"

#include <sched.h>
#include <unistd.h>

int main(int argc, const char ** argv)
{
    mash::CommandList commandList("mash");
    
    commandList.addCommand(new mash::CommandSketch());
    //commandList.addCommand(new CommandFind());
    commandList.addCommand(new mash::CommandDistance());
    commandList.addCommand(new mash::CommandScreen());
    commandList.addCommand(new mash::CommandTriangle());
#ifdef COMMAND_WITHIN
    commandList.addCommand(new mash::CommandContain());
#endif
#ifdef COMMAND_FIND
	commandList.addCommand(new mash::CommandFind());
#endif
    commandList.addCommand(new mash::CommandInfo());
    commandList.addCommand(new mash::CommandPaste());
    commandList.addCommand(new mash::CommandBounds());
	commandList.addCommand(new mash::CommandDumptri());
	commandList.addCommand(new mash::CommandDumpdist());

	//cpu_set_t cpus;
	//CPU_ZERO(&cpus);
	//CPU_SET(0, &cpus);
	//sched_setaffinity(0, sizeof(cpus), &cpus);

    return commandList.run(argc, argv);
}
