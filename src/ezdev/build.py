
import os
import os.path
import sys

sources = ["Random.cpp","BLAS.cpp","Matrix.cpp"]

compiler = "mpicxx"
include_arg = ""

output = "libezdev.a"

object_extension = ".o"

print_only = False

def NeedToCompile(sourcefile,objectfile):
	# check if object file exists
	if(os.path.isfile(objectfile)):
		# check its time signature
		sourcetime = os.path.getmtime(sourcefile)
		objecttime = os.path.getmtime(objectfile)
		# leave a tolerance of 10 seconds to avoid any inaccurate time reporting
		if(sourcetime < (objecttime - 10)):
			return False
	return True

def Execute(command):
	if(print_only):
		print(command)
	else:
		print(command)
		os.system(command)

def GetObjectName(sourcename):
	objectname = sourcename[0:sourcename.rfind('.')] + object_extension
	return objectname

def CompileSource(sourcefile,debug):
	if(NeedToCompile(sourcefile,GetObjectName(sourcefile))):
		if(debug):	compile_string = compiler + "  -Wall -Wextra -pedantic -ansi -g " + include_arg + " -c " + sourcefile
		else:		compile_string = compiler + "  -Wall -Wextra -pedantic -ansi -O2 " + include_arg + " -c " + sourcefile
		Execute(compile_string)
		return True
	return False

def CompileSources(debug):
	for source in sources:
		CompileSource(source,debug)

def Compile(debug):
	CompileSources(debug)

def Link():
	# delete any existing libraries before creating a new one
	remove_command = "rm -f " + output
	Execute(remove_command)
	link_command = "ar -cvq " + output + " "
	for source in sources:
		link_command = link_command + GetObjectName(source) + " "
	Execute(link_command)

clean = False
debug = False

for i in range(1,len(sys.argv)):
	if(sys.argv[i] == "clean"):
		clean = True
	elif(sys.argv[i] == "debug"):
		debug = True

if(clean):
	Execute("rm *.o")
	Execute("rm " + output)

Compile(debug)
Link()


