import os
import subprocess
openfile = open('remove.txt')
for line in openfile:
	badfile = line.strip('\n')
	if os.path.isfile(badfile):
		rmcmd = "rm " + badfile
		subprocess.check_output(rmcmd.split())
