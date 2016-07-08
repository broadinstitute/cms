## relaunches cosi command until a selection trajectory is created successully
#IS THERE A WAY TO SILENCE COSI OUTPUT?
## last updated: 07.07.16 	vitti@broadinstitute.org

import sys
import subprocess
output, cosibuild, params, maxAttempts = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

commandstring = "env COSI_NEWSIM=1 COSI_MAXATTEMPTS=" + str(maxAttempts) + " COSI_SAVE_TRAJ=" + output + " " + cosibuild + " -p " + params 

itWorked = False
nAttempts = 0

while itWorked == False:
	nAttempts +=1
	try:
		subprocess.check_output(commandstring.split())
	except:
		continue
	itWorked = True	

#print("found a trajectory in " + str(nAttempts) + " attempts.")