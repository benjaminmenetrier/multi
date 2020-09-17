#!/usr/bin/env python 

#imported packages:
import sys
import subprocess

################################################################################
def run_model(path_to_code,code_args):
	
	# """ This function allows the user to run a code as if it was done in
	#     command line in a shell, according to its arguments. This function is
	#     usefull if one wants to call several times the same code with different
	#     arguments.

	#     Parameters:
	#     path_to_code (str): command line to run the executable file of the code
	#     without its arguments.
	#     code_args (list): ordered list containing the arguments of the code
	    
	#     Returns:
	#     out_list (list): List containing all the output lines of the code.
	# """
	
        # Creates the command-line to run the code:
        arguments=''
        for par in code_args:
                arguments+=' {} '.format(par)
		#print(arguments)
        exec_command ='echo '+ arguments + ' | ' + path_to_code
	#print(exec_command)
        # Run the code:
        process = subprocess.Popen(exec_command,shell=True,stdout=subprocess.PIPE,
				   stderr=subprocess.STDOUT)
	# Store the outputs of the code:
        out_list=[]
        bad_params=[]
        while True:
                nextline=process.stdout.readline()
                if nextline== '' and process.poll() is not None:
                        break
                if (process.returncode != None):
                        sys.exit(-1)
                        break
                out_list.append(nextline)
        return out_list
################################################################################

# quick test:
n=128
no=3
ni=5
lmp_mode='ritz'
full_res='f'
new_seed='T'
sigma_obs=0.1
sigmabvar=0.
Lb=0.

code_args=[n,no,ni,lmp_mode,full_res,new_seed,sigma_obs,sigmabvar,Lb]
path_to_code='../run/multi '

a=run_model(path_to_code, code_args)

print(a)
