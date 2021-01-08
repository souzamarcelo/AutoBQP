#!/usr/bin/python3
import os
import sys
import subprocess
import shutil

if len(sys.argv) < 5:
    print("\nUsage: ./target-runner.py <configuration_id> <instance_id> <seed> <instance_path_name> <list of parameters>\n")
    sys.exit(1)

# Get the parameters as command line arguments.
configuration_id = sys.argv[1]
instance_id = sys.argv[2]
seed = sys.argv[3]
instance = sys.argv[4]
conf_params = ' '.join(sys.argv[5:])
index = conf_params.index('--cAlg')
instance_params = conf_params[:index]
conf_params = conf_params[index:]
fixed_params = ''

# Create directory and file names
directory = './src-' + str(configuration_id) + '-' + str(instance_id) + '-' + str(seed)
algorithm = './alg-' + str(configuration_id) + '-' + str(instance_id) + '-' + str(seed) + '.cpp'

# Call grammar to generate algorithm according to configuration parameters
os.system('./grammar ' + conf_params + ' > ' + algorithm)

# Create directory with the code and copy algorithm
shutil.copytree('./src', directory)
os.unlink(directory + '/algorithm.cpp')
shutil.move(algorithm, directory + '/algorithm.cpp')

# Compile the generated algorithm into the solver
os.system('cd ' + directory + '; cmake . 1> /dev/null 2> /dev/null')
os.system('cd ' + directory + '; make runner 1> /dev/null 2> /dev/null')

# Create execution command
exe = directory + '/runner'
command = exe + ' ' + fixed_params + ' --ins ' + instance + ' ' + instance_params + ' --seed ' + seed

# Define the stdout and stderr files.
out_file = directory + '/c' + str(configuration_id) + '-' + str(instance_id) + str(seed) + '.stdout'
err_file = directory + '/c' + str(configuration_id) + '-' + str(instance_id) + str(seed) + '.stderr'

# Execute
outf = open(out_file, "w")
errf = open(err_file, "w")
return_code = subprocess.call(command, stdout = outf, stderr = errf, shell = True)
outf.close()
errf.close()

# Check return code
if return_code != 0:
    print('Command returned code ' + str(return_code) + '.')
    sys.exit(1)

# Check output file
if not os.path.isfile(out_file):
    print('Output file <' + out_file  + '> not found.')
    sys.exit(1)

# Get result and print it
result = open(out_file).read().split()[-1]
print(result)

# Clean files and exit
os.remove(out_file)
os.remove(err_file)
shutil.rmtree(directory)
sys.exit(0)


#!/usr/bin/python3
import os
import sys
import subprocess

if len(sys.argv) < 5:
    print("\nUsage: ./target-runner.py <configuration_id> <instance_id> <seed> <instance_path_name> <list of parameters>\n")
    sys.exit(1)

# Get the parameters as command line arguments.
configuration_id = sys.argv[1]
instance_id = sys.argv[2]
seed = sys.argv[3]
instance = sys.argv[4]
conf_params = ' '.join(sys.argv[5:])

# Create execution command
fixed_params = '-r 1 --quiet -t 10'
command = './acotsp ' + fixed_params + ' -i ' + instance + ' --seed ' + seed + ' ' + conf_params


# Define the stdout and stderr files.
out_file = './c' + str(configuration_id) + '-' + str(instance_id) + str(seed) + '.stdout'
err_file = './c' + str(configuration_id) + '-' + str(instance_id) + str(seed) + '.stderr'

# Execute
outf = open(out_file, "w")
errf = open(err_file, "w")
return_code = subprocess.call(command, stdout = outf, stderr = errf, shell = True)
outf.close()
errf.close()

# Check return code
if return_code != 0:
    print('Command returned code ' + str(return_code) + '.')
    sys.exit(1)

# Check output file
if not os.path.isfile(out_file):
    print('Output file <' + out_file  + '> not found.')
    sys.exit(1)
    os.remove(out_file)
    os.remove(err_file)
    sys.exit(0)

# Check if an error occured
if os.stat(out_file).st_size != 0:
    print('Inf')
        
# Get result and print it
result = open(out_file).read().split()[-1].split(';')[1]
print(result)

# Clean files and exit
os.remove(out_file)
os.remove(err_file)
sys.exit(0)
