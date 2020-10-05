#!/usr/bin/python3
import os
import sys
import subprocess
import shutil
import argparse
import datetime

desc = '''
------------------------------------------------------------------------------------
|                                     AutoBQP                                      |
|  A black-box solver for binary problems based on the automatic configuration of  |
|     heuristic algorithms for the unconstrained bynary quadratic programming      |
|                                                                                  |
| Version: 1.0                                                                     |
| Copyright (C) 2020                                                               |
| Marcelo de Souza         <marcelo.desouza@udesc.br>                              |
| Marcus Ritt              <marcus.ritt@inf.ufrgs.br>                              |
|                                                                                  |
| This is free software, and you are welcome to redistribute it under certain      |
| conditions.  See the GNU General Public License for details. There is NO         |
| WARRANTY; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.      |
------------------------------------------------------------------------------------
'''

directory = ''

def create_execution(insfile, insdir, budget, parallel):
    global directory
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    directory = 'execution_' + now
    os.mkdir(directory)
    shutil.copy(insfile, directory + '/instances.txt')
    shutil.copytree(insdir, directory + '/instances')
    shutil.copy('./aac/scenario.txt', directory)
    with open(directory + '/scenario.txt', 'a') as file:
        file.write('maxExperiments = ' + str(budget) + '\n')
        file.write('parallel = ' + str(parallel) + '\n')
    shutil.copy('./aac/parameters.txt', directory)
    shutil.copytree('./src', directory + '/src')
    os.system('cd ' + directory + '/src; cmake . >> /dev/null')
    os.system('cd ' + directory + '/src; make grammar >> /dev/null')
    shutil.copy(directory + '/src/grammar', directory)
    shutil.copy('./aac/target-runner.py', directory)


def irace():
    command = 'cd ' + directory + '; irace 1> results-irace.txt 2> error-irace.txt'
    subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT).stdout.read()
    if os.stat(directory + '/error-irace.txt').st_size != 0:
        return str(open(directory + '/error-irace.txt', 'r').read())
    return ''


def read_algorithm():
    result_file = open(directory + '/results-irace.txt', 'r')
    read = False
    for line in result_file:
        algorithm = line
        if read: break
        read = ('# Best configurations as commandlines' in line)
    result_file.close()
    algorithm = algorithm.replace('\n', '')
    while '  ' in algorithm: algorithm = algorithm.replace('  ', ' ')
    return algorithm[algorithm.index(' ') + 1:]


def organize_algorithm():
    algorithm = read_algorithm()
    algorithm_dir = directory.replace('execution', 'algorithm')
    shutil.copytree('./src', algorithm_dir)
    os.system('cd ./' + algorithm_dir + '; mkdir build >> /dev/null; cd build; cmake .. >> /dev/null; make grammar >> /dev/null; cp grammar .. >> /dev/null;')
    os.system('cd ./' + algorithm_dir + '; ./grammar ' + algorithm + ' > algorithm.cpp; cd build; cmake .. >> /dev/null; make runner >> /dev/null; cp runner .. >> /dev/null;')
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-v', '--version', help = 'show description and exit', action = 'store_true')
    optional.add_argument('--insfile', help = 'file with each instance (name) and specific arguments (if any)', metavar = '<file>', type = str, default = './aac/instances.txt')
    optional.add_argument('--insdir', help = 'directory with problem instances', metavar = '<folder>', type = str, default = './aac/instances')
    optional.add_argument('--budget', help = 'maximum number of executions for the automatic design process', metavar = '<budget>', type = int, default = 2000)
    optional.add_argument('--parallel', help = 'number of cores to be used for parallelization', metavar = '<cores>', type = int, default = 1)
    args, other = parser.parse_known_args()
    if args.version: print(desc); exit()

    print(desc)

    insfile = os.path.abspath(args.insfile)
    insdir = os.path.abspath(args.insdir)
    print('> Using instances listed in this file: ' + insfile)
    print('> Reading instance files from this directory: ' + insdir)
    print('> Constructing the algorithm with budget: ' + str(args.budget))
    print('> Using the following number of cores: ' + str(args.parallel))
    print()

    print('Creating execution directory', end = '', flush = True)
    create_execution(insfile, insdir, args.budget, args.parallel)
    print(' .................... DONE')

    print('Executing the automatic design process', end = '', flush = True)
    result = irace()
    if result == '':
        print(' .......... DONE')
    else:
        print(' .......... ERROR')
        print('\nirace returned the following message:\n' + result)
        print('Aborting.')
        print('------------------------------------------------------------------------------------')
        exit(1)

    print('Organizing the resulting algorithm', end = '', flush = True)
    organize_algorithm()
    print(' .............. DONE')

    print()
    print('Algorithm is in the following directory: ' + directory.replace('execution', 'algorithm'))

    print()
    print('------------------------------------------------------------------------------------')