import os

def run(command_file, out_file):
    command = "charmm < {} > {}".format(command_file, out_file)
    print('Running ', command)
    os.system(command)

