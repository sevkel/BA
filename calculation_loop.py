import os
import sys

# dir path
dirpath = sys.argv[1]

# iterate through folders
for dirname in os.listdir(dirpath):
    folder_path = os.path.join(dirpath, dirname)

    # check if it's a folder
    if os.path.isdir(folder_path):
        ini_files = []
        xyz_files = []

        # iterate through files in current folder
        for filename in os.listdir(folder_path):
            datei_path = os.path.join(folder_path, filename)

            # check file endings
            if filename.endswith('.ini'):
                ini_files.append(datei_path)
            elif filename.endswith('.xyz'):
                xyz_files.append(datei_path)

        # put files into DFT.py and run it
        for ini_datei in ini_files:
            for xyz_datei in xyz_files:
                command = f'sbatch /alcc/gpfs2/home/u/kellerse/jobsub.sh {xyz_datei} {ini_datei}'
                #print(command)
                os.system(command)







