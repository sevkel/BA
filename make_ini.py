import configparser
import sys
import os

#define config
config = configparser.ConfigParser()
config['SYSTEM'] = {'molecule': sys.argv[3]}
config['DFT'] = {'functional': sys.argv[1],
                 'basis': sys.argv[2]
                 }

#define output folder
molecule = config['SYSTEM']['molecule']
output_folder = './testsystems/' + molecule +'/'+ sys.argv[1] + '_' + sys.argv[2]

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#save file
configfile_path = os.path.join(output_folder, sys.argv[1] + '_' + sys.argv[2] + '.ini')
with open(configfile_path, 'w') as cfg:
    config.write(cfg)

