import os
import sys

if len(sys.argv) != 2:
    print("Please insert the folder path for the systems to iterate as argv-argument")

main_path = sys.argv[1]

for roots, dirs, files in os.walk(main_path):

    for file in files:
        if file.endswith('.xyz'):
            xyz = file
        elif file.endswith('.ini'):
            ini = file
        os.system('DFT.py ' + xyz + ' ' + ini)







