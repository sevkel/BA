import sys, os
import configparser


# check number of arguments as input
if len(sys.argv) < 2:
    print("Please insert the filename as argument.")
    sys.exit(1)

# Define sections and parameters
sections = {
    'functionals': ['func1', 'func2', 'func3'],
    'basis': ['basis1', 'basis2'],
    'system': ['molecule']
    # Add more sections and parameters if wanted
}

# create an instance of the configparser
config = configparser.ConfigParser()


# Add sections and parameters
for section, params in sections.items():
    config.add_section(section)
    for param in params:
        value = input(f"Please put in a value for {param} in {section}: ")
        config.set(section, param, value)

#define output folder
system = config['system']['molecule']
output_folder = './testsystems/' + system

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# save INI-file
configfile_path = os.path.join(output_folder, system + '.ini')
with open(configfile_path, 'w') as configfile:
    config.write(configfile)

print(f"The INI-file for '{system}' was created successfully.")