import sys, os
import configparser

# check number of arguments as input
if len(sys.argv) != 2:
    print("Please insert the system- or molecule name as an argument.")
    sys.exit(1)

# Define sections and parameters for functionals and basis sets
sections = {
    'functionals': [],
    'basis': [],
    # Add more sections and parameters if wanted
}

# create an instance of the configparser
config = configparser.ConfigParser()

# Add sections and parameters
for section, params in sections.items():
    config.add_section(section)
    if section == 'functionals':
        i = 1
        while True:
            param = 'functional' + str(i)
            value = input(f"Please put in a value for {param} (or '*' to stop): ")
            if value == '*':
                break  # Exit the loop if '*' is entered
            params.append(param)
            config.set(section, param, value)
            i += 1
    elif section == 'basis':
        i = 1
        while True:
            param = 'basis' + str(i)
            value = input(f"Please put in a value for {param} (or '*' to stop): ")
            if value == '*':
                break  # Exit the loop if '*' is entered
            params.append(param)
            config.set(section, param, value)
            i += 1
config['system'] = {'molecule': sys.argv[1]}

# Define output folder
system = config['system']['molecule']
output_folder = './testsystems/' + system

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Save INI-file
configfile_path = os.path.join(output_folder, system + '.ini')
with open(configfile_path, 'w') as configfile:
    config.write(configfile)

print(f"The INI-file for '{system}' was created successfully.")
