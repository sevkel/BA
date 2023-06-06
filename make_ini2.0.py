import sys, os
import configparser


# Überprüfe, ob ausreichend Argumente übergeben wurden
if len(sys.argv) < 2:
    print("Bitte geben Sie den Dateinamen als Argument an.")
    sys.exit(1)

# Definiere die Sektionen und Parameter
sections = {
    'functionals': ['func1', 'func2', 'func3'],
    'basis': ['basis1', 'basis2'],
    'system': ['molecule']
    # Füge weitere Sektionen und Parameter hinzu
}

# Erstelle eine Instanz des configparser
config = configparser.ConfigParser()


# Füge Sektionen und Parameter hinzu
for section, params in sections.items():
    config.add_section(section)
    for param in params:
        value = input(f"Geben Sie den Wert für {param} in {section} ein: ")
        config.set(section, param, value)

#define output folder
system = config['system']['molecule']
output_folder = './testsystems/' + system

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# Speichere die INI-Datei
configfile_path = os.path.join(output_folder, system + '.ini')
with open(configfile_path, 'w') as configfile:
    config.write(configfile)

print(f"Die INI-Datei für '{system}' wurde erfolgreich erstellt.")