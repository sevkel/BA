import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import os
import numpy as np
import scienceplots

prop = fm.FontProperties(fname='/alcc/gpfs2/home/u/kellerse/fonts/fira-sans/FiraSans-Regular.ttf')
font_name = prop.get_name()
plt.style.use(['science','notebook','no-latex'])
"""plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=10)"""

def get_ref(molec):
    with open("/alcc/gpfs2/home/u/kellerse/results/reference_W4-11.txt", 'r') as ref:
        lines = ref.readlines()
        for line in lines:
            if molec == line.split()[1]:
                refval = line.split()[2]
                ref.close()
                return float(refval)

def get_func_data(func, path):
    func_data = {}
    for molecdir in os.listdir(path):
        molec_dir_path = os.path.join(path, molecdir)
        molec_name = molecdir
        ref = get_ref(molec_name)
        for func_dir in os.listdir(molec_dir_path):
            if func_dir.startswith(func):
                func_path = os.path.join(molec_dir_path, func_dir)
                for file in os.listdir(func_path):
                    if file == 'energies.txt':
                        file_path = os.path.join(func_path, file)
                        with open(file_path, 'r') as f:
                            lines = f.readlines()
                            for line in lines:
                                if line.startswith('Reaction energy:'):
                                    e_react = float(line.split(':')[1].split()[0])
                                    func_data[molec_name] = e_react - ref
                                    break
                            f.close()
                break
    return func_data

def find_outliers(data):
    q1 = sorted(data)[int(len(data) * 0.25)]
    q3 = sorted(data)[int(len(data) * 0.75)]
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return [(key, val) for key, val in data_dict.items() if float(val) < lower_bound or float(val) > upper_bound]

maindir = '/alcc/gpfs2/home/u/kellerse/results/molecs_W4-11-v3'

funcs = ['DM21', 'DM21m', 'DM21mc', 'DM21mu', 'GGA_XC_PBE', 'MGGA_XC_SCAN', 'lda', 'pbe', 'pbe0',
                         'r2scan', 'tpss']

DM21 = get_func_data('DM21', maindir)
print(DM21)

for func in funcs:
    data_dict = get_func_data(f'{func}', maindir)

    data_values = [float(val) for val in data_dict.values()]
    print(func)
    print(np.mean(data_values))
    outliers = find_outliers(data_values)
    print(outliers)
    # Erstelle den Boxplot
    plt.figure(figsize=(5, 5))
    plt.grid()
    plt.xticks(fontproperties=prop)
    plt.yticks(fontproperties=prop)
    bp = plt.boxplot(data_values, patch_artist=True, positions=[1], medianprops={'color': 'black'},
                     showmeans=True, meanline=True, sym='')
    plt.legend([bp['medians'][0], bp['means'][0]], ['Median', 'Mittelwert'], prop=prop)

    # Färbe den Boxplot
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    # Zeichne die Outlier als Punkte und beschrifte sie mit ihren Keys neben den Reaktions-Indizes
    flag = True
    sorted_outliers = sorted(outliers, key=lambda item: item[1])
    for outlier_key, outlier_value in sorted_outliers:
        outlier_x = 1
        outlier_y = float(outlier_value)
        print((outlier_key, outlier_value))
        plt.scatter(outlier_x, outlier_y, color='red', label='Outlier')
        if flag == True:
            plt.annotate(outlier_key.upper(), (outlier_x, outlier_y), textcoords='offset points', xytext=(15, 0.5), ha='left',
                         va='center', fontproperties=prop, size=10)
            flag = False

        else:
            plt.annotate(outlier_key.upper(), (outlier_x, outlier_y), textcoords='offset points', xytext=(-45, 0.5), ha='left',
                         va='center', fontproperties=prop, size=10)
            flag = True


    # Berechne den Mittelwert
    mean_value = np.mean(data_values)


    # Setze die x-Achse mit deinem gewünschten String
    if func in ['lda', 'pbe', 'pbe0', 'r2scan', 'tpss']:
        func = func.upper()
    plt.xticks([1], [f'{func}'], fontproperties=prop, size=14)
    plt.yticks(fontproperties=prop, size=10)
    plt.ylabel('Reaktionsenergie (Ergebnis - Referenz) (kcal/mol)', fontproperties=prop, size=14)
    #plt.title('')
    plt.tight_layout()
    plt.savefig(f'/alcc/gpfs2/home/u/kellerse/results/figures/box_plots-W4-11/boxplot_W4-11_{func}.pdf',
                format='pdf', bbox_inches='tight')

