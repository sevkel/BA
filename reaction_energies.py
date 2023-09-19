import os

def get_reaction(line):
    stoich_dict = {}
    stoich_dict['idx'] = line.split(":")[0]
    idx = stoich_dict['idx']
    molecs = line.split(":")[1].split()
    stoich_coeffs = line.split(":")[2].split()

    if (len(molecs) != len(stoich_coeffs)):
        raise Exception(f"number of reactants and products doesnt match in line {idx} -> check stoichiometry!")
    i = 0
    for molec in molecs:
        stoich_dict[molec] = stoich_coeffs[i]
        i += 1
    stoich_dict['e_react'] = 0
    return stoich_dict


if __name__ == '__main__':

    """
    change folders how to your results
    """

    main_dir = "./resultsDFT/molecs_YOURDATASET"
    result_dir = "./resultsDFT/e_react-YOURDATASET"
    stoich_dir = "./resultsDFT/stoichio_YOURDATASET.txt"

    functionals = ['DM21', 'DM21m', 'DM21mc', 'DM21mu', 'GGA_XC_PBE','MGGA_HM','GGA_HM', 'MGGA_XC_SCAN', 'lda', 'pbe', 'pbe0', 'r2scan', 'tpss']



    #for dir in os.listdir(main_dir):
    with open(stoich_dir, 'r') as file:
        lines = file.readlines()
        for line in lines:
            #print(line)
            dict = get_reaction(line)
            molecs = line.split(":")[1].split()
            #os.makedirs(result_dir + "/" + f"reaction-{dict['idx']}")
            outputfile = result_dir + "/" + f"reaction_energies-{dict['idx']}.txt"
            with open(outputfile, 'w') as otpfile:
                for func in functionals:
                    e_react = 0
                    for molec in molecs:
                        for moldir in os.listdir(main_dir):
                            if molec == moldir:
                                mol_path = os.path.join(main_dir, moldir)
                                for funcdir in os.listdir(mol_path):
                                    functional = funcdir.rsplit('_', 1)[0]
                                    if func == functional:
                                        func_path = os.path.join(mol_path, funcdir)
                                        en_filepath = os.path.join(func_path, "energies.txt")
                                        if os.path.isfile(en_filepath):
                                            with open(en_filepath, "r") as f_en:
                                                lines = f_en.readlines()
                                                for line in lines:
                                                    if "Total energy:" in line:
                                                        e_tot = float(line.strip('Total energy:').split()[0])
                                                        e_react = e_react + (int(dict[molec]) * e_tot)
                                                        break
                                                f_en.close()
                                        break
                                break
                    otpfile.write(f"{func}" + ": " + "{:.15f} kcal/mol\n".format(e_react * 627.5096080305927))
                otpfile.close()

        file.close()