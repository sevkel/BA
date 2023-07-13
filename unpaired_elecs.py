import os
import re


parent_directory = "/alcc/gpfs2/home/u/kellerse/testsystems/molecs_W4-11"
output_file = "unpaired_elecs.txt"


# get last number in line
def extract_number_after_minus(line):
    number = re.search(r'-(\d+)', line)
    if number:
        return int(number.group(1))
    else:
        return None


# open output file
with open(output_file, "w") as f_out:
    # run through all folders of mother folder
    for root, dirs, files in os.walk(parent_directory):
        for folder_name in dirs:
            folder_path = os.path.join(root, folder_name)
            control_file_path = os.path.join(folder_path, "control")

            # check if control exists
            if os.path.isfile(control_file_path):
                with open(control_file_path, "r") as f_control:
                    lines = f_control.readlines()
                    #print(folder_name)

                    for i in range(0, len(lines)):
                        line = lines[i]
                        found_unp = False

                        if "$closed shells" in line:
                            f_out.write(f"{folder_name}: 0\n")
                            break

                        if "$alpha shells" in line:
                            found_unp = True
                            alpha_line = lines[i + 1]
                            beta_line = lines[i + 3]
                            alpha_n = extract_number_after_minus(alpha_line)
                            beta_n = extract_number_after_minus(beta_line)



                        if found_unp == True:
                            difference = abs(beta_n - alpha_n)
                            f_out.write(f"{folder_name}: {difference}\n")



            else:
                f_out.write(f"{folder_name}: Datei 'control' nicht gefunden\n")
