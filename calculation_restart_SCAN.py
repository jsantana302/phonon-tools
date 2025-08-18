#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os
import shutil
import argparse

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Restart CP2K calculation with specified parameters.")
    parser.add_argument('-xi', '--xc_init', required=False, default='PBE', help="Initial exchange correlation functional (e.g., PBE, SCAN, etc.)")
    parser.add_argument('-xf', '--xc_final', required=False, default='R2SCAN', help="Final exchange correlation functional (e.g., R2SCAN, SCAN, etc.)")
    parser.add_argument('-vi', '--vwd_init', required=False, default='D3', help="Initial van der Waals functional (e.g., D3, etc.)")
    parser.add_argument('-vf', '--vwd_final', required=False, default='D3', help="Final van der Waals functional (e.g., D4, etc.)")
    args = parser.parse_args()

    # Define paths and filenames based on input arguments
    xc_init = args.xc_init.upper()
    xc_final = args.xc_final.upper()
    vwd_init = args.vwd_init.upper()
    vwd_final = args.vwd_final.upper()


    source_wfn = f"../{xc_init}-{vwd_init}/{xc_init}-RESTART.wfn"
    source_restart = f"../{xc_init}-{vwd_init}/{xc_init}-1.restart"
    dest_wfn_file = f"{xc_final}-RESTART.wfn"
    dest_inp_file = f"{xc_final}.inp"
    job_file = f"/user/j.santanaandreo/u12658/files_CP2K/cp2k_{xc_final}.job"


    try:
        # Copy the required files
        shutil.copy(source_wfn, dest_wfn_file)
        print(f"Copied {source_wfn} to {dest_wfn_file}")
        shutil.copy(source_restart, dest_inp_file)
        print(f"Copied {source_restart} to {dest_inp_file}")
        shutil.copy(job_file, '.')

        # Log restart information to a file
        with open('restart_log.txt', 'w') as log_file:
            log_file.write(f"Restarting calculation from: {source_restart}\n")
            log_file.write(f"Using WFN file: {source_wfn}\n")
            log_file.write(f"New input file: {dest_inp_file}\n")

        # Modify the .inp file
        with open(dest_inp_file, 'r') as file:
            lines = file.readlines()

        wfn_line_found = False
        in_xc_functional = False
        modified_lines = []

        for line in lines:
            if line.strip().startswith("PROJECT_NAME"):
                new_line = f'   PROJECT_NAME "{xc_final}"\n'
                modified_lines.append(new_line)
            elif line.strip().startswith("RUN_TYPE"):
                modified_lines.append(f'   {line.strip()}\n')
            elif line.strip().startswith("REFERENCE_FUNCTIONAL"):
                if xc_final == "R2SCAN":
                    new_line = f'           REFERENCE_FUNCTIONAL "SCAN"\n'
                else:
                    new_line = f'           REFERENCE_FUNCTIONAL "{xc_final}"\n'
                modified_lines.append(new_line)
            elif line.strip().startswith("WFN_RESTART_FILE_NAME"):
                new_line = f'     WFN_RESTART_FILE_NAME {dest_wfn_file}\n'
                modified_lines.append(new_line)
                wfn_line_found = True
            elif line.strip().startswith("&DFT"):
                modified_lines.append(line)
                if not wfn_line_found:
                    new_line = f'   WFN_RESTART_FILE_NAME {dest_wfn_file}\n'
                    modified_lines.append(new_line)
                    wfn_line_found = True
            elif line.strip().startswith("&XC_FUNCTIONAL"):
                in_xc_functional = True
                if xc_final == "R2SCAN":
                    modified_lines.append(
                        "       &XC_FUNCTIONAL NO_SHORTCUT\n"
                        "         &MGGA_X_R2SCAN T\n"
                        "         &END MGGA_X_R2SCAN\n"
                        "         &MGGA_C_R2SCAN T\n"
                        "         &END MGGA_C_R2SCAN\n"
                        "       &END XC_FUNCTIONAL\n"
                    )
            elif in_xc_functional and line.strip().startswith("&END XC_FUNCTIONAL"):
                in_xc_functional = False
            elif not in_xc_functional:
                modified_lines.append(line)

        with open(dest_inp_file, 'w') as file:
            file.writelines(modified_lines)

        print("Files copied and .inp file modified successfully.")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()
