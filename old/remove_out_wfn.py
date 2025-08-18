#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os

def clean_directories(base_dir):
    for root, dirs, files in os.walk(base_dir):
        # Skip the base directory
        if root == base_dir:
            continue

        # Check if any file contains 'forces' in its name
        if not any('forces' in file_name for file_name in files):
            for file_name in files:
                if file_name.endswith(('wfn', 'wfn.bak-1', 'out')):
                    file_path = os.path.join(root, file_name)
                    try:
                        os.remove(file_path)
                        print(f"Removed: {file_path}")
                    except Exception as e:
                        print(f"Failed to remove {file_path}: {e}")

if __name__ == "__main__":
    base_directory = os.getcwd()  # Set to current directory
    clean_directories(base_directory)

