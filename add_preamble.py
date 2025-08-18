#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import os


def add_preamble_to_python_files(directory, preamble):
    for root, _, files in os.walk(directory):
        for name in files:
            if not name.endswith(".py"):
                continue
            path = os.path.join(root, name)
            with open(path, "r", encoding="utf-8") as f:
                lines = f.readlines()

            i = 0
            while i < len(lines) and lines[i].startswith("#!"):
                i += 1
            while i < len(lines) and lines[i].strip() == "":
                i += 1

            new_lines = [preamble + "\n"] + lines[i:]
            if lines != new_lines:
                with open(path, "w", encoding="utf-8") as f:
                    f.writelines(new_lines)
                print(f"Updated {path}")


if __name__ == "__main__":
    add_preamble_to_python_files(
        os.getcwd(),
        "#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python",
    )

