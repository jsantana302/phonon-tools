#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import yaml
import numpy as np

def extract_acoustic_frequencies(file_path):
    with open(file_path, "r") as f:
        data = yaml.safe_load(f)

    # Find the Gamma point (q = [0,0,0])
    for point in data["phonon"]:
        if np.allclose(point["q-position"], [0.0, 0.0, 0.0], atol=1e-3):
            frequencies = [mode["frequency"] for mode in point["band"]]
            return sorted(frequencies)[:3]  # Get the lowest 3 frequencies (acoustic bands)

    return None  # Return None if Gamma point is not found

# Example usage
file_path = "band.yaml"  # Replace with your actual file path
acoustic_frequencies = extract_acoustic_frequencies(file_path)

if acoustic_frequencies:
    print("Acoustic band frequencies at Gamma point:", acoustic_frequencies)
else:
    print("Gamma point not found in the file.")

