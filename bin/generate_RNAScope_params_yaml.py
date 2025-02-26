import os
import yaml
import fire
from glob import glob


def generate_params_yaml(regex:str, out_name: str = 'output.yaml'):
    # Directory containing your .ome.tif files

    # Get a list of all .ome.tif files in the directory
    # files = [f for f in glob(regex) if f.endswith('.ome.tif')]
    files = [f for f in glob(regex) if not f.endswith('mask.tif')]
    root_dir = os.path.dirname(regex)

    # Create a list of dictionaries, each representing an image
    images = []
    for file in files:
        id = os.path.splitext(file)[0].split("/")[-1]  # Remove the file extension to get the id
        images.append([
            {'id': id.replace(" ", "_")},
            [os.path.join(root_dir, file)]
        ])

    # Create a dictionary representing the YAML file
    data = {'images': images, 'cell_diameters': [30], 'chs_to_call_peaks': [1,2]}

    # Write the data to a YAML file
    with open(out_name, 'w') as f:
        yaml.dump(data, f)

if __name__ == '__main__':
    fire.Fire(generate_params_yaml)  # This allows the script to be run from the command line