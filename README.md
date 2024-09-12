# Protein Membrane Visualization

This project aims to identify and visualize potential membrane regions in protein structures by analyzing hydrophobicity and spatial arrangements. It utilizes the PyMOL visualization tool to display these regions.

## Project Structure

The project consists of the following main components:

1. **`protein.py`**: Defines the `Protein` class for loading and processing protein structures from PDB files.
2. **`membrane_finder.py`**: Contains the `MembraneFinder` class for identifying the best membrane planes in a protein structure.
3. **`geometry_and_axis.py`**: Includes definitions for geometric entities like `Point`, `Vector`, `Plane`, and `Axis`.
4. **`visualization.py`**: Provides functions for visualizing the protein structure and planes in PyMOL.
5. **`main.py`**: The main script to run the analysis and visualization.

## Dependencies

Ensure you have the following dependencies installed:

- `numpy`
- `biopython`
- `pymol` (PyMOL is required to be installed separately; this project assumes PyMOL is available through its Python API)

You can create an environment with the required packages using the provided `environment.yml`.

## Installation

1. **Clone the repository:**
    ```bash
    git clone git@github.com:SanaGUEDOUAR/PM_Visualization.git
    cd PM_Visualization/
    ```

2. **Create and activate the conda environment:**
    ```bash
    conda env create -f environment.yml
    conda activate tm_parts-env
    ```

## Usage

To analyze a protein structure and visualize potential membrane regions, follow these steps:

1. **Run the main script with your PDB file:**
    ```bash
    python main.py path/to/your_protein.pdb
    ```

    Replace `path/to/your_protein.pdb` with the path to your protein structure file in PDB format.

2. **View results in PyMOL:**
    - The session file with a `.pse` extension will be saved in the `results` directory. You can open this file in PyMOL to visualize the protein structure and identified planes.

## Files

- **`protein.py`**: Contains the `Protein` class for loading and analyzing protein structures.
- **`membrane_finder.py`**: Implements the `MembraneFinder` class for identifying optimal membrane planes.
- **`geometry_and_axis.py`**: Defines geometric and axis-related classes.
- **`visualization.py`**: Provides functions for visualizing planes and protein structures in PyMOL.
- **`main.py`**: Main script for executing the analysis and generating visualizations.
- **`environment.yml`**: Conda environment configuration file.

## Example

After running the main script with a PDB file, you will see output indicating the best hydrophobicity score and the found axis. The PyMOL session will be saved in the `results` directory and can be loaded for visualization.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

- **Biopython**: For tools to handle biological data.
- **PyMOL**: For visualization of molecular structures.
