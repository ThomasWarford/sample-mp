from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read, write
from urllib.request import urlretrieve
from pymatgen.core import Element

def create_toy_dataset(list_of_ase_atoms, num_materials, element_strings=None, low_count=None):
    if not element_strings:
        element_strings = [str(Element.from_Z(i+1)) for i in range(18)] # two rows of periodic table (to Argon)
    if not low_count:
        low_count = num_materials // len(element_strings) # this may be undesirable if num_materials is close to len(list_of_ase_atoms)
    
    idx_to_elements = {i: set(atoms.get_chemical_symbols()) for i, atoms in enumerate(list_of_ase_atoms)}

    element_to_idxs = {}
    for element in element_strings:
        element_to_idxs[element] = {i for i, atoms in idx_to_elements.items() if element in atoms}
    
    ## Sampling
    # add "rare" (low count) elements to training set
    sampled_idxs = set()
    for element in element_strings:
        if len(element_to_idxs[element]) < low_count:
            sampled_idxs.update(element_to_idxs[element])
    # add the remainder of the elements
    remaining_idxs = set(range(len(list_of_ase_atoms))).difference(sampled_idxs)
    sampled_idxs.update(np.random.choice(list(remaining_idxs), size=num_materials-len(sampled_idxs), replace=False))
    
    sampled_atoms = [list_of_ase_atoms[i] for i in sampled_idxs]

    return sampled_atoms
    
path_data = Path("./data")
path_data.mkdir(exist_ok=True)
url = "https://github.com/ACEsuit/mace-mp/releases/download/mace_mp_0b/mp_traj_combined.xyz"
path_mp_traj = path_data/url.split("/")[-1] 
if not path_mp_traj.exists():  
    urlretrieve(url, path_mp_traj)
    
mp_traj = read(path_mp_traj, index="0:")

sampled_atoms =create_toy_dataset(mp_traj, 10_000)

write("toy_dataset.xyz", sampled_atoms)