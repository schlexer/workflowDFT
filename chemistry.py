import ase
import numpy as np
import os
import sys

# Global variables.
transition_metals = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                     'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                     'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
metals = ['Li', 'Be', 'Na', 'Mg', 'K', 'Ca', 'Rb', 'Sr', 'Cs', 'Ba', 'Fr',
          'Ra', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
          'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th',
          'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
          'Cf', 'Es', 'Fm', 'Md', 'No', 'Lw', 'B', 'Al', 'Si', 'Ga', 'Ge',
          'In', 'Sn', 'Tl', 'Pb']
metals.extend(transition_metals)

def stoichiometry(atoms):
    """Return a number of layers given a slab.

    Parameters
    ----------
    atoms : object
        ASE atoms object.

    Returns
    -------
    num_dict : dictionary
        First entry is total number of atoms.
        Then key =  element and entry = number
    """
    symbols = atoms.get_chemical_symbols()
    num_dict = {}
    # num_dict['total'] = len(symbols)
    elements = list(set(symbols))
    for element in elements:
        count = symbols.count(element)
        num_dict[element] = count
    return(num_dict)

def is_metal(chemical_symbol):
    """Checks whether string is a metal elementary symbol.

    Parameters
    ----------
    chemical_symbol : string
        The element name.

    Returns
    -------
    metal : Boolean
        Whether it's a metal.
    """
    metal = False
    if chemical_symbol in metals:
        metal = True
    return(metal)

def is_oxide(atoms):
    """Checks whether atms object is an oxide.

    Parameters
    ----------
    atoms : object
        ASE atoms object.

    Returns
    -------
    oxide : Boolean
        Whether it is likely an oxide.
    """
    oxide = False
    num_dict = stoichiometry(atoms)
    elements = list(num_dict.keys())[1:]
    all_count = list(num_dict.values())[0]
    for element in elements:
        if element in transition_metals:
            if num_dict[element]/all_count*1.0 > 0.25:
                try:
                    if num_dict['O']/num_dict[element] >= 1:
                        # print("Materials class likely a metal oxide.")
                        oxide = True
                except KeyError as e:
                    pass
        elif element in metals:
            if num_dict[element]/all_count*1.0 > 0.25:
                try:
                    if num_dict['O']/num_dict[element] >= 1:
                        # print("Materials class likely a metal oxide.")
                        oxide = True
                except KeyError as e:
                    pass
    return(oxide)

if __name__ == '__main__':
    # TODO: Add testing functionality.
    pass
