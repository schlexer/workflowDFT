import ase
import ast
import json
import pandas as pd
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import random
import re

import seaborn as sns
import sqlite3 as sql

from ase.build import molecule
from ase.thermochemistry import IdealGasThermo
from ase.thermochemistry import HarmonicThermo
from matplotlib.lines import Line2D

# Global parameters.
num_dict = {'0': '$_{0}$', '1': '$_{1}$', '2': '$_{2}$', '3': '$_{3}$', '4': '$_{4}$',
            '5': '$_{5}$', '6': '$_{6}$', '7': '$_{7}$', '8': '$_{8}$', '9': '$_{9}$'}
SUB = str.maketrans(num_dict)
# usage of num_dict example: "H2SO4".translate(SUB)

# Matplotlib settings.
sns.set_style('whitegrid')
plt.rc('text', usetex=True)
font = {'size':14}
plt.rc('font',**font)

colors = ['#068587', '#F2B134', '#ED553B', '#C36894', '#46698D', '#a6cee3', '#fdbf6f',
      '#b2df8a', '#1f78b4', '#e31a1c', '#fb9a99', '#33a02c', '#112F41']
edge_colors = []
for i in range(4):
    tmp = colors.copy()
    random.shuffle(tmp)
    for x in tmp:
        edge_colors.append(x)
colors = colors*4

# Physical constants.
cm2ev = 0.0001239841
electronicenergy = 0.0
default_pressure = 0.3235 #mbar
T_std = 273.15 # K
p_std = 1013.25 # mbar

# Dictionaries of molecular properties
molecule_dict = {
    # Pressures in mbar, electronic energies in eV, taken from DOI: 10.1039/C0EE00071J
    # Vibrations: TangRevised2018 dataset
        'H2': {'electronic_energy': 0.00,
               'geometry': 'linear',
               'pressure': 302.96,
               'spin': 1,
               'symmetrynumber': 2,
               'vibrations': [4329]},
        'O2': {'electronic_energy': 0.00,
                'geometry' : 'linear',
                'pressure': p_std,
                'spin' : 1,
                'symmetrynumber' : 2,
                'vibrations' : [1580]},
        'CO': {'electronic_energy': 1.75,
                'geometry' : 'linear',
                'pressure': 55.62,
                'spin' : 0,
                'symmetrynumber' : 2,
                'vibrations' : [2170]},
        'CO2': {'electronic_energy': 0.90,
                'geometry' : 'linear',
                'pressure': 1013.25,
                'spin' : 0,
                'symmetrynumber' : 2,
                'vibrations' : [1333.0, 2349.0, 667.0, 667.0]},
        'H2O': {'electronic_energy': 0.03,
                'geometry': 'nonlinear',
                'pressure': 35.34, # liquid water
                'spin': 0,
                'symmetrynumber': 2,
                'vibrations': [3843, 3721, 1625]},
        'CH4': {'electronic_energy': -1.22,
                'geometry' : 'nonlinear',
                'pressure': 20467,
                'spin' : 0,
                'symmetrynumber' : 2,
                'vibrations' : [2917, 1534, 1534, 3019, 3019, 3019, 1306, 1306, 1306]}}

# Vibrations of molecules on Cu(211) surface (TangRevised2018 data set)
mol_Cu211_dict = {'C': [206.5, 437.3, 667.1],
'CH': [375.2, 388.0, 508.7, 597.0, 604.4, 3064.2],
'CO': [150.2, 158.6, 189.6, 314.6, 1879.2],
'H2COH': [134.7, 199.3, 329.8, 385.9, 565.1, 837.4, 1083.4, 1111.6, 1316.6, 1438.7, 3011.9, 3091.1, 3698.4],
'OCC': [206.5, 437.3, 667.1],
'CH2': [324.2, 348.1, 430.7, 551.3, 636.3, 1335.6, 2978.9, 3056.6],
'CHO': [86.3, 144.1, 207.5, 249.5, 417.8, 679.4, 1244.0, 1527.3, 2810.2],
'HCOO': [104.2, 171.3, 258.3, 291.8, 337.2, 734.3, 1004.1, 1307.4, 1344.9, 1509.6, 2967.0],
'HCOOH': [37.6, 76.8, 131.3, 165.7, 192.1, 622.5, 693.1, 1022.6, 1086.1, 1279.1, 1377.6, 1737.1, 3017.0, 3619.5],
'CH3': [101.5, 198.0, 312.1, 561.0, 570.6, 1124.0, 1405.5, 1406.2, 2951.1, 3011.3, 3021.4],
'HCOH': [24.7, 110.9, 171.7, 230.2, 438.8, 503.1, 952.9, 1193.8, 1248.7, 1424.7, 3003.4, 3438.9],
'COH': [105.4, 181.7, 252.7, 308.0, 341.6, 417.5, 1110.6, 1230.2, 3599.0],
'COOH': [137.9, 157.8, 184.6, 240.3, 375.2, 625.4, 633.4, 942.1, 1233.0, 1663.1, 3510.4]}

class Adsorbate:
    def __init__(self, molecule_name='C'):
        self.molecule_dict = mol_Cu211_dict
        self.name = molecule_name
        self.vibrations = self.molecule_dict[self.name]
        self.vib_energies = np.asarray(self.vibrations) * cm2ev

    def get_helmholtz_energy(self, temperature=T_std, electronic_energy=0, verbose=False):
        """Returns the Helmholtz energy of an adsorbed molecule.

        Parameters
        ----------
        temperature : numeric
            temperature in K
        electronic_energy : numeric
            energy in eV
        verbose : boolean
            whether to print ASE thermochemistry output

        Returns
        -------
        helmholtz_energy : numeric
            Helmholtz energy in eV
        """
        thermo_object = HarmonicThermo(vib_energies=self.vib_energies,
                                       potentialenergy=electronic_energy)
        self.helmholtz_energy = thermo_object.get_helmholtz_energy(temperature=temperature, verbose=verbose)
        return (self.helmholtz_energy)

    def get_internal_energy(self, temperature=T_std, electronic_energy=0, verbose=False):
        """Returns the internal energy of an adsorbed molecule.

        Parameters
        ----------
        temperature : numeric
            temperature in K
        electronic_energy : numeric
            energy in eV
        verbose : boolean
            whether to print ASE thermochemistry output

        Returns
        -------
        internal_energy : numeric
            Internal energy in eV
        """
        thermo_object = HarmonicThermo(vib_energies=self.vib_energies,
                                       potentialenergy=electronic_energy)
        self.internal_energy = thermo_object.get_internal_energy(temperature=temperature, verbose=verbose)
        return (self.internal_energy)

class GasMolecule:
    def __init__(self, molecule_name='O2'):
        self.molecule_dict = molecule_dict
        self.name = molecule_name
        self.electronic_energy = molecule_dict[self.name]['electronic_energy']
        self.atom_object = ase.build.molecule(self.name)
        self.pressure = molecule_dict[self.name]['pressure']
        self.geometry = molecule_dict[self.name]['geometry']
        self.symmetrynumber = molecule_dict[self.name]['symmetrynumber']
        self.spin = molecule_dict[self.name]['spin']
        self.vibrations = molecule_dict[self.name]['vibrations']

    def get_free_energy(self, temperature = T_std, p_mbar = 'Default', electronic_energy = 'Default'):
        """Returns the internal energy of an adsorbed molecule.

        Parameters
        ----------
        temperature : numeric
           temperature in K
        electronic_energy : numeric
           energy in eV
        p_mbar : numeric
           pressure in mbar

        Returns
        -------
        internal_energy : numeric
           Internal energy in eV
        """
        if not temperature or not p_mbar: # either None or 0
            return(0)
        else:
            if electronic_energy == 'Default':
                electronic_energy = molecule_dict[self.name]['electronic_energy']
            else:
                pass
            if p_mbar == 'Default':
                p_mbar = molecule_dict[self.name]['pressure']
            else:
                pass
            pressure = p_mbar * 100  # gives Pa
            ideal_gas_object = IdealGasThermo(vib_energies=self.get_vib_energies(),
                                              potentialenergy=electronic_energy,
                                              atoms=self.atom_object,
                                              geometry=molecule_dict[self.name]['geometry'],
                                              symmetrynumber=molecule_dict[self.name]['symmetrynumber'],
                                              spin=molecule_dict[self.name]['spin'])
            # Returns the Gibbs free energy, in eV, in the ideal gas
            # approximation at a specified temperature (K) and pressure (Pa).
            energy = ideal_gas_object.get_gibbs_energy(temperature=temperature,
                                                       pressure=pressure, verbose=False)
            self.free_energy = energy
            return(self.free_energy)

    def get_enthalpy(self, temperature = T_std, electronic_energy = 'Default'):
        """Returns the internal energy of an adsorbed molecule.

        Parameters
        ----------
        temperature : numeric
        temperature in K
        electronic_energy : numeric
        energy in eV


        Returns
        -------
        internal_energy : numeric
        Internal energy in eV
        """
        if not temperature: # either None or 0
            return(0, 0, 0)
        if electronic_energy == 'Default':
            electronic_energy = molecule_dict[self.name]['electronic_energy']
        else:
            ideal_gas_object = IdealGasThermo(vib_energies=self.get_vib_energies(),
                                              potentialenergy=electronic_energy,
                                              atoms=self.atom_object,
                                              geometry=molecule_dict[self.name]['geometry'],
                                              symmetrynumber=molecule_dict[self.name]['symmetrynumber'],
                                              spin=molecule_dict[self.name]['spin'])
            energy = ideal_gas_object.get_enthalpy(temperature=temperature, verbose=False)
            self.enthalpy = energy
            return(self.enthalpy)

    def get_vib_energies(self):
        """Returns a list of vibration in energy units eV.

        Returns
        -------
        vibs : list of vibrations in eV
        """
        vibs = self.molecule_dict[self.name]['vibrations']
        vibs = np.array(vibs) * cm2ev
        return (vibs)

def get_ZPE(viblist):
    """Returns the zero point energy from a list of frequencies.

    Parameters
    ----------
    viblist : List of numbers or string of list of numbers.

    Returns
    -------
    ZPE : Zero point energy in eV.
    """
    if type(viblist) is str:
        l = ast.literal_eval(viblist)
    else:
        l = viblist
    l = [float(w) for w in l]
    ZPE = 0.5*sum(l)*cm2ev
    return(ZPE)

def auto_labels(df):
    """Transforms atomic system information into well-formatted labels.

    Parameters
    ----------
    df : Pandas DataFrame.

    Returns
    -------
    labels : list of system labels.
    """
    systems = list(df.system)
    facets = list(df.facet)
    systems_labels = [w.replace('_', '\ ') for w in systems]
    systems_labels = [w.translate(SUB) for w in systems_labels]
    systems_labels = [w.replace('}$$_{', '') for w in systems_labels]
    systems_labels = [w.replace('$', '') for w in systems_labels]
    systems_labels = ['$' + w + '$' for w in systems_labels]

    facets_label = [w.replace('_', '\ ') for w in facets]
    facets_label = ['(' + w + ')' for w in facets_label]

    labels = []
    for i, sys in enumerate(systems_labels):
        labels.append(sys + facets_label[i])
    # labels = list(set(labels))
    return(labels)

def proton_hydroxide_free_energy(pH=0, temperature=T_std):
    """Returns the Gibbs free energy of proton in bulk solution.

    Parameters
    ----------
    pH : pH of bulk solution
    temperature : numeric
        temperature in K
    p_mbar : numeric
       pressure in mbar

    Returns
    -------
    G_H, G_OH : Gibbs free energy of proton and hydroxide.
    """
    R = 8.314472 # J/(mol K)
    F = 9.64853399*(10**4) # C/mol
    z = 1.0 # n_electrons
    ln10 = 2.30258509 # log10 for pH
    H2 = GasMolecule('H2')
    H2O = GasMolecule('H2O')
    G_H2 = H2.get_free_energy(temperature = temperature, p_mbar = 1013.25)
    G_H2O = H2O.get_free_energy(temperature = temperature)
    G_H = (0.5*G_H2) - ((R*temperature)/(z*F))*ln10*pH
    G_OH = G_H2O - G_H  # Do not need Kw when water equilibrated
    return(G_H, G_OH)

def get_FEC(molecule_list, temperature, pressure_mbar='Default', electronic_energy='Default'):
    """Returns the Gibbs free energy corrections to be added to raw reaction energies.

    Parameters
    ----------
    molecule_list : list of strings
    temperature : numeric
        temperature in K
    p_mbar : numeric
       pressure in mbar

    Returns
    -------
    G_H, G_OH : Gibbs free energy of proton and hydroxide.
    """
    if not temperature or not pressure_mbar:
        return(0)
    else:
        molecule_list = [m for m in molecule_list if m != 'star']
        # print(molecule_list)
        FEC_sum = []
        for molecule in molecule_list:
            if 'gas' in molecule:
                mol = GasMolecule(molecule.replace('gas', ''))
                if pressure_mbar == 'Default':
                    p = mol.pressure
                else:
                    p = pressure_mbar
                if electronic_energy == 'Default':
                    ee = mol.electronic_energy
                else:
                    ee = electronic_energy
                FEC = mol.get_free_energy(temperature=temperature, p_mbar=p, electronic_energy = ee)
                FEC_sum.append(FEC)
            if 'star' in molecule:
                FEC = Adsorbate(molecule.replace('star', ''))
                FEC = FEC.get_helmholtz_energy(temperature=temperature)
                FEC_sum.append(FEC)
        FEC_sum = sum(FEC_sum)
    return (FEC_sum)


def reaction_scheme(df, reaction_intermediates: list = ['COgas', 'COstar', 'CHOstar'],
                    potential=None, pH=None, temperature=T_std, pressure_mbar='Default'):
    """Returns a dataframe with Gibbs free reaction energies.

    Parameters
    ----------
    df : Pandas DataFrame generated by db_to_df
    temperature : numeric
        temperature in K
    p_mbar : numeric
       pressure in mbar
    pH : PH in bulk solution
    potential : Electric potential vs. SHE in eV

    Returns
    -------
    df : DataFrame suitable for plotting.
    """
    # set reaction scheme
    reactions = reaction_intermediates

    # set reaction labels for plotting
    reaction_labels = [w.replace('gas', '$_g$') for w in reactions]
    reaction_labels = [w.replace('star', '*') for w in reaction_labels]
    reaction_scheme = reactions[1:]
    reaction_number = np.arange(len(reactions))

    # get unique substrates
    labels = list(df.labels)

    # Gather and summarize data from data base.
    unique_system_label = []
    surface = []
    facet = []
    intermediate_label = []
    reaction_coordinate = []
    energies = []

    # compute energetics
    G_PE = 0.0
    if potential is not None and pH is not None:
        ne = 1.0
        G_PE = -ne * potential + \
               proton_hydroxide_free_energy(pH=pH, temperature=temperature)[0]
        # print('G_PE: '+ str(G_PE))
    for lab in list(set(labels)):
        proton_transfers = 0
        energies_system = [0.0]
        reaction_coordinate_system = [[x, x + 0.5] for x in reaction_number]
        reaction_coordinate.append(reaction_coordinate_system)
        for step, products in enumerate(reaction_scheme):
            df_tmp = df[df['labels'].apply(lambda x: lab in x)]
            df_tmp = df_tmp[df_tmp['products'].apply(lambda x: products in x)]

            if step == 0:
                unique_system_label.append(df_tmp.iloc[0]['labels'])
                surface.append(df_tmp.iloc[0]['system'])
                facet.append(df_tmp.iloc[0]['facet'])
                intermediate_label.append(reaction_labels)

            all_reactants = list(df_tmp['reactants'])[0]
            all_products = list(df_tmp['products'])[0]
            reactant_FEC = get_FEC(all_reactants, temperature, pressure_mbar)
            product_FEC = get_FEC(all_products, temperature, pressure_mbar)
            FEC = product_FEC - reactant_FEC
            # print('Total FEC: '+str(FEC))

            if 'H2gas' in df_tmp.iloc[0]['reactants']:
                proton_transfers += 1
            reaction_energy = df_tmp.iloc[0]['reaction_energy'] + FEC - proton_transfers * G_PE
            # print(df_tmp.iloc[0]['reaction_energy'])
            # print(reaction_energy)
            energies_system.append(reaction_energy)
        energies.append(energies_system)

    # Dump data into data frame.
    df_new = pd.DataFrame()
    df_new['system_label'] = unique_system_label
    df_new['system'] = surface
    df_new['facet'] = facet
    df_new['intermediate_labels'] = intermediate_label
    df_new['reaction_coordinate'] = reaction_coordinate
    df_new['reaction_energy'] = energies
    df_new = df_new.sort_values(by=['facet', 'system'])
    df_new = df_new.reset_index(drop=True)
    return (df_new)


def plot_reaction_scheme(df, show=False, potential = None, pH = None,
                         temperature=T_std, pressure_mbar='Default', e_lim=None):
    """Returns a matplotlib object with the plotted reaction path.

    Parameters
    ----------
    df : Pandas DataFrame generated by reaction_network
    temperature : numeric
        temperature in K
    p_mbar : numeric
       pressure in mbar
    pH : PH in bulk solution
    potential : Electric potential vs. SHE in eV
    e_lim: Limits for the energy axis.

    Returns
    -------
    df : DataFrame suitable for plotting.
    """
    ncols = int((df.shape[0]/20)) +1
    fig_width = ncols + 2*len(df['intermediate_labels'][0])
    figsize = (fig_width, 5)
    fig, ax = plt.subplots(figsize=figsize)

    lines = []
    for j, energy_list in enumerate(df['reaction_energy']):
        R = df['reaction_coordinate'][j]
        E = [[x, x] for x in energy_list]
        labels = df['system_label']

        for i, n in enumerate(R):
            if i == 0:
                line = Line2D([0], [0], color=colors[j], lw=4)
                lines.append(line)
                ax.plot(n, E[i], ls='-', color=colors[j], linewidth=3.25, solid_capstyle='round',
                        path_effects=[pe.Stroke(linewidth=6, foreground=edge_colors[j]), pe.Normal()],
                        label=labels[j])
            else:
                ax.plot(n, E[i], ls='-', color=colors[j], linewidth=3.25, solid_capstyle='round',
                        path_effects=[pe.Stroke(linewidth=6, foreground=edge_colors[j]), pe.Normal()])
            if i < len(R) - 1:
                ax.plot([n[1], n[1] + 0.5], [E[i], E[i + 1]], ls='--',
                        dashes=(3, 2), color=colors[j], linewidth=1.)

    ax.legend(handlelength=0.4, ncol=ncols, loc=2,
              bbox_to_anchor=(1.05, 1), borderaxespad=0., fontsize=12)
    # ax.set_xticklabels(reaction_number,reaction_labels, rotation=45)
    if e_lim:
        ax.set_ylim(e_lim)
    ax.set_xlabel('Reaction coordinate')
    ax.set_ylabel('Reaction free energy (eV)')
    reaction_labels = df['intermediate_labels'][0]
    reaction_labels = [w.translate(SUB) for w in reaction_labels]
    plt.xticks(np.arange(len(reaction_labels)) + 0.25, tuple(reaction_labels), rotation=45)
    # plt.tight_layout()
    if potential is not None and pH is not None:
        ax.set_title('U = '+str(potential)+' eV vs. SHE; pH = '+str(pH)+'; \n T = '+str(temperature)+' K; p = '+str(pressure_mbar)+' mbar')
        # ax.set_title('U = -0.5 eV vs. RHE; pH = '+str(pH)+'; \n T = '+str(temperature)+' K; p = '+str(pressure_mbar)+' mbar')

    else:
        ax.set_title('T = '+str(temperature)+'; p = '+str(pressure_mbar)+' mbar')
    if show:
        plt.show()
    return(fig)


def select_data(db_file, slab=None, facet=None):
    """Gathers relevant data from SQL database generated by CATHUB.

    Parameters
    ----------
    db_file : Path to database
    slab : Which metal (slab) to select.
    facet : Which facets to select.

    Returns
    -------
    data : SQL cursor output.
    """
    con = sql.connect(db_file)
    cur = con.cursor()
    if slab and facet:
        select_command = 'select chemical_composition, facet, reactants, products, reaction_energy ' \
                     'from reaction where facet='+str(facet)+' and chemical_composition LIKE "%'+slab+'%";'
    elif slab and not facet:
        select_command = 'select chemical_composition, facet, reactants, products, reaction_energy ' \
             'from reaction where chemical_composition LIKE "%'+slab+'%";'
    else:
        select_command = 'select chemical_composition, facet, reactants, products, reaction_energy from reaction;'
    cur.execute(select_command)
    data = cur.fetchall()
    return(data)


def db_to_df(db_file, slabs=None, facet=None):
    """Transforms database to data frame.

    Parameters
    ----------
    db_file : Path to database
    slabs : Which metals (slabs) to select.
    facet : Which facets to select.

    Returns
    -------
    df : Data frame.
    """
    systems = []
    data = []
    if slabs:
        for slab in slabs:
            data_tmp = select_data(db_file, slab=slab, facet=facet)
            data.append(data_tmp)
            subsystem = [tup[0] for i, tup in enumerate(data_tmp)]
            systems.append(list(set(subsystem))[0])
    else:
        data_tmp = select_data(db_file)
        data.append(data_tmp)

    df = pd.DataFrame()
    system, facet, reactants, products, reaction_energy = [], [], [], [], []
    for entry in data:
        for reaction in entry:
            system.append(str(reaction[0]))
            facet.append(str(reaction[1]))
            reactants_i = [molecule for molecule in ast.literal_eval(reaction[2]).keys()]
            reactants.append(reactants_i)
            products_i = [molecule for molecule in ast.literal_eval(reaction[3]).keys()]
            products.append(products_i)
            reaction_energy.append(float(reaction[4]))

    df[0] = system
    df[1] = facet
    df[2] = reactants
    df[4] = products
    df[5] = reaction_energy

    df.columns = ['system', 'facet', 'reactants', 'products', 'reaction_energy']
    labs = auto_labels(df)
    df['labels'] = labs
    df = df.sort_values(by=['facet', 'system'])
    df = df.reset_index(drop=True)
    return(df)

def unique_reactions(df):
    """Identifies unique elementary reactions in data frame.

    Parameters
    ----------
    df : Data frame.

    Returns
    -------
    reaction_list : List of unique elementary reactions.
    """
    reaction_list =[]
    for idx, entry in enumerate(df['reactants']):
        reaction = []
        for x in entry:
            reaction.append(x)
        reaction.append('-->')
        for y in df['products'][idx]:
            reaction.append(y)
        reaction_list.append(reaction)
    string_list = [str(reaction) for reaction in reaction_list]
    string_list = sorted(list(set(string_list)))
    reaction_list = [ast.literal_eval(entry) for entry in string_list]
    return(reaction_list)
