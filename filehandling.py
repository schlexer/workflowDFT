import pandas as pd
import os

def delete_large_files(workdir='cwd', nfiles=10):
    """Delete the n largest files in the directory.

    Parameters
    ----------
    workdir : path to starting dir

    """
    # Check which files are large.
    if workdir == 'cwd':
        workdir = os.getcwd()
    file_names = [os.path.join(path, name) for path, _, filenames in os.walk(workdir)
            for name in filenames]
    file_sizes = [(os.path.getsize(name), name) for name in file_names]
    file_sizes = sorted(file_sizes, key=lambda tup: tup[0])
    largest = file_sizes[-nfiles:]
    # delete files.
    for tup in largest:
        #print(tup[1])
        os.remove(tup[1])
    return(None)

def file_to_df(file_name):
    """Read in file and return pandas data_frame.

    Parameters
    ----------
    file_name : file must be in the same folder

    Returns
    -------
    df : pandas data frame
    """
    filename, file_extension = os.path.splitext(file_name)
    if file_extension=='.csv':
        df = pd.read_csv(file_name, sep=',', header=0)
    elif file_extension=='.tsv':
        df = pd.read_csv(file_name, sep='\t', header=0)
    else:
        print('Please provide csv or tsv file format.')
    return(df)

def find_dirs(workdir='cwd', dir_contains=''):
    """Return a list of paths that contain desired folder keyword
    Parameters
    ----------
        workdir : path to starting dir
        dir_contains : Filter only for dirs that contain keyword
    """
    if workdir == 'cwd':
        workdir = os.getcwd()
    dirs = []
    for dirpaths, dirnames, fnames in os.walk(workdir):
        if dir_contains in dirpaths:
            dirs.append(dirpaths+'/')
    return(dirs)

def find_files(workdir='cwd', dir_contains='', filename=None, file_contains=''):
    """Return a list of paths that contain desired file
    Parameters
    ----------
        workdir : path to starting dir
        dir_contains : Filter only for dirs that contain keyword
        filename: exact match for filename
        file_contains: the file name contains this keyword
    """
    if workdir == 'cwd':
        workdir = os.getcwd()
    fpaths, dirs = [], []
    for dirpaths, dirnames, fnames in os.walk(workdir):
        dirs.append(dirpaths)
        for fname in fnames:
            if filename is not None:
                if fname == filename:
                    fpaths.append(dirpaths + '/' + fname)
            else:
                if file_contains in fname:
                    fpaths.append(dirpaths+'/'+fname)
    fpaths = [w for w in fpaths if dir_contains in w]
    return(fpaths)

def save_slab(atoms, slab_element=None, facet=None, adsorbate=None, isomer=0, local_path=None):
    """ Save the structure and automatically find out where.

    Parameters
    ----------
        atoms : object
        slab_element : string of slab name or element
        facet : Miller indices of surface
        adsorbate : strin of adsorbate name
        local path : costumize your folder name
    """
    if slab_element and facet and adsorbate:
        if not isomer:
            traj = adsorbate+'_'+slab_element+str(facet)+'.json'
            folder = os.getcwd()+'/'+slab_element+'/'+str(facet)+'/'+adsorbate+'/'
        else:
            traj = adsorbate+'_'+slab_element+str(facet)+'.json'
            folder = os.getcwd()+'/'+slab_element+'/'+str(facet)+'/'+adsorbate+'/isomer'+str(isomer)+'/'
    elif slab_element and facet and not adsorbate:
        traj = slab_element+str(facet)+'.json'
        folder = os.getcwd()+'/'+slab_element+'/'+str(facet)+'/clean/'
    else:
       traj = "".join(list(set(atoms.get_chemical_symbols())))+'.json'
       folder = os.getcwd()+'/'+local_path
    print(folder)
    if not os.path.exists(folder):
        os.makedirs(folder)
    ase.io.write(folder+traj,atoms)
    return(atoms)

if __name__ == '__main__':
    cwd = os.getcwd()+'/'
    file_to_df(cwd+'test/data.tsv')
    find_files()
    find_dirs()
