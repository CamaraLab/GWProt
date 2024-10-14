import os
import math
import numpy as np
from subprocess import Popen, PIPE, STDOUT, DEVNULL
import sys




from .GW_scripts import *
from .read_pdb import *
from .GW_protein import *
from .GW_stress import *
from .weighted_alignment import *



class my_pymolPy3:
    # pretty much the one from online
    # https://github.com/carbonscott/pymolPy3
    def __init__(self):



        ### "pymol --version" gives the version
        V =  Popen(
            ["pymol", "--version"], 
            shell=False,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE,
            universal_newlines=True,
        )
        self.version = V.stdout.read().split()[1]
        #print(self.version)
        V.communicate()
        V.terminate()

        self.pymolpy3 = Popen('pymol -pc', shell=True, stdin=PIPE,stdout=DEVNULL, stderr = DEVNULL)
        #print('made')
        
        # self.pymolpy3 = Popen(
        #     ["pymol", "-pc"], # -p means get commands from stdin #  -K to continue going
        #     shell=False,
        #     stdin=PIPE,
        #     stdout=PIPE,
        #     stderr=PIPE,
        #     universal_newlines=True,
        # )




    def __del__(self): 
        # Turn off stdin...
        self.pymolpy3.stdin.close()
        
        #print('closed')
        # Wait until program is over...
        self.pymolpy3.wait()
        #self.pymolpy3.communicate() #for stdout = PIPE
        #print('waited')

        # Terminate the subprcess...
        self.pymolpy3.terminate()
        #print('terminated')


    def __call__(self, s):
        # Keep reading the new pymol command as a byte string...
        #self.pymolpy3.stdin.write( s + "\n")
        self.pymolpy3.stdin.write( bytes( (s + '\n').encode() ) )

        # Flush it asap...
        #self.pymolpy3.stdin.flush()


        return 0


    def quit(self):
        self("quit")

    def get_err(self):
        return self.pymolpy3.stderr.read()



    def communicate(self, obj: str):
        # self.pymolpy3.stdin.write(f"print {obj} \n")
        # self.pymolpy3.stdin.flush()
        # str = self.pymolpy3.stdout.read()
        # err = self.pymolpy3.stderr.read()
        str, err = self.pymolpy3.communicate(f"print {obj} \n")
        return str, err


def compare_proteins_in_pymol(file1:str , file2:str, output_file:str, chain1:str = 'A', chain2: str = 'A',
     transport_plan: np.array = None, threshold:float =0.5)->None:
    """
    This loads two pdb files and display them in Pymol and aligns them with a transport plan, then saves the scene to a .pse file.
    A rigid alignment is created minimizing weighted RSMD. Note that if Pymol 2 is used it uses ``cmd.cealign`` instead. 
    For a pair of aligned residues, a line will connect them if over ``threshold`` of each of their mass is connected. The proteins are also colored by the stress levels.
    
    :param file1: Filepath to the first protein.
    :param file2: Filepath to the second protien.
    :param output_file: Filepath where the resulting file should be saved to.
    :param chain1: Which chain of the first protein to use, default is ``A``.
    :param chain2: Which chain of the second protein to use, default is ``A``.
    :param transport_plan: A transport plan to align the two proteins. If none is provided one will be calculated with ``GW_protein.run_GW``.
    :param threshold: The threshold for displaying aligned residues. 

    """
    p1 = GW_protein.make_protein_from_pdb(file1,chain_id = chain1)
    p2 = GW_protein.make_protein_from_pdb(file2 ,chain_id = chain2)

    if transport_plan is  None:
        c, transport_plan = GW_protein.run_GW(p1,p2, transport_plan=True) 


    assert transport_plan.shape == (len(p1), len(p2))
    assert math.isclose(np.sum(transport_plan),1)


    ps = get_pairing(transport_plan, threshold0=threshold, threshold1=threshold)

    pret, rot, trans = weighted_RMSD(p1.coords, p2.coords, transport_plan)

    ll = pymol_transform(pretrans=pret, rot=rot, posttrans=trans)

    stress1, stress2 = GW_protein.GW_stress(p1,p2, transport_plan)
 
    pm = my_pymolPy3()

    pm(f"cmd.load('{file1}', 'prot1')")
    pm(f"cmd.load('{file2}', 'prot2')")
    pm(f"cmd.hide('cartoon','prot1 and not /prot1//{chain1}' )")
    pm(f"cmd.hide('cartoon','prot2 and not /prot2//{chain2}' )")

    pm(f"residue_numbers1 = list(set(atom.resi_number for atom in cmd.get_model('/prot1//{chain1}').atom))")
    #print(f"residue_numbers1 = list(set(atom.resi_number for atom in cmd.get_model('/prot1//{chain1}').atom))")
    pm(f"residue_numbers2 = list(set(atom.resi_number for atom in cmd.get_model('/prot2//{chain2}').atom))")
    #print(f"residue_numbers2 = list(set(atom.resi_number for atom in cmd.get_model('/prot2//{chain2}').atom))")
    pm(f"stress1 = {str(stress1)}")
    pm(f"stress2 = {str(stress2)}")
    pm(f"cmd.alter( '/prot1//{chain1}',  'b = stress1[residue_numbers1.index(int(resi))]' )")
    pm(f"cmd.alter( '/prot2//{chain2}',  'b = stress2[residue_numbers2.index(int(resi))]' )")
    pm(f"cmd.spectrum(expression='b', selection='/prot1//{chain1}', palette='red_lime', byres=1)")
    pm(f"cmd.spectrum(expression='b', selection='/prot2//{chain2}', palette='red_marine', byres=1)")


    pm(f"ps = {str(ps)}")
    pm("chain1 = '" + chain1 +"'")
    pm("chain2 = '" + chain2 +"'")


    path =  os.path.abspath(__file__).removesuffix('pymol_protein_viewer.py')
    script = os.path.join(path, '_compare_proteins_in_pymol_helper.py')
    pm(f"run {script}")
    

    pm("cmd.set('dash_gap', '0')")
    pm("cmd.set('dash_width', '3')")
    pm("cmd.zoom('all')")
    pm("cmd.bg_color('grey80')")
    if "2." in pm.version:
        raise ValueError("Pymol 2 can segmentation fault when running transform_selection, cmd.cealign used instead")
        pm("cmd.cealign('prot1', 'prot2')")

    pm(f"cmd.transform_selection( 'prot2' , matrix =  {ll})")
    pm(f"cmd.center('/prot1//{chain1} and /prot2//{chain2}')")
    pm(f"cmd.zoom('/prot1//{chain1} and /prot2//{chain2}')")
    pm(f"cmd.save( '{output_file}') ")
    #return ps


def show_proteins_with_values(infiles: list[str], chain_ids : list[str],  data_lists: list[float], output_file: str, hide: bool = True)->None:
    """    
    This loads pdb files and display them in Pymol with colors based on the ``data_lists``, then saves the scene to a .pse file.
    
    :param infiles: Filepaths to the first protein.
    :param chain_ids: Which chains of the proteins to use, a value must be entered for each.
    :param output_file: Filepath where the resulting file should be saved to.
    :param hide: Whether to hide the chains that aren't selected.
    
    """
    if len(infiles) != len(chain_ids):
        raise ValueError('The number of input files must match the number of chain_ids given. The latter can contain None')
    if len(infiles) != len(data_lists):
        raise ValueError('The number of input files must match the number of lists given.')

    
    n = len(data_lists)
    prots = [GW_protein.make_protein_from_pdb(infiles[i], chain_id = chain_ids[i]) for i in range(n)]

    pm = my_pymolPy3()
    #print('pm created')


    for i in range(n):
        #print(i,' started')
        if len(prots[i]) != len(data_lists[i]):
            raise ValueError(f'length of protein {i} does not match length of data {i}. This could be caused by missing data in a PDB file.')
        pm(f"cmd.load('{infiles[i]}', 'prot' + str({i}))")
        #print(i,' loaded')
        pm(f"residue_numbers{i} = list(set(atom.resi_number for atom in cmd.get_model('/prot{str(i)}//{chain_ids[i]}').atom))")
        #print(i, 'resnums set')
        pm(f"new_b{i} = {str(list(data_lists[i]))}")
        #print(i, 'b set')
        pm(f"cmd.alter( selection='/prot' + str({i}) + '//{chain_ids[i]}', expression='b = new_b' + str({i}) + '[residue_numbers{i}.index(int(resi))]')")
        #print(i, 'altered')
        pm(f"cmd.spectrum(expression='b', selection='/prot' + str({i}) + '//{chain_ids[i]}', palette='yellow_red', byres=1)")
        #print(i, 'colored')
        if i != 0:
            pm(f"cmd.cealign( '/prot0//{chain_ids[0]}' , '/prot' + str({i}) + '//{chain_ids[i]}')")
        if hide:
            #pm(f"cmd.hide('cartoon','not chain {chain_ids[i]} and prot{str(i)}' )")
            pm(f"cmd.hide('everything','prot{str(i)}')")
            pm(f"cmd.show('cartoon','/prot{str(i)}//{chain_ids[i]} and polymer')")
        #print(i, 'hidden') #hiding is not working

    pm("cmd.zoom('visible')")
    #print('zoomed')
    pm("cmd.center('visible')")
    #print('centered')
    pm("cmd.bg_color('grey80')")
    #print('colored')
    pm(f"cmd.save('{output_file}')")
    #print('saved',  f"cmd.save('{output_file}')")

    #print('done?')



