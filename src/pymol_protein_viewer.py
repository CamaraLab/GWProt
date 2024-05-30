
import os
import pymolPy3
import math
import numpy as np
from subprocess import Popen, PIPE, STDOUT
import pickle
import time
import sys
sys.path.insert(0,'../PGC020.a12/src')
sys.path.insert(0,'../PGC020.a3')
sys.path.insert(0,'../PGC020.a15')

from pymol import cmd
import pymol

import GW_scripts
import read_pdb
import FGW_protein
import IdInit
import GWstress
import weighted_alignment

# copied to src 5/30/2024 from pymol_protein_viewer1



class my_pymolPy3:
    #pretty much the one from online
    def __init__(self):

   #     self.pymolpy3 = Popen(['pymol' ,'-pc'], shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text = True, universal_newlines=True)
        self.pymolpy3 =Popen(["pymol", "-pc"], shell=False, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        
    def __del__(self):

        # Turn off stdin...
        self.pymolpy3.stdin.close()
        # Wait until program is over...
        self.pymolpy3.wait()

        # Terminate the subprcess...
        self.pymolpy3.terminate()

    def __call__(self, s, pause = False):
        # Keep reading the new pymol command as a byte string...
        self.pymolpy3.stdin.write( s + '\n' )

        # Flush it asap...
        self.pymolpy3.stdin.flush()

        if pause:
            self.pymolpy3.stdin.write("with open('lock.tmp', 'w') as f: f.writeline('lock') \n")
            self.pymolpy3.stdin.flush()
            while 'lock.tmp' not in os.listdir('.'):
                time.sleep(0.1)
            os.remove('lock.tmp')
        return 0

    def quit(self):
        self('quit')
        
    def get_err(self):
        return(  self.pymolpy3.communicate())

    def get_data(self, obj:str, tempfile = None):

        if tempfile is None:
            tempfile = str(time.time()) + '.tmp'
        
        #gets the object named obj from the pymol kernel
        if tempfile in os.listdir('.'):
            print(f'{tempfile} already exists, aborting')
            return -1
        
        self('import pickle')
        self(f"with open('{tempfile}', 'wb') as file: pickle.dump({obj}, file)")
        while tempfile not in os.listdir('.') or os.path.getsize(tempfile) == 0:
            time.sleep(0.1)
        with open(tempfile, 'rb') as file:
            local_obj = pickle.load(file)
        os.remove(tempfile)
        return local_obj

    def communicate(self, obj:str):
        # self.pymolpy3.stdin.write(f"print {obj} \n")
        # self.pymolpy3.stdin.flush()
        # str = self.pymolpy3.stdout.read()
        # err = self.pymolpy3.stderr.read()
        str, err = self.pymolpy3.communicate(f"print {obj} \n")
        return str, err
        


def compare_proteins_in_pymol(file1, file2, output_file, threshold = 0.5):
    p1 = FGW_protein.FGW_protein.make_protein_from_pdb(file1)
    p2 =FGW_protein.FGW_protein.make_protein_from_pdb(file2)
    
    c , stress1, stress2 , T= GWstress.GW_stress_from_prots(p1,p2, transport_plan = True)
    ps = GWstress.get_pairing(T, threshold0 = threshold, threshold1 = threshold)
    
    pret, rot, trans = weighted_alignment.weighted_RMSD(p1.coords,p2.coords,T)
    
    ll = weighted_alignment.pymol_transform( pretrans = pret, rot = rot, posttrans = trans)

    cmd.delete('all')
    cmd.load(file1, 'prot1')
    cmd.load(file2, 'prot2')
    #cmd.cealign('prot1', 'prot2')
    residue_numbers1 = list(set(atom.resi_number for atom in cmd.get_model('prot1').atom))
    residue_numbers2 = list(set(atom.resi_number for atom in cmd.get_model('prot2').atom))
    my_namespace = {'new_b1' : stress1, 'new_b2' : stress2, 'res_num1' : residue_numbers1 , 'res_num2' : residue_numbers2} 
    cmd.alter(selection = 'prot1', expression = 'b = new_b1[res_num1.index(int(resi))]', space = my_namespace)
    cmd.alter(selection = 'prot2', expression = 'b = new_b2[res_num2.index(int(resi))]', space = my_namespace)
    cmd.spectrum( expression = "b" , selection = 'prot1', palette = "yellow_red",byres = 1)
    cmd.spectrum( expression = "b" , selection = 'prot2', palette = "blue_red",byres = 1)
    
    
    # Iterate through each pair of residues
    for p in ps:
        i = residue_numbers1[p[0]]
        j = residue_numbers2[p[1]]
           
        # Construct the selection strings for the atoms of interest
        atom_selection1 = f'prot1//A/{i}/CA'
        atom_selection2 = f'prot2//A/{j}/CA'
        
        # Calculate and draw the distance between the selected atoms
        cmd.distance( f'd{i}_{j}' ,selection1=atom_selection1, selection2=atom_selection2)
        cmd.show("dashes",f'd{i}_{j}')
        cmd.hide( 'labels', f'd{i}_{j}')
    cmd.set('dash_gap' , '0')
    cmd.set('dash_width', '1')
    cmd.zoom('all')
    cmd.bg_color('grey80')
    cmd.save('temp.pse')
    pm = my_pymolPy3()
    pm("cmd.load( 'temp.pse') ")
    pm(f"cmd.transform_selection( 'prot1' , matrix =  {ll})")
    pm("cmd.center('all')")
    pm("cmd.zoom('all')")
    pm(f"cmd.save( '{output_file}') ")

def show_proteins_with_values(file,  output_file, *argv):
    #final param is list of data
    cmd.delete('all')
    data = list(argv)
    n = len(data)
    prots = {}
    
    for i in range(n):
        prots[i] = FGW_protein.FGW_protein.make_protein_from_pdb(file1)
        cmd.load(file, 'prot' + str(i))
        if i ==0:
            residue_numbers = list(set(atom.resi_number for atom in cmd.get_model('prot0').atom))
            assert len(residue_numbers) == len(prots[i])
            my_namespace = { 'res_num' : residue_numbers } 
            
        assert len(data[i]) == len(residue_numbers)
        my_namespace['new_b' + str(i) ] = data[i]
        cmd.alter(selection = 'prot' + str(i), expression = 'b = new_b' + str(i)+'[res_num.index(int(resi))]', space = my_namespace)
        cmd.spectrum( expression = "b" , selection = 'prot'+str(i), palette = "yellow_red",byres = 1) 


    cmd.zoom('all')
    cmd.center()
    cmd.bg_color('grey80')
    cmd.save(output_file)
   
    
    


if __name__ == "__main__":
    pymol.finish_launching(['pymol', '-pc'])