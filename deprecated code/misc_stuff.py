"""
This file is for storing code we don't need but I didn't want to totally delete



"""


# This is just a method for comparing the earth mover distances stresses of different structural alignments of proteins

def get_pymol_transport(file1, file2, aligner = 'cealign', TM_exe = '/home/elijah/pymol/bin/TMalign' ): #filepaths to pbds
    # TM_exe filepath to TM align
    cmd.delete('all')
    cmd.load(file1, 'prot1')
    cmd.load(file2, 'prot2')

    coords01 = []
    coords02 = []
    cmd.iterate_state(1, "prot1 and name CA", "coords01.append((x, y, z))", space={'coords01': coords01})
    cmd.iterate_state(1, "prot2 and name CA", "coords02.append((x, y, z))", space={'coords02': coords02})
    coords01 = np.stack(coords01)
    coords02 = np.stack(coords02)

    match aligner:
        case 'cealign':
            cmd.cealign('prot1', 'prot2')
        case 'align':
            cmd.align('prot1', 'prot2')
        case 'super':
            cmd.super('prot1', 'prot2')     
        case 'tmalign':
            pymol_tmalign_wrapper_Copy1.tmalign('prot1', 'prot2', quiet=1, exe = TM_exe, return_alignment=False) 
            #cmd.tmalign('prot1', 'prot2') 
        case _:
            raise ValueError("valid arguments for aligner are 'align', 'cealign', 'super', tmalign")
    coords1 = []
    coords2 = []
    cmd.iterate_state(1, "prot1 and name CA", "coords1.append((x, y, z))", space={'coords1': coords1})
    cmd.iterate_state(1, "prot2 and name CA", "coords2.append((x, y, z))", space={'coords2': coords2})
    coords1 = np.stack(coords1)
    coords2 = np.stack(coords2)
    
    
    #print(coords1.shape)
    #print(coords2.shape)
    D = ot.dist(coords1, coords2)
    #print(D.shape)
    a = np.ones(coords1.shape[0])/coords1.shape[0]
    b = np.ones(coords2.shape[0])/coords2.shape[0]

    T = ot.emd(a,b, D)
    stress = np.einsum('ij,ij ->ij', D,T)
    cost = np.sum(stress)
    stress1 = np.sum(stress, axis = 1)
    stress2 = np.sum(stress, axis = 0)

    return T, cost, stress1, stress2 

   #warning -  the pymol coords could be returned in a different order
# or try cmd.get_coords(str selection)
# not totally sure how ordering works

# pymol_tmalign_wrapper_Copy1.tmalign('prot1', 'prot2', quiet=1, return_alignment=True)



# 1 get distance matrix

def dist_csv_to_dist_matrix_no_doubles(
    infile: str, #name of csv file 
    pdb_file_list: list[str], #list of all the pdb files we should expect, matrix is output in that order
    doubled_ids: list[str] = [], # list of ids we should remove and exclude, may or may not be in pdb_file_list
    check: bool = True,
    fun = lambda x: x, #applying a function to the data to get distances
    pdb_files = True #whether it assumes the file names are .pdb and thus drops all other lines
    ): 
    
    #print("test")
    #this function also will deal with cases where the csv has distances between a protein and itself
    
    assert len(doubled_ids) == len(set(doubled_ids))
    for id in doubled_ids:
        if id in pdb_file_list:
            pdb_file_list.remove(id)

    entry_dict = {}
    # protein: index in pdb_file_list
    for i in range(len(pdb_file_list)):
        entry_dict[pdb_file_list[i]]=i

    assert len(pdb_file_list) == len(entry_dict.keys())
    assert set(range(len(pdb_file_list))) == set(entry_dict.values())

    dist_mat = np.zeros((len(pdb_file_list), len(pdb_file_list)))


    
    counter = 1
    doubled_prots = set()#for debugging
    with open(infile, "r", newline="") as infile:
        csv_reader = csv.reader(infile, delimiter=",")

        for line in csv_reader:
            if '.pdb' not in line[0] and pdb_files:
                continue            
            if line[0] not in pdb_file_list or line[1] not in pdb_file_list:
                continue    
            if line[0] in doubled_ids or line[1] in doubled_ids:
                continue
            if line[0] == line[1]:
                doubled_prots.add(line[0])
            counter +=1
            if float(line[2]) <= 0:
                dist_mat[entry_dict[line[0]]][entry_dict[line[1]]] = 0

            else:
                
                dist_mat[entry_dict[line[0]]][entry_dict[line[1]]] = max(fun(float(line[2])),0)
                dist_mat[entry_dict[line[1]]][entry_dict[line[0]]] = max(fun(float(line[2])),0)

    files = set(pdb_file_list)       
    #print(files.difference(doubled_prots))
    #print(doubled_prots.difference(files))
#     print(counter)
#     print(len(pdb_file_list))
#     print(len(list(doubled_prots)))
    if check:
        assert 2* counter == len(pdb_file_list)**2 - len(pdb_file_list)
        assert scipy.spatial.distance.is_valid_dm(dist_mat)
    
    #assert len(dist_mat) ==5127
    
    return dist_mat
    


def GW_from_coords(prot1_coords, prot2_coords, n = np.inf, scaler = lambda x : x):
    # assumes uniform distributions
    # prot1_coords - (n1,3) array or nested list of the CA coords
    # prot2_coords - (n2,3) array or nested list of the CA coords
    # n - downsampling 

    new_indices1=np.linspace(0, len(prot1_coords), num=min(n, len(prot1_coords)), endpoint=False, dtype=np.int_) # samples n of them, evenly spaced
    new_indices2=np.linspace(0, len(prot2_coords), num=min(n, len(prot2_coords)), endpoint=False, dtype=np.int_) # samples n of them, evenly spaced

    downsampled_p1 = prot1_coords[new_indices1, :]
    downsampled_p2 = prot2_coords[new_indices2, :]

    ipdm1 = squareform(pdist(downsampled_p1))
    ipdm2 = squareform(pdist(downsampled_p2))


    GW_cell1 = gw_cython.GW_cell( dmat = np.vectorize(scaler)(ipdm1), distribution = unif(min(n, len(prot1_coords))))
    GW_cell2 = gw_cython.GW_cell( dmat = np.vectorize(scaler)(ipdm2), distribution = unif(min(n, len(prot2_coords))))


    return GW_identity_init(GW_cell1, GW_cell2)



def GW_from_ipdms(ipdm1, ipdm2,n=np.inf, scaler = lambda x : x):

    n1 = len(ipdm1)
    if n1 > n:
        ipdm1 = submat(ipdm1, np.linspace(0, n1, num=n, endpoint=False, dtype=np.int_))

    n2 = len(ipdm2)
    if n2 > n:
        ipdm2 = submat(ipdm2, np.linspace(0, n2, num=n, endpoint=False, dtype=np.int_))

    ipdm1 = np.vectorize(scaler)(ipdm1)
    ipdm2 = np.vectorize(scaler)(ipdm2)

    GW_cell1 = gw_cython.GW_cell( dmat = (ipdm1), distribution = unif(len(ipdm1)))

    GW_cell2 = gw_cython.GW_cell( dmat = (ipdm2), distribution = unif(len(ipdm2)))


    return GW_identity_init(GW_cell1, GW_cell2)

def random_permutation_initial_coupling_unif(m,n, seed = None):
    if seed:
        np.random.seed(seed)
    P = id_initial_coupling_unif(m,n)
    Q = np.random.permutation(P)
    return(Q)
    
