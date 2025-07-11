
for p in ps:
    i = str(residue_numbers1[p[0]])
    j = str(residue_numbers2[p[1]])
    atom_selection1 = 'prot1//'+chain1+'/'+i+'/CA'
    atom_selection2 = 'prot2//'+chain2+'/'+j+'/CA'
    cmd.distance('d' + i + '_' + j, selection1=atom_selection1, selection2=atom_selection2)
    cmd.show('dashes', 'd'+i +'_'+ j)
    cmd.hide('labels', 'd'+i +'_'+ j)
    #print(f"cmd.hide('labels', 'd{i}_{j}')")
	#cmd.color_deep('blue', 'd'+i +'_'+ j,0)