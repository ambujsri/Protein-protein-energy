# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 11:31:46 2021

@author: Ambuj
"""


import numpy as np
import re
import pandas as pd
import math
from Bio import PDB
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning


class interaction_energy:
    def __init__(self, pdb_file,prot1_chain='A',prot2_chain='B'):
        self.pdb_file = pdb_file
        self.prot1_chain = prot1_chain
        self.prot2_chain = prot2_chain
        self.pattern ='^ATOM.{16}'
#        self.coords_line = []
#        self.coords = []
#    coords_xyz = np.zeros((3))
# f1_header function returns the header of PDB file.
    def f1_header(self):
        record = open(self.pdb_file,'r').readlines()[0]
        dict1 = {'title':record[10:50].rstrip(),'date':record[51:60].rstrip(),'ID':record[62:66]}
        return dict1
# f2_atom_coords return the ATOM records from the PDB file for the desired chain
    def f2_atom_coords(self, chain='A'):
        coords_line = []
        if len(chain) == 1:
            chain = '.'+chain
#            print chain
            pdb_ptr = open(self.pdb_file,'r')
            
            for line in pdb_ptr:
                
                pat_srch = re.search(self.pattern+chain, line)
#                chx = line[21:23]
#                print chx
#                break
                if pat_srch or line[21:22] == 'chain':
#                    print line
                    coords_line.append(line)
#                    self.coords.append(list([line[31:38], line[39:46], line[47:54]]))
#                    self.coords[i].append(float(line[39:46]))
#                    self.coords[i].append(float(line[47:54]))
#            print self.coords
        return coords_line
# f3_extract_seqres extract the SEQRES records from the PDB file 
    def f3_extract_seqres(self, chain='A'):
        seq_line = []
        pattern2 = '^SEQRES.{4}'
#        pattern3 = '^
        if len(chain) == 1:
            chain = '.'+chain
            pdb_ptr = open(self.pdb_file,'r')
            for line in pdb_ptr:
                pat_srch = re.search(pattern2+chain, line)
                if pat_srch:
                    seq_line.append(line[19:])
                    
        #            spt_line = line[19:].rstrip().split()
                    
#        check = molecule_type(pdb_file, chain)        
#        seq=inst.fasta_format(seq_line)
        return seq_line
# f4_seqres_to_fasta is used to extract sequence in fasta format from SEQRES records of PDB file 
    #Output of f3_extract_seqres is input for this function
    def f4_seqres_to_fasta(self, seq_line, header):    
        aa = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M','PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
        nt = {'DA':'A','DG':'G','DU':'U','DC':'C','DT':'T', 'DI':'I'}
        rnt = ['A','G','C','U','I']
        newseq = ''
        finalseq =[]
        finalseq.append('>'+header['ID'])
        for line in seq_line:
            spt_line = line.rstrip().split()
            for res in spt_line:
                res = res.upper()
                if len(res) == 3:
                    if res in aa:
                        newseq += aa[res]
                    else:
                        print (res+' is not a standard amino acid residue')
                if len(res) == 2:
                    if res in nt:
                        newseq += nt[res]
                    else:
                        print (res+' is not a standard DNA nucleotide')
                if len(res) == 1:
                    if res in rnt:
                        newseq += res;
                    else: 
                        print (res+' is not a standard RNA nucleotide')
        finalseq.append(newseq)
        return finalseq
# f5_seq_comp function calculate the sequence composition from SEQRES record
#Output of f3_extract_seqres is input for this function
    def f5_seq_comp(self, seq_line):          
        aa = {'ALA':0,'ARG':0,'ASN':0,'ASP':0, 'CYS':0, 'GLN':0, 'GLU':0, 'GLY':0, 'HIS':0, 'ILE':0, 'LEU':0, 'LYS':0, 'MET':0,'PHE':0, 'PRO':0, 'SER':0, 'THR':0, 'TRP':0, 'TYR':0, 'VAL':0}
        dnt = {'DA':0,'DG':0,'DU':0,'DC':0,'DT':0, 'DI':0}
        nt = {'A': 0, 'G':0, 'C':0, 'U':0, 'I':0}
        x1 = 0
        x2 = 0
        x3 = 0
        
        for line in seq_line:
            spt_line = line.rstrip().split()
#            print (spt_line)
            for res in spt_line:
                res = res.upper()
#                print len(res)
                if len(res) == 3:
                    if res in aa:
                        aa[res] += 1
                        x1 += 1
                if len(res) == 2:
                    if res in dnt:
                        dnt[res] += 1
                        x2 = x2 + 1
#                        print x2
                if len(res) == 1:
                    nt[res] += 1
                    x3 += 1
                    
#        print dnt
        if x1 >= x2 and x1 >= x3:
            return aa
        elif x2 >= x1 and x2 >= x3:
            return dnt
        elif x3 >= x1 and x3 >= x2:
            return nt

# f6_protein_protein_contact2 calculate the atom to atom cotact under the distance cutoff (default is 5 Angstrom)
# The PDF file path, protein chain and RNA chain will be used for the calculation which are initialized during variable declaration
    
    def f6_protein_protein_contact2(self, dist_cutoff=5):
        
        warnings.simplefilter('ignore', PDBConstructionWarning)
        int_df1 = pd.DataFrame()
        flag = 0
        parser = PDB.PDBParser()
        structure = parser.get_structure("pdb", self.pdb_file)
        model = structure[0]
        prot1_chain = model[self.prot1_chain] 
        prot2_chain = model[self.prot2_chain]
        # nal1 = ['A','G','C','U']
        pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
        # c1 = [str(x).split()[1] for x in list(prot1_chain.get_residues())]
        # c2 = [str(x).split()[1] for x in list(prot2_chain.get_residues())]
        # checkp = set(c1).intersection(set(pal1))
        # checkn = set(c2).intersection(set(pal1))
        # if (len(checkp) > 0):
        #     print("check the chains of pdb file "+self.pdb_file+" error in protein1 or protein2 molecule\n See protein1 set "+str(checkp))
        #     flag = 1
        #     return int_df1
        # if (len(checkn) >= 2):
        #     print("check the chains of pdb file "+self.pdb_file+" error in protein1 or protein2 molecule\n See protein2 set "+(str(checkn)))
        #     flag = 1
        #     return int_df1
        for prot1_res in prot1_chain:
            for prot1_atoms in prot1_res:
                for prot2_res in prot2_chain:
                    prot2_resname = prot2_res.resname
                    for prot2_atoms in prot2_res:
                        distance = prot1_atoms-prot2_atoms
                        # print (distance)
                        if (distance< dist_cutoff):
                            dict1 = {'distance':distance, 'prot2_atm':prot2_atoms.get_full_id()[4][0], 'prot2_atmno':prot2_atoms.get_serial_number(),
                                     'prot2_res':prot2_res.resname, 'prot2_resno':prot2_atoms.get_full_id()[3][1], 'prot2_coord':prot2_atoms.get_coord(), 'prot1_atm':prot1_atoms.get_full_id()[4][0], 
                                     'prot1_atmno':prot1_atoms.get_serial_number(), 'prot1_res':prot1_res.resname, 'prot1_resno':prot1_atoms.get_full_id()[3][1], 
                                     'prot1_coord':prot1_atoms.get_coord() }
                            
                            int_df1 = int_df1.append(dict1,ignore_index=True)
        
        b1 = [x.strip() in pal1 for x in int_df1['prot1_res']]
        temp1 = int_df1[b1]
        b2 = [x.strip() in pal1 for x in temp1['prot2_res']]
        df_inter = temp1[b2]
        return df_inter
    
# f7_protein_energy2 is used to calculate the energy value using optimized van der Waals parameters for amino acids and RNA nucleotides
# Output of f6_proteinRNAcontact2 is to be used as df_inter for this function
# aa_param and na_param are the parameters optimized for amino acids and RNA nucleotides, respectively
    def f7_protein_energy2(self, df_inter, aa_param):     
                                                            
        df1 = df_inter.copy(deep=True)
        
        
        vdwrad = 'Vdwradrna'
        vdweps = 'Vdwepsrna'
        
        total_energy = 0
        df1['prot1_atmtype'] ='NA'
        df1['prot1_charge'] = 0
        df1['prot1_vdwradius'] =0
        df1['prot1_vdweps'] =0
        df1['prot2_atmtype'] ='NA'
        df1['prot2_charge'] =0
        df1['prot2_vdwradius'] =0
        df1['prot2_vdweps'] =0
        df1['Vdw_energy'] =0
        
#        numcols = len(df1.columns)
        for i in df1.index:
            idx1 = aa_param[(str(df1['prot1_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot1_res'][i]).strip() == aa_param['Res_name'])].index
            # print (idx1)
            if (len(idx1) == 0 ):
                print ("Unavailable in protein parameter: Residue - "+ str(df1['prot1_res'][i]).strip()+ " , Atom - "+ str(df1['prot1_atm'][i]).strip())
            else:
                df1.loc[i,'prot1_atmtype'] = list(aa_param['Atm_type'][idx1])[0]
                df1.loc[i,'prot1_charge'] = list(aa_param['Charge'][idx1])[0]
                df1.loc[i,'prot1_vdwradius'] = list(aa_param['Vdwrad'][idx1])[0]
                df1.loc[i,'prot1_vdweps'] = list(aa_param['Vdweps'][idx1])[0]
            idx2 = aa_param[(str(df1['prot2_atm'][i]).strip() == aa_param['Atm_name']) & (str(df1['prot2_res'][i]).strip() == aa_param['Res_name'])].index
            if (len(idx2) == 0 ):
                # print (list(df1.columns))
                print ("Unavailable in RNA parameter: Residue - "+ str(df1['prot2_res'][i]).strip()+ " , Atom - "+ str(df1['prot2_atm'][i]).strip())
            else:
                df1.loc[i,'prot2_atmtype'] = list(aa_param['Atm_type'][idx2])[0]
                df1.loc[i,'prot2_charge'] = list(aa_param['Charge'][idx2])[0]
                df1.loc[i,'prot2_vdwradius'] = list(aa_param['Vdwrad'][idx2])[0]
                df1.loc[i,'prot2_vdweps'] = list(aa_param['Vdweps'][idx2])[0]
            eps = float(df1.loc[i,'prot1_vdweps']) *float(df1.loc[i,'prot2_vdweps'])
            vdwr = float(df1.loc[i,'prot1_vdwradius']) +float(df1.loc[i,'prot2_vdwradius'])
            chrg = float(df1.loc[i,'prot1_charge']) *float(df1.loc[i,'prot2_charge'])
            vdwr2 = vdwr*vdwr
            vdwr6 = vdwr2*vdwr2*vdwr2
            vdwr12 = vdwr6*vdwr6
            epssqrt = math.sqrt(eps)
            
            aec12ab = epssqrt*vdwr12
            aec6ab = 2*epssqrt*vdwr6
            rijs = float(df1.loc[i,'distance'])*float(df1.loc[i,'distance'])
            rs = float(1)/rijs
            rij2=rs
            
            rij6 = rij2*rij2*rij2
            rij12 = rij6*rij6
            
            
            # print (rij12, aec12ab)
            v1 = (aec12ab*rij12)-(aec6ab*rij6)
            v2 = chrg*rij2
            v12 =v1+v2
            # print (v1,v2,v12)
            total_energy += v12
    #            print v12
            if (v12 < 0):
                df1.loc[i,'Vdw_energy'] = v12
            else:
                df1.loc[i,'Vdw_energy'] = 0
        
    
#            v1 = (aec12ab*rij12)-(aec6ab*rij6)
#            v2 = chrg*rij2
#            v12 =v1+v2
        return df1
# f8_interaction_type function calculate the composition of interacting atoms (e.g. 'CO' represent a 'C' atom from protein interacting with 'O' atom of RNA) 
# Output of f6_proteinRNAcontact2 is input for this function
    def f8_interaction_type (self, df_inter):                     
        int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
        temp_ptr = open('temp.txt','w')
        for index, dfrow in df_inter.iterrows():
            temp_ptr.write('\n'+dfrow["prot1_atm"]+'\t'+dfrow["prot2_atm"])
            in_t = dfrow["prot1_atm"][0:1]+dfrow["prot2_atm"][0:1]
            if in_t in int_dict:
                int_dict[in_t] += 1
            else:
                int_dict.update({in_t:1})
        temp_ptr.close()
        return int_dict




aa_param = pd.read_csv('aa_20vdrch.csv')
# na_param = pd.read_csv('rna_4vdrch.csv')

inst2 = interaction_energy('6wpt.pdb',prot1_chain='A',prot2_chain='C')
h1 = inst2.f1_header()
atom_coords = inst2.f2_atom_coords(chain='A')
seq_line=inst2.f3_extract_seqres()
seq = inst2.f4_seqres_to_fasta(seq_line,h1)
seq_comp = inst2.f5_seq_comp(seq_line)
df_inter = inst2.f6_protein_protein_contact2(dist_cutoff=3.5)
df_inter.to_csv('atom_interaction.csv')
energy_df = inst2.f7_protein_energy2(df_inter,aa_param)
energy_df.to_csv('interaction_energy.csv')
atoms_involve = inst2.f8_interaction_type(df_inter)



