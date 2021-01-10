import re
import numpy as np
import os
import sys

def make_symmetric(matrix):
    return matrix+matrix.T - np.diag(np.diag(matrix))

def triangle_fill(matrix,istrt,iend,jstrt,count,element):
    for i in range(istrt,iend):
        for j in range(jstrt,i+1):
            matrix[i,j] = element[count]
            count += 1
    return matrix, count

def block_fill(matrix,istrt,nbf,jstrt,jend,count,element):
    for i in range(istrt,nbf):
        for j in range(jstrt,jend):
            matrix[i,j] = element[count]
            count += 1
    return matrix, count

def create_matrix(elements, nbf):
  count  = 0 	# count is our index
  matrix = np.zeros((nbf,nbf))
  for i in range(0,nbf-nbf%5,5):
    matrix,count = triangle_fill(matrix,i,i+5,i,count,elements)
    matrix,count = block_fill(matrix,i+5,nbf,i,i+5,count,elements)
  matrix,count = triangle_fill(matrix,nbf-nbf%5,nbf,nbf-nbf%5,count,elements)
  return make_symmetric(matrix)

def main(g09file):
    file_name, file_extension = os.path.splitext(g09file)
    if os.path.exists(file_name):
        sys.exit('Error: directory exists! Delete it and re-run gauss_parse')
    else:
        os.makedirs(file_name)
    np.set_printoptions(precision=2) 

    done_get_overlap = False
    done_get_ke      = False
    done_get_pe      = False
    done_get_eri     = False

    logfile = open(g09file,'r')
    for text in logfile:
        words = text.split()
        if all(x in words for x in ['Overlap']):
            done_get_overlap = True
        if all(x in words for x in ['Kinetic', 'Energy']):
            done_get_ke = True
        if all(x in words for x in ['Potential', 'Energy']):
            done_get_pe = True
        if all(x in words for x in ['Dumping','Two-Electron','integrals']):
            done_get_eri = True
        if all(x in words for x in ['primitive','gaussians,','basis','functions,']):
            nbf = int(words[0])   # number basis functions 
        if all(x in words for x in ['nuclear','repulsion','energy','Hartrees.']):
            enuc = float(words[3]) # nuclear repulstion energy in Hartrees
        if all(x in words for x in ['alpha','beta','electrons']):
            nelec = [int(words[0]), int(words[3])] # number alpha elec and beta elec 
    logfile.close()

    logfile = open(g09file,'r')
    data    = logfile.read()
    # grab all text between  "Overlap ***" and "*** Kinetic"
    raw_overlap_string = re.findall(r'Overlap \*\*\*(.*?)\*\*\* Kinetic',data,re.DOTALL)
    raw_overlap_string = raw_overlap_string[0].replace('D','E')
    raw_overlap_elements = raw_overlap_string.split()
    matrix_elements = []
    for overlap_value in raw_overlap_elements:
        if 'E' in overlap_value:
            matrix_elements.append(overlap_value)
    overlap = create_matrix(matrix_elements, nbf)
    logfile.close()


    logfile = open(g09file,'r')
    data    = logfile.read()
    raw_ke_string = re.findall(r'Kinetic Energy \*\*\*(.*?)Entering OneElI...',data,re.DOTALL)
    raw_ke_string = raw_ke_string[0].replace('D','E')
    raw_ke_elements = raw_ke_string.split()
    matrix_elements = []
    for ke_value in raw_ke_elements:
        if 'E' in ke_value:
            matrix_elements.append(ke_value)
    ke_matrix = create_matrix(matrix_elements, nbf)
    logfile.close()

    logfile = open(g09file,'r')
    data = logfile.read()
    raw_pe_string = re.findall(r'Potential Energy \*\*\*\*\*(.*?)\*\*\*\*\*\* Core Hamiltonian',data,re.DOTALL)
    raw_pe_string = raw_pe_string[0].replace('D','E')
    raw_pe_elements = raw_pe_string.split()
    matrix_elements = []
    for pe_value in raw_pe_elements:
        if 'E' in pe_value:
            matrix_elements.append(pe_value)
    pe_matrix = create_matrix(matrix_elements,nbf)
    logfile.close()

    logfile = open(g09file,'r')
    eri_val_list = []
    eri_index_list = []
    count = 0
    for text in logfile:
        words = text.split()
        if 'I=' and 'J=' and 'K=' and 'L=' in words:
            eri_index_list.append([int(words[1]),int(words[3]),int(words[5]),int(words[7])])
            eri_val_list.append([float(words[9].replace('D','E'))])
    eri_val_array   = np.array(eri_val_list).T
    eri_index_array = np.array(eri_index_list)
    logfile.close()

    if done_get_overlap:
        np.savetxt(file_name + '/overlap.dat',overlap,fmt='% 20.12e',delimiter = ',')
    if done_get_ke:
        np.savetxt(file_name + '/kinetic_energy.dat',ke_matrix,fmt='% 20.12e',delimiter = ',')
    if done_get_pe:
        np.savetxt(file_name + '/potential_energy.dat',pe_matrix,fmt='% 20.12e',delimiter = ',')
    if done_get_eri:
        np.savetxt(file_name + '/eri_val.dat',eri_val_array.T,fmt='% 20.12e',delimiter = ',')
        np.savetxt(file_name + '/eri_index.dat',eri_index_array,fmt='%d, %d, %d, %d',delimiter = ',')
        
    np.savetxt(file_name + '/nuclear_repulsion.dat',np.array([enuc]),fmt='% 20.12f')
    np.savetxt(file_name + '/number_basis_functions.dat',np.array([nbf]),fmt='%d')
    np.savetxt(file_name + '/number_electrons.dat',np.array([nelec]),fmt='%d', delimiter = ',')

if __name__ == '__main__':
    main(sys.argv[1])
