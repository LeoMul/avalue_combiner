#ultimate a value combiner and what not.
#tested, kind of. I made an effort to make it idiot proof.
#However - since I, Leo, AM an idiot - this leaves the user liable to check what they are doing regardless.

import numpy as np
import argparse

def add_transition_matrix(transition_matrix_old,new_data):
    #new data = from readin functions
    #transition_matrix_old = currently read in and processed data.

    lower = new_data[:,0].astype(int) - 1 #subtracting as python is 0 indexed, as God intended it
    upper = new_data[:,1].astype(int) - 1
    avalues = new_data[:,3]

    for kk in range(0,len(avalues)):
        i = lower[kk]
        j = upper[kk] 
        avalue = avalues[kk]
        transition_matrix_old[i,j]+= avalue
    return transition_matrix_old

def finalise(transition_matrix):
    si = np.shape(transition_matrix)
    for jj in range(0,si[0]):
        for ii in range(0,si[1]):
            if transition_matrix[jj,ii] == 0.0:
                transition_matrix[jj,ii] = 1.0e-30
    return transition_matrix

def reorder(transition_matrix,new_order_file):
    new_transition_matrix = np.zeros_like(transition_matrix)
    skip = 0
    #read the dord-esque file
    dord = open(new_order_file,'r')
    check_string_array = dord.readline().split()
    if check_string_array[0] == "&OMORDER":
        print("dord format detected, skipping namelist")
        skip = 1
    dord.close()
    new_order_data = np.loadtxt(new_order_file,skiprows=skip)

    old_levels = new_order_data[:,0].astype(int)-1
    new_levels = new_order_data[:,1].astype(int)-1
    #print(old_levels)
    num = len(old_levels)
    expected = np.shape(transition_matrix)[0]
    if (num != expected):
        print("INCOMPLETE REORDERING INDEX.")
        print("Expected ",expected," found ",num)
        exit()
    for jj in range(0,num):
        for ii in range(0,num):
            old_i = old_levels[ii]
            old_j = old_levels[jj]
            new_i = new_levels[ii]
            new_j = new_levels[jj]
            new_transition_matrix[new_i,new_j] = transition_matrix[old_i,old_j]

    return new_transition_matrix

def output(transition_matrix,output_file,num_levels):
    f = open(output_file,'w')
    print("Writing to ",output_file)
    f.write("#Upper Lower A Value\n")
    for ii in range(0,num_levels):
        for jj in range(ii+1,num_levels):
            lower = ii+1#python is 0 indexed, but we are 1-indexed
            upper = jj+1
            avalue = transition_matrix[ii,jj]
            #this clunky formatting is to make it the same as the grasprad code
            #i could use fortranformat to make it exactly the same,
            #but this is a dependency not many will have, or be bothered to intsall.
            string_to_be_written = '{:5}'.format(upper)
            string_to_be_written += '{:5}'.format(lower)
            avalue_to_be_written = round(avalue,ndigits=2)
            string_to_be_written += '{:}'.format(" ")
            avalue_string = '{:.2E}'.format(avalue)
            avalue_string = avalue_string.replace("E","")
            string_to_be_written += avalue_string +"\n"
            f.write(string_to_be_written    )


    return 0
def readin_until_finds_aij(file):
    
    #this function reads in a data file of path 'file' and outputs any electric AND
    #magnetic transitions it finds.
    #it reads in the file iwht the regular python api and locates the positions of the data.
    #these positions are then re-read in by numpy with the data straight into ndarrays.
    #while this isn't particularly elegant - it is the least amount of work and less error prone 
    #than converting the data myself.
    
    #that said - this is already pretty bug prone as I repeat pretty much the same code twice.
    #for electric and magnetic. possibly refactor a wee bit.
    #it works anyway i guess

    f = open(file,'r')
    content = f.readlines()
    numlines = len(content)
    print("Length of file ", file, "is ",numlines)
    found_array = []
    for i in range(0,numlines):
        current_line = content[i]
        current_line_list = current_line.split()
        if current_line_list != []: 
            if current_line_list[0] == "Electric":
                current_line = current_line.replace("\n","")
                print("found electric")
                print(current_line)
                found_array.append(current_line)
                break
    if i >= (numlines-1): 
        print("no Electric transitions found, check input")
        return
    first_line_with_aij_index_electric = i + 12
    last_line_electric = first_line_with_aij_index_electric
    for j in range(first_line_with_aij_index_electric,numlines):
        if content[j].split() != []:
            last_line_electric+=1 
        else:
            break 
    #print(last_line_electric)

    print("Electric data begins at ",first_line_with_aij_index_electric," and ends at ",last_line_electric)
    
    for i in range(last_line_electric,numlines):
        current_line = content[i]
        current_line_list = current_line.split()
        if current_line_list != []: 
            if current_line_list[0] == "Magnetic":
                current_line = current_line.replace("\n","")
                print("found magnetic")
                print(current_line)
                found_array.append(current_line)
                break
    if i >= (numlines-1): 
        print("no magnetic transitions found, check input")
        return 
    first_line_with_aij_index_magnetic = i + 12
    last_line_magnetic = first_line_with_aij_index_magnetic
    for j in range(first_line_with_aij_index_magnetic,numlines):
        if content[j].split() != []:
            last_line_magnetic+=1 
        else:
            break 
    #print(last_line_magnetic)

    print("Magnetic data begins at ",first_line_with_aij_index_magnetic," and ends at ",last_line_magnetic)
    f.close()

    electric_data = np.loadtxt(file,usecols=[0,1,2,3],skiprows=first_line_with_aij_index_electric,max_rows=last_line_electric-first_line_with_aij_index_electric)
    magnetic_data = np.loadtxt(file,usecols=[0,1,2,3],skiprows=first_line_with_aij_index_magnetic,max_rows=last_line_magnetic-first_line_with_aij_index_magnetic)
    num_levels_1 = max(electric_data[:,1])
    num_levels_2 = max(magnetic_data[:,1])
    #print(num_levels_1,num_levels_2)
    return electric_data,magnetic_data,found_array,int(max(num_levels_1,num_levels_2))

#def globbingfiles(globstring):
#    files = []
#    for file in glob.glob(globstring):
#        files.append(file)
#    return files

def manualdatafiles(datafiles):
    files = []
    for file in args.datafiles:
        files.append(file)
    return files  



#path1 = 'GRASP.OUT'
parser = argparse.ArgumentParser()
# Adding optional argument
parser.add_argument('-n', '--name',  help='Specify desired file name. Leave blank for default: MERGED_A_VALUES.OUT')
#parser.add_argument('-g', '--globbing',  help='Selects datafiles by globbing')
parser.add_argument('-d', '--datafiles', nargs='+', help='Paths of input files. Multiple files put a space between them.')
parser.add_argument('-r', '--reorder', help='Path of reorder file - i.e translates unshifted state indices to their shifted ')
args = parser.parse_args()
#num_levels = 40

def main(output_file_name):
    found_array = []
    found_num_levels = []
    obtained_transition_arrays = []


    for file in input_array:
        electric,magnetic,found,num_levels = readin_until_finds_aij(file)
        found_num_levels.append(num_levels)
        obtained_transition_arrays.append(electric)
        obtained_transition_arrays.append(magnetic)
        found_array += (found)
    num_levels = max(found_num_levels)
    print("I think the number of levels is",num_levels)

    transition_matrix = np.zeros([num_levels,num_levels])
    for matrix in obtained_transition_arrays:
        transition_matrix = add_transition_matrix(transition_matrix,matrix)

    transition_matrix = finalise(transition_matrix)
    
    
    if len(found_array) != len(set(found_array)):
        print("DUPLICATE TRANSITIONS FOUND CHECK INPUT")
        print("Located transitions: ")
        print(found_array)
        print("aborting")
    else:
        if args.reorder:
            print("REORDERING LEVELS ACCORDING TO INPUT FILE ",args.reorder)
            transition_matrix = reorder(transition_matrix,args.reorder)

        output(transition_matrix,output_file_name,num_levels=num_levels) 



if args.datafiles:
    input_array = manualdatafiles(args.datafiles) 
    output_file_name = args.name if args.name else "MERGED_A_VALUES.OUT"
    main(output_file_name)
else:
    print("no data files specified - stopping and printing help")
    parser.print_help()
