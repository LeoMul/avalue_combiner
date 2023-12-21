#ultimate a value combiner and what not.
#rn - reads in from files e1, e2, e3 etc. howver will eventually make it grasp.outs
import numpy as np
num_levels = 40
path1 = 'GRASP.OUT'
transition_matrix = np.zeros([num_levels,num_levels])
input_array = [path1]


def add_transition_matrix(transition_matrix_old,new_data):
    #new data = from readin functions
    #transition_matrix_old = currently read in and processed data.

    lower = new_data[:,0].astype(int) - 1 #subtracting as python is 0 indexed, god save us
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

def output(transition_matrix,output_file,num_levels):
    #f = open(output_file,'w')
    print("Outputing:")
    print("Upper Lower A Value")
    for ii in range(0,num_levels):
        for jj in range(ii+1,num_levels):
            lower = ii+1#python is 0 indexed. god save us.
            upper = jj+1
            avalue = transition_matrix[ii,jj]
            #print(upper,lower,avalue)
            string_to_be_written = '{:5}'.format(upper)
            string_to_be_written += '{:5}'.format(lower)
            avalue_to_be_written = round(avalue,ndigits=2)
            string_to_be_written += '{:}'.format(" ")
            avalue_string = '{:.2E}'.format(avalue)
            avalue_string = avalue_string.replace("E","")
            string_to_be_written += avalue_string
            print(string_to_be_written)


    return 0
def readin_until_finds_aij(file):
    
    #this function reads in a data file of path 'file' and outputs any electric AND
    #magnetic transitions it finds.
    #it reads in the file iwht the regular python api and locates the positions of the data.
    #these positions are then re-read in by numpy with the data straight into ndarrays.
    #while this isn't particularly elegant - it is the least amount of work and less error prone 
    #than converting the data myself.

    f = open(file,'r')
    content = f.readlines()
    numlines = len(content)
    print("Length of file ", file, "is ",numlines)

    for i in range(0,numlines):
        current_line = content[i]
        current_line_list = current_line.split()
        if current_line_list != []: 
            if current_line_list[0] == "Electric":
                print("found electric")
                print(current_line)
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
                print("found magnetic")
                print(current_line)
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
    return electric_data,magnetic_data


def main():

    for file in input_array:

        electric,magnetic = readin_until_finds_aij(file)
        transition_matrix = add_transition_matrix(transition_matrix,electric)
        transition_matrix = add_transition_matrix(transition_matrix,magnetic)

    transition_matrix = finalise(transition_matrix)

    output(transition_matrix,'we',num_levels=10) 

main()