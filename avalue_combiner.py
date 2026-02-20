#ultimate a value combiner and what not.
#tested, kind of. I made an effort to make it idiot proof.
#However - since I, Leo, AM an idiot - this leaves the user liable to check what they are doing regardless.

import numpy as np
import argparse

from vr_lib import * 


elements=['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg'
        ,'Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn'
        ,'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr'
        ,'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb'
        ,'Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd'
        ,'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir'
        ,'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th'
        ,'Pa','U ']


AXELROD_CONST_FOR_NOW = 0.004 
RYDBERG_CONSTANT_CM = 109737.316


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

    transition_matrix = transition_matrix + np.transpose(transition_matrix)

    for jj in range(0,si[0]):
        for ii in range(0,si[1]):
            if transition_matrix[jj,ii] == 0.0:
                transition_matrix[jj,ii] = 1.0e-30
    return transition_matrix

def reorder(transition_matrix,new_order_file,dd1=np.array([]),dd2=np.array([])):
    new_transition_matrix = np.zeros_like(transition_matrix)
    skip = 0
    #read the dord-esque file
    if len(dd1) == 0:
        dord = open(new_order_file,'r')
        check_string_array = dord.readline().split()
        if check_string_array[0] == "&OMORDER":
            print("dord format detected, skipping namelist")
            skip = 1
        dord.close()
        new_order_data = np.loadtxt(new_order_file,skiprows=skip)
        old_levels = new_order_data[:,0].astype(int)-1
        new_levels = new_order_data[:,1].astype(int)-1
    else:
        print('will need to use old levels new levels here')
        old_levels = dd1 
        new_levels = dd2 
        print(new_levels[71],old_levels[71])
        
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
            #print(ii,jj,new_i,new_j)
            new_transition_matrix[new_i,new_j] = transition_matrix[old_i,old_j]

    return new_transition_matrix

def output(transition_matrix,output_file,num_levels,exp_energies):
    f = open(output_file,'w')
    print("Writing to ",output_file)

    header = "#Upp  Low  AVal      "

    if len(exp_energies) != 0 :
        header = header + 'UpperWN  UppJ   LowerWN  LowJ       WLnm'

    #f.write(header+"\n")
    for ii in range(0,num_levels):
        for jj in range(ii+1,num_levels):
            lower = ii+1#python is 0 indexed, but we are 1-indexed
            upper = jj+1
            avalue = transition_matrix[ii,jj]
            #this clunky formatting is to make it the same as the grasprad code
            #i could use fortranformat to make it exactly the same,
            #but this is a dependency not many will have, or be bothered to intsall.
            string_to_be_written = '{:4}'.format(upper)
            string_to_be_written += '{:4}'.format(lower)
            avalue_to_be_written = round(avalue,ndigits=2)
            string_to_be_written += '{:}'.format(" ")
            avalue_string = '{:.2E}'.format(avalue)
            avalue_string = avalue_string.replace("E","")
            string_to_be_written += avalue_string 

            if len(exp_energies) != 0:
                energy_string = '{:10}'.format(exp_energies[jj,0])
                energy_string += " "+'{:5}'.format(exp_energies[jj,1])

                energy_string += '{:10}'.format(exp_energies[ii,0])
                energy_string += " "+'{:5}'.format(exp_energies[ii,1])

                wavelength = 1e7 / np.abs(exp_energies[ii,0]-exp_energies[jj,0])
                energy_string += " "+'{:10}'.format(round(wavelength,4))
                 
                string_to_be_written += energy_string
                
            string_to_be_written += "\n"
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
    num_levels_1 = []
    found_array = []
    first_line_levels = 0
    def getNelec():
        for i in range(0,numlines):
            current_line = content[i]
            current_line_list = current_line.split()
            if len(current_line_list) > 0:
                if current_line_list[0] == 'Input':
                    break 
        flag = True 
        i+=4
        nelec =0 
        while flag:
            current_line = content[i].split()
            if current_line[0] != 'ANG':
                nelec += int(current_line[1])
            else:
                break 
            i+=1 
        i = 0
        return nelec
            
    nelec = getNelec()    
        
                
    for i in range(0,numlines):
        current_line = content[i]
        current_line_list = current_line.split()
        
        
        if current_line_list != []: 
            if current_line_list[0] == 'Z':
                atomicnumber = int(float(current_line_list[-1]))
            if current_line_list[-1] == "lowest":
                current_line = current_line.replace("\n","")
                print("found levels")
                print(current_line)
                first_line_levels = i
                break
    first_line_levels = first_line_levels + 5
    last_line_levels = first_line_levels 

    for i in range(first_line_levels,numlines):
        if content[i].split() != []:
            last_line_levels+=1 
        else:
            break         
    
    print(first_line_levels,last_line_levels)
    
    for i in range(last_line_levels,numlines):
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
    
    level_data    = np.loadtxt(file,skiprows=first_line_levels,max_rows=last_line_levels-first_line_levels,dtype=str)
    #print(level_data)
    electric_data = np.loadtxt(file,usecols=[0,1,2,3],skiprows=first_line_with_aij_index_electric,max_rows=last_line_electric-first_line_with_aij_index_electric)
    magnetic_data = np.loadtxt(file,usecols=[0,1,2,3],skiprows=first_line_with_aij_index_magnetic,max_rows=last_line_magnetic-first_line_with_aij_index_magnetic)
    num_levels_1 = 0
    if len(electric_data[:,1] != 0):
        num_levels_1 = max(electric_data[:,1])
    num_levels_2 = 0
    if len(magnetic_data[:,1] != 0):
        num_levels_2 = max(magnetic_data[:,1])
    #print(num_levels_1,num_levels_2)
    return electric_data,magnetic_data,found_array,int(max(num_levels_1,num_levels_2)),level_data,atomicnumber,nelec

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
parser.add_argument('-l', '--levels',  help='max level of output transitions, is the number of found levels by default (untested)')
#parser.add_argument('-g', '--globbing',  help='Selects datafiles by globbing')
parser.add_argument('-d', '--datafiles', nargs='+', help='Paths of input files. Multiple files put a space between them.')
parser.add_argument('-r', '--reorder', help='Path of reorder file - i.e translates unshifted state indices to their shifted ')
parser.add_argument('-e', '--energies', help='Path of energies (in exp. order for now, i might change this) ')

args = parser.parse_args()
#num_levels = 40
def writeadf04(level_data,avalues,max_n,levels_found,charge,nelec,args=np.array([]),shifted=np.array([])):
    from pathlib import Path
    ROOT_DIR = Path(__file__).parent
    TEXT_FILE = ROOT_DIR / 'ion_energy.dat'
    ip = np.loadtxt(TEXT_FILE)
    TEXT_FILE_terms = ROOT_DIR / 'ion_terms.dat'
    tt = np.loadtxt(TEXT_FILE_terms,dtype='>U2')
    ion_stage_index = charge - nelec 
    atomic_number_index = charge - 1 
    ion_pot = ip[ion_stage_index,atomic_number_index]
    ion_term = tt[ion_stage_index,atomic_number_index]
    
    g = open('adf04.vr.axel.dat','w')
    g.write('{:2}+{:2}        {:2}        {:2}    {:11.1f}.({:2})\n'.format(
        elements[charge-1],
        charge-nelec,
        charge,
        charge-nelec,
        ion_pot,
        ion_term
    ))
    
    nlevels = len(level_data[:,0])
    max_n = min(levels_found,max_n,nlevels)
    #print(levels_found,max_n,nlevels)
    parities    = np.ones(nlevels)
    statweights = np.zeros_like(parities)
    temps = np.array([1.00E+03, 
                      1.50E+03, 
                      1.80E+03, 
                      2.00E+03, 
                      2.50E+03, 
                      5.00E+03, 
                      7.50E+03, 
                      1.00E+04, 
                      1.50E+04, 
                      1.80E+04, 
                      2.00E+04, 
                      3.00E+04, 
                      4.00E+04, 
                      5.00E+04, 
                      6.00E+04, 
                      7.00E+04, 
                      8.00E+04, 
                      9.00E+04, 
                      1.00E+05])
    temps = 10 ** np.linspace(2,6,19)
    
    ups = np.ones_like(temps)
    
    if len(shifted)>0:
        for jj in range(0,len(shifted)):
            level_data[jj,-1] = float(shifted[jj])
    
    if len(args > 0):
        level_data = level_data[args,:]
    print(level_data)
        
        
    
    
    for jj in range(0,levels_found):
        
        jvalue = level_data[jj,1]
        #print(len(level_data[:,1]))
        if ((len(jvalue) >2)  and (jvalue[-2] == '/')):
            statweights[jj] = float (jvalue[0:-2]) + 1
        else:
            statweights[jj] = 2.0 * float (jvalue) + 1
        
        parity = level_data[jj,2]
        
        if parity == 'odd':
            parities[jj] = -1.0 
        elif parity != 'even':
            print('parity error in level ',jj)
        
    form = '{:3} Level{:3} (0)0({:4.1f}) {:13.4f}\n';
    
    for ii in range(0,levels_found):
        en = float(level_data[ii,-1])
        en *= RYDBERG_CM
        g.write(form.format(ii+1,
                            ii+1,
                            (statweights[ii]-1)*0.5,
                            en)) 
        
    firstPart = '{:5}{:5} {:7.2e}'
    middle = ' {:7.2e}'
    end = '{}{:7.2e}\n'
    
    
    
    g.write('   -1\n')
    temperature_header = ' 3.00    3'
    temperature_header += ' '*8
    for ii in range(0,len(temps)):
        temperature_header += middle.format(temps[ii])
    temperature_header+= '\n'
    temperature_header=temperature_header.replace('e','')
    #print(temperature_header)
    g.write(temperature_header)
    for ii in range(0,levels_found):

        for jj in range(0,ii):
            stringgg = ''
            allowed = True 
            
            weightdiff = abs(statweights[ii] - statweights[jj])
            paritydiff = parities[jj] * parities[ii]
            
            parallowed = paritydiff == -1.0 
            angallowed = (weightdiff == 0.0) or (weightdiff == 2.0)
            notzerozero   = not ( (statweights[ii] == 1.0) and (statweights[jj] == 1.0) )

            if (not notzerozero): 
                allowed = False 
            else:
                if (not angallowed):
                    allowed = False
                else:
                    paritydiff = parities[jj] * parities[ii]
                    if (paritydiff != -1.0):
                        allowed = False 
            inf = 0
            if allowed:
                wl = abs( float(level_data[ii,-1]) -  float(level_data[jj,-1]))  * RYDBERG_CONSTANT_CM 
                wl = 1e7/wl #nm 
                ups[:] = van_regemorter_single_ionised(
                    temps,
                    abs( float(level_data[ii,-1]) -  float(level_data[jj,-1])),
                    statweights[ii],
                    avalues[ii,jj]
                    )
                
                #inf = avalues[ii,jj]
                inf = 4 * wl**3 * statweights[ii] * avalues[ii,jj]
                inf /= 3*2.03e15
                #print('vr')
                marker='-'
            else:
                inf = 0 
                marker = ' '
                ups[:] = AXELROD_CONST_FOR_NOW * statweights[ii] * statweights[jj]
            
            for kk in range(0,len(temps)):
                stringgg += middle.format(ups[kk])
            
            g.write(
                 ( firstPart.format(ii+1,jj+1,avalues[ii,jj]) + stringgg + end.format(marker,inf)).replace('e','')
                )
    g.write('  -1\n')
    g.write('  -1  -1\n')
    g.write('C generated by Leo Patick Mulhollands GraspToVRAcode')
    g.close()
    
    return 0 

def main(output_file_name):
    found_array = []
    found_num_levels = []
    obtained_transition_arrays = []

    
    for file in input_array:
        electric,magnetic,found,num_levels,level_data,atomicnumber,nelec = readin_until_finds_aij(file)
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
        if args.levels:
            if int(args.levels) > num_levels:
                print("too many levels requested, printing default = ",num_levels)
            else:
                
                print("truncating original",num_levels,"levels to ",int(args.levels),"levels")
                num_levels = int(args.levels)
        exp_energies = []
        shift = np.array([]) 
        arg =   np.array([])
        if args.energies:
            shift = np.loadtxt(args.energies)
            arg = np.argsort(shift)
            #print(shift[arg])
            #print(arg)
            #ol = np.arange(0,len(arg),1,dtype=int)
            #print("REORDERING LEVELS ACCORDING TO INPUT FILE ")
            #transition_matrix = reorder(transition_matrix,args.reorder,arg,ol)
            
            #print(len(ol))
            
        else:
            exp_energies = []
        #print(level_data)
        writeadf04(level_data,transition_matrix,30,num_levels,atomicnumber,nelec,arg,shift)

        output(transition_matrix,output_file_name,num_levels=num_levels,exp_energies =exp_energies) 


if args.datafiles:
    input_array = manualdatafiles(args.datafiles) 
    output_file_name = args.name if args.name else "MERGED_A_VALUES.OUT"
    main(output_file_name)
else:
    print("no data files specified - stopping and printing help")
    parser.print_help()


