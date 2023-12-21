Takes many GRASP.OUTS , reads the A values and adds them. Detects if duplicate data is parsed.

usage: 

python3 avalue_combiner.py -d grasp.out1 grasp.out2 ...

needs spaces between the file names.

optional arguments: -n desired_file_name - outputs to a different file.
                    -r reorder_file - reads in a file in the dord format from the other codes.
                                      reorders a_values according to shift. All that is necessary is the index adjustments, the namelist in the old dord will be ignored if it is present. presently needs all desired levels to be re-indexed. (if they do not change you can just do something like i   i in the dord-esque file.)