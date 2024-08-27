import os

def create_source_file_list(src_dir):
    path = src_dir+'/'
    # obtain list of source files
    source_file_roots = []
    # from pflotran object file text file
    for line in open(path+'pflotran_object_files.txt','r'):
        w = line.split('}')
        if len(w) == 2:
            w2 = w[1].split('.o')
            source_file_roots.append(path+w2[0])

    # pflotran files specifying 'program'
    source_file_roots.append(path+'pflotran')
    source_file_roots.append(path+'pflotran_rxn')
    source_file_roots.append(path+'pflotran_derivative')

    # clm-pflotran files
    path = src_dir+'../clm-pflotran/'
    source_file_roots.append(path+'clm_pflotran_interface_data')
    source_file_roots.append(path+'mapping')
    source_file_roots.append(path+'pflotran_interface_main')
    source_file_roots.append(path+'pflotran_model')

    # unit tests
    path = src_dir+'/unittests/'
    for file in os.listdir(path):
        if file.endswith('.pf'):
            source_file_roots.append(path+file)

    # alphabetize
    source_file_roots.sort()
    #print(source_file_roots)
    #os.remove('pflotran_provenance.F90')
    file_count = 0
    source_file_list = []
    for root in source_file_roots:
        file_count += 1
        if not root.endswith('pf'):
            filename = root+'.F90'
        else:
            filename = root
        source_file_list.append(filename)
    #print(source_file_list)
    print('{} sources files in list'.format(file_count))
    return source_file_list

