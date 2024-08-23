import sys
import shutil
import os
import fnmatch
import re


myflags = re.I

routine_list = []
routine_list.append(['MatSetValuesLocal','PetUtilMatSVBL'])

def compile_regex(list_):
    new_list = []
    for entry in list_:
        new_list.append([re.compile(r'\s*call\s*'+str(entry[0])+r'\s*\('),entry[0],entry[1]])
    return new_list

#regex_matsetval = re.compile(r'(?P<pre>\s*call\s*)MatSetValuesLocal\s*\(\s*(?P<arg>[^(&]+)(?P<ending>[\)\&])\s*\w*',flags=myflags)
#regex_matsetval = re.compile(r'(?P<pre>\w*)\s*\)MatSetValuesBlockedLocal\s*\(\s*(?P<arg>\w*)',flags=myflags)

def refactor_routine(routine_info,line,fin,fout):
    old_name = routine_info[1]
    new_name = routine_info[2]
    regex_string = r'(?P<pre>\s*call\s*)'+\
                   str(old_name)+\
                   r'\s*\(\s*(?P<arg>[^(&]+)(?P<ending>[\)\&])\s*\w*'
    regex_ = re.compile(regex_string,flags=myflags)
    m = re.search(regex_,line)
    if not m:
        print('error')
        sys.exit(0)
    arguments = []
    for entry in m.group('arg').strip().split(','):
        if len(entry) > 0:
            arguments.append(entry)
    last_character = m.group('ending')
    while last_character == '&':
        another_line = fin.readline().strip()
        if another_line.startswith('!'):
            continue
        if len(another_line.split('&')) > 1:
               last_character = '&'
        else:
            last_character = ''
        args = re.split(r'!|\)\s*;',another_line)[0]
        args = re.split(r'\s*,\s*',args)
        for entry in another_line.split(')')[0].split(','):
            if len(entry) > 0:
                arguments.append(entry)

    new_pre = m.group('pre').rstrip()+' '+new_name+'('
    indentation = ' '*len(new_pre)
    line_out = new_pre
    num_args = len(arguments)
    for iarg in range(num_args):
        entry = arguments[iarg]
        if (iarg < num_args-1):
            new_segment = entry+','
        else:
            new_segment = entry+');CHKERR('+entry+')'
        if len(line_out+new_segment) > 80:
            fout.write(line_out+' &\n')
            line_out = indentation+new_segment
        else:
            line_out += new_segment
    fout.write(line_out+'\n')


def refactor(filename,regex_list):
    fin = open(filename,'r')
    fout = open(filename+'.tmp','w')
    prev_line = ''
    for line in fin:
        found = False
        for routine in regex_list:
            m = re.search(routine[0],line)
            if m:
                refactor_routine(routine,line,fin,fout)
                found = True
                break
            # arguments = []
            # for entry in m.group('arg').strip().split(','):
            #     if len(entry) > 0:
            #         arguments.append(entry)
            # last_character = m.group('ending')
            # while last_character == '&':
            #     another_line = f.readline().strip()
            #     if another_line.startswith('!'):
            #         continue
            #     if len(another_line.split('&')) > 1:
            #            last_character = '&'
            #     else:
            #         last_character = ''
            #     args = re.split(r'!|\)\s*;',another_line)[0]
            #     args = re.split(r'\s*,\s*',args)
            #     for entry in another_line.split(')')[0].split(','):
            #         if len(entry) > 0:
            #             arguments.append(entry)

            # new_pre = m.group('pre').rstrip()+' PetUtilMatSVBL('
            # indentation = ' '*len(new_pre)
            # line_out = new_pre
            # num_args = len(arguments)
            # for iarg in range(num_args):
            #     entry = arguments[iarg]
            #     if (iarg < num_args-1):
            #         new_segment = entry+','
            #     else:
            #         new_segment = entry+');CHKERR('+entry+')'
            #     if len(line_out+new_segment) > 80:
            #         fout.write(line_out+' &\n')
            #         line_out = indentation+new_segment
            #     else:
            #         line_out += new_segment
            # fout.write(line_out+'\n')
        if not found:
            fout.write(line)
    fin.close()
    fout.close()



def main():
    print('\nRefactoring MatSetValues calls\n')
    regex_list = compile_regex(routine_list)
    suffix = '*.F90'
    for root, dirnames, filenames in os.walk('.'):
        for filename in fnmatch.filter(filenames,suffix):
            filename = os.path.join(root,filename)
#            print(filename)
            refactor(filename,regex_list)


if __name__ == "__main__":
    try:
        suite_status = main()
        print("success")
        sys.exit(suite_status)
    except Exception as error:
        print(str(error))
#        if cmdl_options.backtrace:
#            traceback.print_exc()
#        traceback.print_exc()
        print("failure")
        sys.exit(1)

print('done')