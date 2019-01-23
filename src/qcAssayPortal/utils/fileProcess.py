import os

def headerModify(inf_tmp, outf_tmp):
    # This function is to remove space in the column name at the first row of inf_tmp.
    outf = open(outf_tmp, 'w')
    with open(inf_tmp, 'r') as inf:
        line = inf.readline()
        outf.write(line.replace(' ', ''))
        line = inf.readline()
        while line != '':
            # replace #N/A with '', so that when R read this file, it will transform '' into NA automatically.
            outf.write(line.replace('#N/A', ''))
            line = inf.readline()
    outf.close()
    os.remove(inf_tmp)