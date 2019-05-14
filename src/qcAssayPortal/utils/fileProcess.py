import os
import pandas as pd
import random

def headerModify(inf_tmp, outf_tmp, experiment_type_tmp):
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
    # Judge is the skyline template old or new if experiment_type_tmp is exp1 or exp2
    if experiment_type_tmp in ["exp1", "exp2"]:
        # Read outf_tmp into a dataframe and keep all the values in the format of string.
        df_tmp = pd.read_csv(outf_tmp, sep='\t', header=0, converters={i: str for i in range(0, 100)})
        if df_tmp.shape[0] < 10:
            sampledList1 = df_tmp['AnalyteConcentration'].values.tolist()
            sampledList2 = df_tmp['InternalStandardConcentration'].values.tolist()
        else:
            sampledList1 = random.sample(df_tmp['AnalyteConcentration'].values.tolist(), 10)
            sampledList2 = random.sample(df_tmp['InternalStandardConcentration'].values.tolist(), 10)
        if any([item != '' for item in sampledList1]) and any([item != '' for item in sampledList2]):
            skylineTemp_type = 'new'
        else:
            skylineTemp_type = 'old'
    else:
        skylineTemp_type = None
    return skylineTemp_type