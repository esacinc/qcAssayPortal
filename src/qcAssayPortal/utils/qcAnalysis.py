import os
import sys
import time
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET

def select_rows_old(df,search_strings):
	# This function is to return the rows which contain all of the strings in search_strings
	unq,IDs = np.unique(df,return_inverse=True)
	unqIDs = np.searchsorted(unq,search_strings)
	df_out1 = df[((IDs.reshape(df.shape) == unqIDs[:,None,None]).any(-1)).all(0)]
	df_out2 = df[~((IDs.reshape(df.shape) == unqIDs[:,None,None]).any(-1)).all(0)]
	return df_out1, df_out2

def string_finder(row, words, waived_list):
	# row is a type of Series.
	row1 = row.tolist()
	row2 = [v for u, v in enumerate(row1) if u not in waived_list]
	if any(word == field for field in row2 for word in words):
		return True
	return False

def select_rows(df, match, waived_col_id_list):
	ids = df.apply(string_finder, words=match, waived_list=waived_col_id_list, axis=1)
	df_out1 = df[ids]
	df_out2 = df[~ids]
	return df_out1, df_out2

def locate_missing(valueList, colNameList, match, waived_list):
	outList = [colNameList[index] for index, value in enumerate(valueList) if (value in match) and (index not in waived_list)]
	return "Error", "Attribute", "Essential attribute(s) has(have) missing values, including "+'; '.join(outList)+'.'

def qcAnalysisGlobal(experiment_type, skyTsvDirList, fileNameList, waived_col_dic):
	# Check all the *.tsv files in skyTsvDirList to make sure all of them exist,
	# Otherwise, it means some *.sky.zip files are not successfully transformed into *.tsv file via SkylineCMd and an error will be yielded.
	failedList = [fileNameList[i] for i, item in enumerate(skyTsvDirList) if not os.path.isfile(item)]
	if len(failedList) >0:
		print >> sys.stderr, "The file(s) which is(are): %s can't be opened by Skyline. Please check the version compatibility between Skyline document and Skyline program."%(", ".join(failedList))
		sys.exit(1)
	# touch the first skyFile
	dfTemplate = pd.read_csv(skyTsvDirList[0], sep='\t', header=0, converters={i: str for i in range(0, 100)})
	# Because the exported .tsv files are based on the skyline report template, they must have the same column information.
	# If one sky document file lacks one column data, the exported .tsv file will assign missing values to this this column.
	errorDf = pd.DataFrame(columns=['SkyDocumentName', 'IssueType', 'IssueSubtype', 'IssueReason']+list(dfTemplate.columns.values))
	normalDf = pd.DataFrame(columns=['SkyDocumentName']+list(dfTemplate.columns.values))
	# Step 1: detect missing values "" without regard to the waived_col of the specific the experiment type
	search_strings = ['']
	waived_col_list = waived_col_dic[experiment_type]
	for i, skyFileDir in enumerate(skyTsvDirList):
		# Read the *.tsv file into dataframe and keep all the values in the format of string. 
		df = pd.read_csv(skyFileDir, sep='\t', header=0, converters={i: str for i in range(0, 100)})
		# Add 'SkyDocumentName' into df as the first column
		col_name1 = df.columns.tolist()
		col_name1.insert(0, 'SkyDocumentName')
		df = df.reindex(columns=col_name1)
		df['SkyDocumentName'] = fileNameList[i]
		waived_col_id_list = [s for s, item in enumerate(df.columns.tolist()) if item in waived_col_list]
		#print waived_col_id_list
		dfTmp1, dfTmp2 = select_rows(df, search_strings, waived_col_id_list)
		#print dfTmp1.shape[0]
		removedPeptideList = []
		if dfTmp1.shape[0] > 0:
			col_name = dfTmp1.columns.tolist()
			col_name.insert(1,'IssueType')
			col_name.insert(2,'IssueSubtype')
			col_name.insert(3,'IssueReason')
			dfTmp1 = dfTmp1.reindex(columns=col_name)
			waived_col_id_list_update = [s for s, item in enumerate(dfTmp1.columns.tolist()) if item in waived_col_list]
			dfTmp1['IssueType'], dfTmp1['IssueSubtype'], dfTmp1['IssueReason']= dfTmp1.apply(lambda row: locate_missing([row[colName] for colName in dfTmp1.columns.values], list(dfTmp1.columns.values), search_strings, waived_col_id_list_update), axis=1)
			errorDf = pd.concat([errorDf, dfTmp1],  ignore_index=True)
			# Deduplicate the 'PeptideModifiedSequence' with the PrecursorCharge
			peptideList = list(set(dfTmp1['PeptideModifiedSequence']+'$$$$'+dfTmp1['PrecursorCharge']))
			#peptideList = list(set(dfTmp1['PeptideModifiedSequence']))
			removedPeptideList = removedPeptideList + peptideList
		# The peptides in dfTmp1[''] will be removed from dfTmp2, because the peptide information is 
		# incomplete and these peptides should not be exported into *.tsv file for downstream R code.
		# Pay attention: If part of data for one peptide sequence have missing value and the rest part don't have missing value, all the rows of this peptide with specific precursor charge will be removed.
		for removedPeptide in removedPeptideList:
			removedPeptide1 = removedPeptide.split('$$$$')[0]
			removedPeptide2 = removedPeptide.split('$$$$')[1]
			dfTmp2 = dfTmp2[~((dfTmp2['PeptideModifiedSequence'] == removedPeptide1) & (dfTmp2['PrecursorCharge'] == removedPeptide2))]
		
		#dfTmp2 = dfTmp2[~dfTmp2['PeptideModifiedSequence'].isin(removedPeptideList)]
		# Detect whether there are PeptideModifiedSequences from dfTmp2 which exist in normalDf
		# Because current normalDf is empty, it's unnecessary to do it now. Comment the bellow codes.
		#normalDfPeptideList = list(set(normalDf['PeptideModifiedSequence']))
		#dfTmp2Peptideist = list(set(dfTmp2['PeptideModifiedSequence']))
		#duplicatedPeptideList = [item for item in dfTmp2Peptideist if item in normalDfPeptideList]
		#dfTmp3 = dfTmp2[dfTmp2['PeptideModifiedSequence'].isin(duplicatedPeptideList)]
		#col_name1 = dfTmp3.columns.tolist()
		#col_name1.insert(1,'IssueType')
		#col_name1.insert(2,'IssueSubtype')
		#col_name1.insert(3,'IssueReason')
		#dfTmp3 = dfTmp3.reindex(columns=col_name1)
		#dfTmp3['IssueType'] = 'Error'
		#dfTmp3['IssueSubtype'] = 'Duplicated peptide'
		#dfTmp3['IssueReason'] = 'Duplicated peptide modified sequences exist'
		# add the duplicated peptideModifidedSequence into errorDf
		#errorDf = pd.concat([errorDf, dfTmp3],  ignore_index=True)
		# remove the duplicated peptideModifidedSequence in dfTmp2
		#dfTmp2 = dfTmp2[~dfTmp2['PeptideModifiedSequence'].isin(duplicatedPeptideList)]
		
		# Add the dfTmp2 into normalDf
		normalDf = pd.concat([normalDf, dfTmp2], ignore_index=True)
	return errorDf, normalDf

def detectIS(skyFileDir, fileName, experiment_type, error_report_path, errorDfColNumber):
	# Parse the *.sky file whose format is XML.
	for event, elem in ET.iterparse(skyFileDir):
			if event == 'end':
				if elem.tag == 'peptide_modifications':
					internal_standard= elem.get('internal_standard', default=None)
					if (internal_standard is None or internal_standard == 'none') and (experiment_type=='exp1' or experiment_type=='exp2'):
					# This means that the user doesn't set the Internal Standard Type when preparing data for upload
					 	#print >> sys.stderr, "Internal_standard value in the peptide_modifications underneath peptide_settings of the *.sky file of %s is unset. Please check it."%(skyFileDir)
						if internal_standard is None:
							errorInfor = '\t'.join([os.path.basename(fileName), 'Error', 'Internal standard','Internal standard type in the peptide_modifications underneath peptide_settings is not set.']+['']*(errorDfColNumber-4))+'\n'
						else:
							errorInfor = '\t'.join([os.path.basename(fileName), 'Error', 'Internal standard','Internal standard type in the peptide_modifications underneath peptide_settings is set to be none.']+['']*(errorDfColNumber-4))+'\n'
						with open(error_report_path, 'a') as outfTmp:
							outfTmp.write(errorInfor)
						internal_standard_type = 'unset'
					elif internal_standard is None or internal_standard == 'none':
						internal_standard_type = 'none'
					else:
						internal_standard_type = internal_standard
			elem.clear()
	return internal_standard_type

def qcAnalysisRcode(experiment_type, error_report_path, dataset_path, fileList_path, mypeptideType_file_path, RscriptBinary, rScript, plot_output, plot_output_dir):
	if experiment_type == 'exp2':
		os.system('"%s" %s %s %s %s %s >> %s'%(RscriptBinary, rScript, dataset_path, fileList_path, plot_output, plot_output_dir, error_report_path))
	if experiment_type == 'exp1':
		os.system('"%s" %s %s %s %s %s %s >> %s'%(RscriptBinary, rScript, dataset_path, fileList_path, plot_output, plot_output_dir, mypeptideType_file_path, error_report_path))