
from __future__ import division
import ROOT as R
import math
import array
import time
import sys
import random
import os
import treesuite as t

'''
Module to read and report on the ht condor job output: log files in the /log directory and error files in the 
\error directory
'''

##################### Functions for checking files #####################

def check_shadow(lastlines):
    '''
    Checks if a shadow exception occured at the end of the log file
    Args, file0, the file to be checked
    Returns, message string, index of line tuple if it found it, otherwise None
    '''

    #get the file lines as a string
    #lines=file0.read().splitlines()

    #lastlines=lines[-40:]

    #check for error phrase
    key='Shadow exception!'
    index=-1

    #iter thru lines
    #want to know which line it was found at
    for i in range(len(lastlines)):

        if key in lastlines[i]:

	    #update so only the last instance is returned
	    index=i

    #check what to return
    if index != -1:
	return 'Job had a '+key+' on line '+str(index)+'\n'

    else:
	return ''




def check_disconnect(lastlines):
    '''
    Checks if a Job disconnected too long error occured at the end of the log file
    Args, file0, the file to be checked
    Returns, message string, index of line tuple if it found it, otherwise None
    '''

    #get the file lines as a string
    #lines=file0.read().splitlines()

    #lastlines=lines[-40:]

    #check for error phrase
    key='Job disconnected too long'
    index=-1

    #iter thru lines
    #want to know which line it was found at
    for i in range(len(lastlines)):

        if key in lastlines[i]:

	    #update so only the last instance is returned
	    index=i

    #check what to return
    if index != -1:
	return 'Job had a '+key+' on line '+str(index)+'\n'

    else:
	return ''



def check_IO(file0):
    '''
    Checks if an IO error occured at the end of the error file
    Args, file0, the file to be checked
    Returns, message string if it found it, otherwise null string
    '''

    #get the file lines as a string
    lines=file0.read().splitlines()

    #screen out empty output files
    if len(lines) == 0:
	return ''

    #error should be on the last line
    lastline=lines[-1]

    #check for error phrase
    key='TTree I/O error'

    if key in lastline:

	return 'Job had a '+key+'\n'

    else:
	return ''


def check_remove(lastlines):
    '''
    Checks if a periodic system remove occured at the end of the log file
    Args, file0, the file to be checked
    Returns, message string, index of line tuple if it found it, otherwise None
    '''

    #get the file lines as a string
    #lines=file0.read().splitlines()

    #lastlines=lines[-40:]

    #check for error phrase
    key='SYSTEM_PERIODIC_REMOVE'
    index=-1

    #iter thru lines
    #want to know which line it was found at
    for i in range(len(lastlines)):

        if key in lastlines[i]:

	    #update so only the last instance is returned
	    index=i

    #check what to return
    if index != -1:
	return 'Job had a '+key+' on line '+str(index)+'\n'

    else:
	return ''


def return_value(lastlines):
    '''
    Finds the return value at the end of the log file
    Args, file0, the file to be checked
    Returns, return value, if found, otherwise 'None found'
    '''

    #get the file lines as a string
    #lines=file0.read().splitlines()

    #lastlines=lines[-40:]

    #check for error phrase
    key='Normal termination (return value'
    index=-1

    #iter thru lines
    #want to know which line it was found at
    for i in range(len(lastlines)):

        if key in lastlines[i]:

	    #update so only the last instance is returned
	    index=i

    #get the line with the return value, if found
    if index != -1:
        myline=lastlines[index]
    else:
	return 'none found'

    #find return val in myline
    mystring = myline.split(key)[-1]

    #iter thru to find numbers

    #numbers as chars for finding numbers
    numlist=['0','1','2','3','4','5','6','7','8','9']

    #search for the return value
    return_val=-1
    mychar=''
    for char in mystring:
	if char in numlist:
	    mychar += char

    #update return value
    if mychar != '':
	return_val=int(mychar)


    return return_val


def time_profile(file0):
    '''
    Finds time profile as given at the end of the output file
    Args, file0, the file to be checked
    Returns, time profile, if found, or null string
    '''

    #get the file lines as a string
    lines=file0.read().splitlines()

    #last 3 lines contain profile
    lastlines=lines[-3:]

    #check for error phrase
    key='Profile'

    #iter thru lines
    found=False
    for line in lastlines:

        if key in line:
	    found=True

    #return none if not found
    if not found:
	return ''

    else:
	return lastlines[-2]+'\n'+lastlines[-1]+'\n'


def file_size(rootfile):

    size = int(os.path.getsize(rootfile))

    return size


def get_good_files(n_trees,batch):
    '''
    Returns a list of the .root files that are ok to use as well as a list of their job numbers
    '''

    print 'finding good files'
    good_files=[]
    good_trees=[]

    #iter thru jobs
    for tree_n in range(n_trees):

        #find the log, output and error files
	filestring='/eos/user/c/cbunker/'
        rootfilepath=filestring+'WJet_batch'+str(batch)+'_tree'+str(tree_n)+'.root'

        #check size of root file
        size=file_size(rootfilepath)

	#print rootfilepath, ' size is ', size

	if size > 2*10**6:
	    good_files.append(rootfilepath)
	    good_trees.append(tree_n)


    print str( len(good_files) )+ ' good files found'
		
    #return good_files

    return good_files, good_trees






##################### Check the files #####################


if __name__=='__main__':

    #get the range of files to check from the terminal
    #remember sys args are always taken as strings
    tree_start=int(sys.argv[1])
    tree_end=int(sys.argv[2])

    #also the part number
    batch=int(sys.argv[3])

    #make counts of the categories
    bad_list=[]
    success_count=0
    return0_count=0
    IO_shadow_count=0
    IO_count=0
    shadow_count=0
    disconnect_remove_count=0
    remove_count=0
    disconnect_count=0


    #iter thru jobs
    for tree in range(tree_start,tree_end+1):

        #find the log, output and error files
        #log=open('log/tree'+str(tree)+'.txt','r')
        #error=open('error/tree'+str(tree)+'.txt','r')
        #output=open('output/tree'+str(tree)+'.txt','r')

        log=open('log/batch'+str(batch)+'_tree'+str(tree)+'.txt','r')
        error=open('error/batch'+str(batch)+'_tree'+str(tree)+'.txt','r')
        output=open('output/batch'+str(batch)+'_tree'+str(tree)+'.txt','r')

        #rootfilepath='/eos/user/c/cbunker/WJet_part'+str(part)+'_job'+str(job)+'.root'
	rootfilepath= '/afs/cern.ch/work/c/cbunker/WJet_batch'+str(batch)+'_tree'+str(tree)+'.root'

        #get last lines of log file
        log_lines=log.read().splitlines()[-30:]

        #start with a string representation of the job
        report= 'WJet_batch'+str(batch)+'_tree'+str(tree)+' report \n'

        #get the return value
        return_val=return_value(log_lines)

        #convert to string
        return_string='Return value: '+str(return_val)+'\n'

        #get time profile
        time_string = time_profile(output)

        #check size of root file
        size=file_size(rootfilepath)
	size_string='The .root file is '+str(size)+' bytes or '+str(size/10**6)+'megabytes \n'

        #update report with this info
	report += return_string
	report += size_string
	report += time_string

        #update with error checks
        report += check_IO(error)
	report += check_shadow(log_lines)
	report += check_remove(log_lines)
	report += check_disconnect(log_lines)

        if not ( return_val == 0 and size > 2*10**6): #cutoff 2 MB
	    bad_list.append(tree)

        #update counts
        if return_val == 0 and size > 7*10**6: #cutoff 7 MB
	    success_count += 1
	elif return_val ==0:
	    return0_count += 1
        elif 'Shadow exception!' in report and 'I/O error' in report:
	    IO_shadow_count += 1
        elif 'Shadow exception!' in report:
	    shadow_count += 1
        elif 'I/O error' in report:
  	    IO_count += 1
	elif 'Job disconnected too long' in report and 'SYSTEM_PERIODIC_REMOVE' in report:
	    disconnect_remove_count += 1
        elif 'SYSTEM_PERIODIC_REMOVE' in report:
	    remove_count += 1
	elif 'Job disconnected too long' in report:
	    disconnect_count += 1



        #print results, if necessary
        print report


    #print overall summary
    print 'Summary report'
    print str(success_count)+' jobs ran successfully'
    print str(return0_count)+' jobs had a return value of 0 but an incomplete .root file'
    print str(IO_shadow_count)+' jobs faulted with a TTree I/O error and a shadow exception'
    print str(shadow_count)+' jobs had a shadow exception only'
    print str(IO_count)+' jobs had a TTree I/O error only'
    print str(disconnect_remove_count)+' jobs were removed with a job disconnected error by SYSTEM_PERIODIC_REMOVE'
    print str(disconnect_count)+' jobs had a job disconnected error only'
    print str(remove_count)+' jobs were removed by SYSTEM_PERIODIC_REMOVE only'

    #check all jobs have been accounted for
    if success_count+return0_count+IO_shadow_count+IO_count+shadow_count++disconnect_remove_count+disconnect_count+remove_count < (tree_end + 1) - tree_start:
        print 'Not all jobs accounted for in this report'

    if success_count+return0_count+IO_shadow_count+IO_count+shadow_count++disconnect_remove_count+disconnect_count+remove_count > (tree_end + 1) - tree_start:
        print 'Too many jobs accounted for in this report'

    print 'The bad files are: '
    print bad_list

    goodlist=''
    for i in range(250):
	if i not in bad_list:
	    goodlist += 'WJet_batch'+str(batch)+'_tree'+str(i)+'.root, '

    print 'The good files are: '
    print goodlist













