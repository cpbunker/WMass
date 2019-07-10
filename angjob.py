
from __future__ import division
import ROOT as R
import math
import array
import time
import sys
import numpy
import random
import treesuite as t


'''
Job script for submission to ht condor

Takes a list of pt and y edges and creates appropriate bins

Then in each binn we will:
	Fill with main and friend tree data
	Put data into hists
	write hists to root file


'''


################################ Get Data from TTrees ################################ 

print 'Get Data from TTrees'

#important variables for job submission I/O
#these come from the command line/wrap job arguments

tree_number=int(sys.argv[1]) #tree number the data is being pulled from

batch_number=int(sys.argv[2]) #resampling batch this job is a part of 

trees_per_file=50 #number of trees of data in each WData_part.root file

part_number= int(math.floor(tree_number/trees_per_file) ) #tells part of WData file

n_resamples=100 #how many resamples to perform for this job

resample_start= batch_number*n_resamples #where to start counting resamples from


#Get TTree
file1=R.TFile('/afs/cern.ch/work/c/cbunker/trees/WData_part'+str(part_number)+'.root')
tree1=file1.Get('WData_tree'+str(tree_number)+'.root')

print 'trees/WData_part'+str(part_number)+'.root','WData_tree'+str(tree_number)+'.root'

################################ Put Data into Binns ################################ 

start_time_1=time.clock()

print 'Put Data into Binns'

#choose edges based on number of bins
n_pt=6
n_y=3

if n_pt==1 and n_y==1:
    ptedges,yedges=[0,1336],[0,8.85]


if n_pt==2 and n_y==1:
    ptedges,yedges=[0, 3.87, 6.75], [0, 0.91]

if n_pt==8 and n_y==4:
    ptedges,yedges=[0, 3.87, 6.75, 10.32, 15.08, 21.88, 32.15, 50.04, 1336], [0, 0.91, 1.84, 2.83, 8.85]

if n_pt==6 and n_y==3:
    ptedges,yedges=[0, 4.78, 9.03, 15.08, 24.82, 42.48, 1356], [0, 1.21, 2.49, 8.85]

if n_pt==5 and n_y==4:
    ptedges,yedges=[0, 5.54, 11.16, 20.29, 37.81, 1356], [0, 0.91, 1.84, 2.83, 8.85]

if n_pt==10 and n_y==8:
    ptedges,yedges=[0, 4.0, 6.0, 9.0, 12.0, 16.0, 21.0, 28.0, 38.0, 57.0, 1247.0], [0, 0.46, 0.91, 1.37, 1.84, 2.32, 2.83, 3.41, 7.88]


#update the edges specifically for this parallelization
#new_pt, new_y = t.change_edges(ptedges,yedges,job_number)
#print new_pt, new_y

#create binn holder
bholder=t.holder(ptedges, yedges,part=tree_number)


#get seed for resampling
seed=abs(hash(sys.argv[3]))
print seed
numpy.random.seed(seed)

#fill binn holder from tree
print 'start ttree loop'
bholder.Fill(tree1)
print 'tree filled'

#iter thru all bgroups to resample
for bgroup in bholder.bgroups():

    print bgroup[0]

    #make resampled bins
    for j in range(1+resample_start,1+n_resamples+resample_start):  
    #it is critical that this count start at 1, not 0

	#create the new, resampled binn
	bgroup[0].resample(j)

    #fill rest of binns in bgroup from first one's resamples
    t.clone_group(bgroup)

end_time_1=time.clock()


################################ Make Hists and Write to .root ################################ 

start_time_2=time.clock()

print 'Make Hists, Write to .root and make ttree'

filename='WJet_batch'+str(batch_number)+'_tree'+str(tree_number)+'.root'

#try making the file out here
myfile=R.TFile(filename,'update')
myfile.Close()

#iter thru all binns
for b in bholder.forin():

    print b

    #make hists
    b.process()

    #write to .root
    b.write(filename)

#make and write TTree
bholder.write_ttree(filename)

end_time_2=time.clock()


################################ Print Timing Profile ################################ 

print 'Timing Profile'

fill_time=(end_time_1 - start_time_1)/60

print 'Filling binns took '+str(fill_time)+' minutes'

write_time=(end_time_2 - start_time_2)/60

print 'Writing hists took '+str(write_time)+' minutes'








	




    



	    
	    









        
       
    


