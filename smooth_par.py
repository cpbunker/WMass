
from __future__ import division
import ROOT as R
import math
import array
import random
import time
import sys
import filecheck
import makeplots
import treesuite as t

'''
Module for the process of finding the smooth curves f=f(A...), and smoothing the raw data histograms to these curves via reweighting
'''

#this iteration is designed to run for only one binn and save the output to a .root


################################ Fill binns with hists ################################ 

start_time_1=time.clock()

#### get batch variables from the command line, define other important variables here #####

var=sys.argv[1] #tells which variable to make plot for (see options in variable list)

#screen incorrect var argument here
variables=['eta','pt','phi','cost','gen2d','reco2d','bos2d','dress2d']
vars1d=['lep_eta','lep_pt','phi','cost']

if (not (var in variables) ):
    # so we want a purposeful error to hold the job
    print 'variable or function entry unsupported'
    d1= 'hello'+5

bgroup_index=int(sys.argv[2]) #tells which pT, y binn to get data for

if sys.argv[3] not in vars1d:
    TwoD=True

else: 
    TwoD=False

n_trees=250 #how many trees to (try to) get data out of (variable)

n_resamples=100  #how many resamples there are per batch 

n_batches= 1 #how many batches of resamples there are

total_resamples= n_resamples*n_batches #how many resamples contribute to the end result in total

#choose edges based on number of bins
n_pt=6 #how many pT bins
n_y=3 #how many y bins

if n_pt==1 and n_y==1:
    ptedges,yedges=[0,1336],[0,8.85]

if n_pt==8 and n_y==4:
    ptedges,yedges=[0, 3.87, 6.75, 10.32, 15.08, 21.88, 32.15, 50.04, 1336], [0, 0.91, 1.84, 2.83, 8.85]

if n_pt==6 and n_y==3:
    ptedges,yedges=[0, 4.78, 9.03, 15.08, 24.82, 42.48, 1356], [0, 1.21, 2.49, 8.85]

if n_pt==5 and n_y==4:
    ptedges,yedges=[0, 5.54, 11.16, 20.29, 37.81, 1356], [0, 0.91, 1.84, 2.83, 8.85]

if n_pt==10 and n_y==8:
    ptedges,yedges=[0, 4.0, 6.0, 9.0, 12.0, 16.0, 21.0, 28.0, 38.0, 57.0, 1247.0], [0, 0.46, 0.91, 1.37, 1.84, 2.32, 2.83, 3.41, 7.88]

#update the edges specifically for this parallelization
new_pt, new_y = t.change_edges(ptedges,yedges,bgroup_index)

#create binn holder
bholder=t.holder(new_pt,new_y)

print '\n '+'*'*40+' \n filling binns with hists \n '+'*'*40+' \n'

#iter thru all the binns in the holder and retrieve hists

#we don't iterate through all bins anymore, but instead pick only the one bgroup
bgroup=bholder.bgroups()[0]

for b in bgroup:	

    #store the histogram data in the binn
    b.gethists(0) #non-resampled hists always come from the 0 batch

    #get all resamples
    b.getresamples(n_batches,n_resamples)

end_time_1=time.clock()

#hello= 9+ 'see'


################################ Smooth trees and make plots ################################ 

start_time_2=time.clock()

print ' \n '+'*'*40+' \n Smooth trees and make plots \n '+'*'*40+' \n'

#smooth the trees, filling reweightor
reweights=bholder.fill_smooth(n_batches,n_resamples)


#get hists of all four data types    
#returns four 1 val lists
d1,d2,d3,d4 = bholder.comp_all(reweights,var,hist=True)

#get what we want out of lists
raw, smooth, rms, boot = d1[0],d2[0],d3[0],d4[0]



#choose whether to make plots and how to make plots or write to a .root
root=False


if (var=='reco2d' or var=='dress2d' or var=='bos2d') and TwoD: #this makes the reco2d and dress2d graphs

    #condense gets rid of very low cross section points, remember these are Agraph objects 
    new_raw= makeplots.condense( makeplots.unroll(raw,var,'raw'), var)
    new_smooth= makeplots.condense( makeplots.unroll(smooth,var,'smooth'), var)
    new_rms= makeplots.condense( makeplots.unroll(rms,var,'rms'), var)
    new_boot= makeplots.condense( makeplots.unroll(boot,var,'boot'), var)

    makeplots.comp_plot(new_raw, new_rms, new_boot,var,bgroup_index,new_pt,new_y,total_resamples)

elif (not root) and TwoD: #this makes the gen2d graphs


    new_raw=  makeplots.unroll(raw,var,'raw')
    new_smooth= makeplots.unroll(smooth,var,'smooth')
    new_rms=  makeplots.unroll(rms,var,'rms')
    new_boot=  makeplots.unroll(boot,var,'boot')

    makeplots.old_comp_plot(new_raw, new_rms, new_boot,var,bgroup_index,new_pt,new_y,total_resamples)

elif not TwoD and not root:

    new_raw, new_rms, new_boot = makeplots.project3(raw,rms,boot,var,x=False)
    makeplots.old_comp_plot(new_raw, new_rms, new_boot,sys.argv[3],bgroup_index,new_pt,new_y,total_resamples)

    for i in range(1,20):
	print new_rms.GetBinContent(i)/10**6

elif root:

    new_raw=  makeplots.unroll(raw,var,'raw')
    new_smooth= makeplots.unroll(smooth,var,'smooth')
    new_rms=  makeplots.unroll(rms,var,'rms')
    new_boot=  makeplots.unroll(boot,var,'boot')

    #make a new .root file and write the results
    myfile=R.TFile(var+str(bgroup_index)+'.root','update')

    raw.Write()
    smooth.Write()
    rms.Write()
    boot.Write()

    new_raw.Write()
    new_smooth.Write()
    new_rms.Write()
    new_boot.Write()

    myfile.Close()


#save the graph to a .root

end_time_2=time.clock()



################################ Print Timing Profile ################################ 

print 'Timing Profile'

fill_time=(end_time_1 - start_time_1)/60

print 'Filling binns took '+str(fill_time)+' minutes'

graph_time=(end_time_2 - start_time_2)/60

print 'Making graphs took '+str(graph_time)+' minutes'











