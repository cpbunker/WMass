

from __future__ import division
import ROOT as R
import math
import array
import random
import numpy
import sys
import filecheck
import makeplots

"""
A module for creating objects that allow for more intuitive manipulation of TTree objects
"""


class ev(object):
    '''
    This is an event object, designed for use with W decay data in which the relevant parameters are:
    mass, pt, y, costcs, phics, weightGen
    
    Has those attributes as well as type= prefsrw or genw
    '''
    
    def __init__(self,event,dtype='prefsrw'):
        '''
        inits the ev with the relevant params from the TTree event
	dtype tells whether genw or prefsrw
        '''

        if dtype == 'genw':
	    self.type='genw'
        
            #def numerical attributes
            self.pt=event.genw_pt #transverse momentum in GeV
            self.y=event.genw_y #rapidity, unitless
            self.cost=event.genw_costcs #cos(theta), unitless, in collins-sopper frame
            self.phi=event.genw_phics #phi, radians, in collins-sopper frame
            self.w=event.weightGen #weight of event, unitless

	    #lep eta and pt need to be converted to tuples for immutability
	    self.lep_pt=tuple(event.LepGood_pt)
	    self.lep_eta=tuple(event.LepGood_eta)
	    self.dress_pt=tuple(event.GenLepDressed_pt)
	    self.dress_eta=tuple(event.GenLepDressed_eta)
            
        elif dtype =='prefsrw':
            self.type='prefsrw'
        
            #def numerical attributes
            self.pt=event.prefsrw_pt #transverse momentum in GeV
            self.y=event.prefsrw_y #rapidity, unitless
            self.cost=event.prefsrw_costcs #cos(theta), unitless, in collins-sopper frame
            self.phi=event.prefsrw_phics #phi, radians, in collins-sopper frame
            self.w=event.weightGen #weight of event, unitless

	    #lep eta and pt need to be converted to tuples for immutability
	    self.lep_pt=tuple(event.LepGood_pt)
	    self.lep_eta=tuple(event.LepGood_eta)
	    self.dress_pt=tuple(event.GenLepDressed_pt)
	    self.dress_eta=tuple(event.GenLepDressed_eta)
	    

    def __str__(self):
	return 'pt = '+str(self.pt)+', y = '+str(self.y)+', phi = '+str(self.phi)+', cost = '+str(self.cost)

    

    
class binn(object):
    '''
    Bin object contains a list of all events that fit into some slice of the 5 parameter event space
    The values that define the slice are called edges
    Attributes:
        evs, list of all events
        edges, a list of all edges
        massmin, massmax, ... phics min, phics mass, the individual edges
        datatype: string, whether to include genw data, prefsrw or both
	i, int, the index of the coefficient that this binn's data aims to calculate
	resamp: int, tells if the binn is a resampled one or not: 0 no, else yes
	number tells numbe rof resampling, this is important so that all hists have different names
	phics, costcs, both_hists (for A0 only) hists of raw ttree data
	minibinns: list of minibinns dividing pt, y space of binn into phi, cost space
	resamples, a list of binns which are the resampled counterparts of the original binn

    In general, the life of a bin looks like:
    Hist making bin: Fill, [resample, making a bunch of clones] process, write to file
    A_i calculating bin: gethists, plot(which utilises calc)

    '''

    #################### Here are the overloaded methods ####################

    def __init__(self,ptmin,ptmax,ymin,ymax,datatype,i,resamp=0,part=9999):

	#basic attributes
        self.evs=[]
        self.type=datatype
	self.i=i
	self.resamp=resamp
	self.part=part

        #only allow for pt, y edges for now!!!

        #pt edges
        self.ptmin=ptmin
        self.ptmax=ptmax

        #y edges: convert to |y|
        self.ymax=max(abs(ymin),abs(ymax))
        self.ymin=min(abs(ymin),abs(ymax))

	#init TH1Fs for data storage
	self.num, self.den, self.num_err, self.den_err = makehists(str(self))

	#phi and cost hists, which hold raw data for later smoothing process
	#these are only necessary for A0 binns
	if self.i ==0 or self.i == 8: #A8 does not exist, used to distinguish the reweightor binns
					#and avoid memory leaks

	    #create boson angle hists
	    self.nbins_phi=20
	    self.nbins_cost=20

	    self.phics=R.TH1F(str(self)+'phics','\phi;\phi;cross section',self.nbins_phi,0,2*math.pi)
	    self.costcs=R.TH1F(str(self)+'costcs','cos(#theta);cos(#theta);cross section',self.nbins_cost,-1,1)
	    self.gen2d=R.TH2F(str(self)+'gen2d','Monte Carlo boson angular data; \phi ; cos(#theta) ; cross section',self.nbins_phi,0,2*math.pi,self.nbins_cost,-1,1) 

	    #create lepton kinematic hists 
	    self.nbins_pt=20
	    self.nbins_y=20

	    self.lep_pt=R.TH1F(str(self)+'lep_pt','Lepton p_{T}',self.nbins_pt,0,120)
	    self.lep_eta=R.TH1F(str(self)+'lep_eta','Lepton \eta',self.nbins_y,-2.4,2.4)
	    self.reco2d=R.TH2F(str(self)+'reco2d','Monte Carlo lepton kinematic data; \eta; p_{T}; cross section',self.nbins_y,-2.4,2.4,self.nbins_pt,0,120)

	    #create dressed lepton kinematic hists
	    self.dress2d=R.TH2F(str(self)+'dress2d','Dressed lepton kinematic data; \eta; p_{T}; cross section',self.nbins_y,-2.4,2.4,self.nbins_pt,0,120)

	    #create boson kinematic hists
	    self.pt=R.TH1F(str(self)+'pt','Boson p_{T}',self.nbins_pt,0,120)
	    self.y=R.TH1F(str(self)+'y','Boson y',self.nbins_y,0,3)
	    self.bos2d=R.TH2F(str(self)+'bos2d','Monte Carlo boson kinematic data; y; p_{T}; cross section',self.nbins_y,-4,4,self.nbins_pt,0,120)

	    #call sumw2 for all these hists
	    self.phics.Sumw2()
	    self.costcs.Sumw2()
	    self.gen2d.Sumw2()

	    self.lep_pt.Sumw2()
	    self.lep_eta.Sumw2()
	    self.reco2d.Sumw2()

	    self.dress2d.Sumw2()

	    self.pt.Sumw2()
	    self.y.Sumw2()
	    self.bos2d.Sumw2()
	

	#init the resamples attribute as a list for storing the resampled version of this binn
	self.resamples=[]   

	#init TH2 to hold reweights if i==8
	if self.i==8:
	    self.reweights=R.TH2F( str(self)+'reweights','reweights',self.nbins_phi,0,2*math.pi,self.nbins_cost,-1,1)

        
    def __len__(self):
        #overload len operator
        return len(self.evs)

    def __str__(self):
	#overload string representation

	#part is included now
	#if self.part != 9999:
	    #returnstring=str(self.part)
	#else:
	returnstring = ''
	#split off representation of resampled binns
	if self.resamp ==0:
	    returnstring += 'A'+str(self.i)+'_'+str(self.ptmin)+'_'+str(self.ptmax)+'_'+str(self.ymin)+'_'+str(self.ymax)
	else:
	    returnstring += 'A'+str(self.i)+'_'+str(self.ptmin)+'_'+str(self.ptmax)+'_'+str(self.ymin)+'_'+str(self.ymax)+'_r'+str(self.resamp)

	return returnstring

    def __eq__(self,other):
	#overload equal to operator
	if self.i == other.i:
	    if self.ptmin== other.ptmin:
		if self.ptmax== other.ptmax:
		    if self.ymin == other.ymin:
			if self.ymax==other.ymax:
			    return True
	return False



    #################### Here are the basic methods ####################

    def nice_str(self):

	return str(self.ptmin)+' < p_{T} < '+str(self.ptmax)+' , '+str(self.ymin)+' < |y| < '+str(self.ymax)

    def contains(self,event):
        #determine if event is within the bin ranges and of the appropriate type
 	if event.prefsrw_pt >= self.ptmin:
	    if event.prefsrw_pt < self.ptmax:
		if event.prefsrw_y >= self.ymin:
		    if event.prefsrw_y < self.ymax:
			return True
        return False

    def contains_trimmed(self,event):
        #determine if event is within the bin ranges and of the appropriate type
	#for events from trimmed trees rather than original presfrw/genw events
 	if event.pt >= self.ptmin:
	    if event.pt < self.ptmax:
		if event.y >= self.ymin:
		    if event.y < self.ymax:
			return True
	return False
    

    def clone(self):
        '''
        returns a new but empty binn obj with the same params as the original
        '''
        return binn(self.ptmin,self.ptmax,self.ymin,self.ymax,self.type,self.i)


    def getweights(self):
        w=[]
        for event in self.evs:
            w.append(e.w)
            
        return w
    
    
    def cost(self,i):
        #gets cost of the ith event in evs
        return self.evs[i].cost
    
    
    def phi(self, i):
        #get phi of the ith event in evs
        return self.evs[i].phi


    #################### Here are the batch methods ####################
       
        
    def Fill(self,TTree):
        '''
        fills the binn with the appropriate events from Ttree, according to edges
        event is a Ttree event
        '''

	#iter thru TTree
	for event in TTree:

            #convert from Ttree event to ev object and place in self.evs, if applicable   
            if self.contains(event):
                self.evs.append(ev(event))
 
    
    def resample(self,n):
        '''
        For bootstrapping: creates a binn of the same size as the original binn by randomly sampling, with replacement, the events of the old binn for a new, distinct data set.
        Note that care is taken to ensure there are the same number of each type of event
        Args: self
	n, int, number of the resampling
	seed, into or none, tells how to seed python random module
        Returns: new binn obj
        '''

        #init new tree
        newbinn=binn(self.ptmin,self.ptmax,self.ymin,self.ymax,self.type,self.i,resamp=n,part=self.part) 
	#n means hists will have different names here

        #if we have both types we use the old method from tree

        if self.type=='both':

            #for getting random indices
            n1=len(self.evs) #even integer, because we fill 2 evs for 1 TTree event
            n2= n1/2 #because n1 is even, this is an integer, and is the number of TTree events present

            #fill new tree with n2 pairs of events by random sampling with replacement
            for i in range(int(n2)):

                # get random index
                index=random.randint(0,int(n2-1))
                #grab pair of events and add to tree
                newbinn.evs.append(self.evs[2*index]) #guaranteed to be a genw event
                newbinn.evs.append(self.evs[2*index+1]) #guaranteed to be a prefsrw event

        #otherwise don't worry about types
        else:
            
	    #size of the resample should be poisson distributed

	    #take the size of binn as expectation value
	    lam=len(self.evs)

	    #size of old binn and size of new binn
	    n1=lam
	    n2=numpy.random.poisson(lam)

	    print 'poisson size is'+str(n2)

            for i in range(n2): #get n2 events

                #get random index
                index=numpy.random.randint(0,n1)
                newbinn.evs.append(self.evs[index])

	    #h1.Draw()
	    #c1.Print('randomtest.pdf')

        #place the resampled binn in the self resamples list
	self.resamples.append(newbinn)

	return
    
    
    def process(self):
        '''
        Processes all the events in the filled bin into histograms, which collect all the info needed to 
        compute the ith coefficient and its error
        
        We want to find two quantities over the many, many data points: A_i and its error, sigma
        The easiest way to do this over many, many root files is using histograms.
        
        For a coefficient related to moment <f>,
        
            <f>= sum(f*weight) / sum (weight)
            
            with error
            
            sigma= sigma (sum (f*weight)^2 , sum (weight)^2)
            
        So we need four histograms, one for each of the sums involved in these calculations
            h_num: for all f*weight quantities
            h_den: for all weight quantities
            h_num_err: for all (f*weight)^2 quantities
            h_den_err: for all (weight)^2 quantities
            
            These have been initiated already, and now need to be filled 
        '''
        
        #iterate thru events

        for e in self.evs:
            
            #control which f to find according to which A_i to find
        
            if self.i==0: #we want to find A0         
                #find f
                f= (1/2)*(1-3*(e.cost)**2)
                          
            elif self.i==1:
                #find f
                f= 2*math.sqrt(1-e.cost**2)*e.cost*math.cos(e.phi)
                          
            elif self.i==2:
                #find f
                f= (1-e.cost**2)*math.cos(2*e.phi)
                
            elif self.i==3:
                f= math.sqrt(1-e.cost**2)*math.cos(e.phi)
                          
            elif self.i==4:
                f= e.cost
                          
            elif self.i==5: 
                f=(1-e.cost**2)*math.sin(2*e.phi)
                          
            elif self.i==6:
                f=2*math.sqrt(1-e.cost**2)*e.cost*math.sin(e.phi)
                          
            elif self.i==7:
                f=math.sqrt(1-e.cost**2)*math.sin(e.phi)
                
            #fill hists
            self.num.Fill(f*e.w/abs(e.w)) #so we need a factor of abs(e.w) later     
            self.den.Fill(e.w/abs(e.w)) #likewise
            self.num_err.Fill(f**2) #so we need a factor of e.w**2
            self.den_err.Fill(1) #likewise

	    #fill MC data hists for A0 binns
	    if self.i == 0:

		#boson angular
		self.phics.Fill(e.phi,e.w)
		self.costcs.Fill(e.cost,e.w)
		self.gen2d.Fill(e.phi,e.cost,e.w)

		#fill lepton kinematics from tuples
		#for eta in e.lep_eta:
		    #self.lep_eta.Fill(eta,e.w)

		#for pt in e.lep_pt:
		    #self.lep_pt.Fill(pt,e.w)

	        for i in range( len( e.lep_eta) ):

	            self.reco2d.Fill(e.lep_eta[i],e.lep_pt[i],e.w )

		#fill dressed lepton kinematics from tuples
		for i in range( len(e.dress_eta)):

		    self.dress2d.Fill(e.dress_eta[i],e.dress_pt[i],e.w)
	
		#boson kinematics
	    	#self.y.Fill(e.y,e.w)
	    	#self.pt.Fill(e.pt,e.w)
	    	self.bos2d.Fill(e.y,e.pt,e.w)


	#recursively process all the resamples
	for re in self.resamples:
	    re.process()


	return None



    def write(self,filename):
	'''
	Writes all of the histograms filled in process() to the new TFile given by filename
	''' 

	#open file
	file1=R.TFile(filename,'update')

	#write hists back to file
	self.num.Write()
	self.den.Write()
	self.num_err.Write()
	self.den_err.Write()

	#write raw data hists as well, if they exist
	if self.i == 0 or self.i ==8:

	    self.phics.Write()
	    self.costcs.Write()
	    self.gen2d.Write()

	    #self.lep_pt.Write()
	    #self.lep_eta.Write()
	    self.reco2d.Write()

	    self.dress2d.Write()

	    #self.pt.Write()
	    #self.y.Write()
	    self.bos2d.Write()



	#close file
	file1.Close()

	#recursively write all resamples
	for re in self.resamples:
	    re.write(filename)


    def write_ttree(self,filename):
	'''
	Writes the events of the binn evs list to a small TTree to be accessed later during smoothing
	'''
	
	#open file
	file1=R.TFile(filename,'update')

	mytree=R.TTree(str(self),str(self) )

	#create arrays for float variables
	pt=array.array('d',[0.])
	y=array.array('d',[0.])
	cost=array.array('d',[0.])
	phi=array.array('d',[0.])
	w=array.array('d',[0.])

	#create arrays for number of leptons, lep_eta and lep_pt (tuples/buffers)
	n_lep=array.array('i',[0])
	lep_eta=array.array('d',50*[0.])
	lep_pt=array.array('d',50*[0.])

	n_dress=array.array('i',[0])
	dress_eta=array.array('d',50*[0.])
	dress_pt=array.array('d',50*[0.])

	#create ttree branches
	mytree.Branch('pt',pt,'pt/D' )
	mytree.Branch('y',y,'y/D' )
	mytree.Branch('cost',cost,'cost/D' )
	mytree.Branch('phi',phi,'phi/D' )
	mytree.Branch('w',w,'w/D' )

	mytree.Branch('n_lep',n_lep,'n_lep/I')
	mytree.Branch('lep_eta',lep_eta,'lep_eta[n_lep]/D' )
	mytree.Branch('lep_pt',lep_pt,'lep_pt[n_lep]/D' )

	mytree.Branch('n_dress',n_dress,'n_dress/I')
	mytree.Branch('dress_eta',dress_eta,'dress_eta[n_lep]/D' )
	mytree.Branch('dress_pt',dress_pt,'dress_pt[n_lep]/D' )


	#fill tree with event attributes
	print 'start events loop'

	for e in self.evs:


	    #float vars
    	    pt[0]=e.pt
    	    y[0]=e.y
	    cost[0]=e.cost
	    phi[0]=e.phi
	    w[0]=e.w

	    #buffer vars
	    n_lep[0]=len(e.lep_eta) #get length of the buffer
	    for j in range(n_lep[0]):
		lep_eta[j]=e.lep_eta[j]
		lep_pt[j]=e.lep_pt[j]
		
	    n_dress[0]=len(e.dress_eta)
	    for j in range(n_dress[0]):
		dress_eta[j]=e.dress_eta[j]
		dress_pt[j]=e.dress_pt[j]

	    #fill the tree
    	    mytree.Fill()



	#write to a .root file
	file1.Write()
	file1.Close()

	#recursively write all resamples
	for re in self.resamples:
	    re.write_ttree(filename)


    #################### Here are the post batch methods ####################


    def gethists(self,batch):
	'''
	takes the histograms associated with this binn out of the TFile filename, then adds them to the 'master'
	hists stored in the binn as attributes
	'''

	print self

	#get list of files to get hists from
	#goodfiles,goodtrees=filecheck.get_good_files(n_trees,batch)

	#iter thru all good files
	if True:

	    filename='/eos/user/c/cbunker/WJet_batch'+str(batch)+'.root'
	    #filename='/afs/cern.ch/work/c/cbunker/WJet_batch'+str(batch)+'.root'
	
	    #open the TFile
	    file1=R.TFile(filename)

	    #get the hists from the tfile
	    #recall that str(self) gives the keyword for accessing the hists
	    num=file1.Get(str(self)+'_num')
	    den=file1.Get(str(self)+'_den')
	    num_err=file1.Get(str(self)+'_num_err')
	    den_err=file1.Get(str(self)+'_den_err')

	    #print str(tree_n)+str(self)+'_num'
	    #print num

	    #add the extracted hists to the self hists
	    self.num.Add(num)
	    self.den.Add(den)
	    self.num_err.Add(num_err)
	    self.den_err.Add(den_err)

	    #also get MC data hists if necessary
	    if self.i == 0:

	    	self.phics.Add(file1.Get(str(self)+'phics'))
	    	self.costcs.Add(file1.Get(str(self)+'costcs'))
	    	self.gen2d.Add(file1.Get(str(self)+'gen2d'))

	    	#self.lep_eta.Add(file1.Get(str(tree_n)+str(self)+'lep_eta'))
	    	#self.lep_pt.Add(file1.Get(str(tree_n)+str(self)+'lep_pt'))
	    	self.reco2d.Add(file1.Get(str(self)+'reco2d'))

		self.dress2d.Add(file1.Get(str(self)+'dress2d'))

	    	#self.pt.Add(file1.Get(str(tree_n)+str(self)+'pt'))
	    	#self.y.Add(file1.Get(str(tree_n)+str(self)+'y'))
	    	self.bos2d.Add(file1.Get(str(self)+'bos2d'))


	    #close file
	    file1.Close()


    def getresamples(self,n_batches,n_resamples):
	'''
	Creates resampled binns corresponding to those created in ang_job, which we now have hists for
	For each of those binns, gets the associated histograms from the Tfile filename and adds them to the attribute hists
	'''

	#iter thru batches
	for batch in range(n_batches):

	    #iter thru resamples in this batch
	    for i in range(1+n_resamples*batch,1+n_resamples*(batch+1)): #needs to start at 1
	    
	        #make newbinn
	        newbinn=binn(self.ptmin,self.ptmax,self.ymin,self.ymax,self.type,self.i,resamp=i) 

	        newbinn.gethists(batch) #stores the resampled data in the binn

	        self.resamples.append(newbinn) #stores the resampled binn in the resamples list attribute

	return None


    def gettree(self,filename):
	'''
	Gets the trimmed tree corresponding to the binn from the .root file filename
	Returns this tree
	'''

	#open the TFile
	file1=R.TFile(filename)

	tree1=file1.Get(str(self))

	return tree1
	


    def calc(self):
	'''
	Takes all the data stored in the histogram attributes to return A_i and the error
	'''

	#abs(event.weightGen)
        weight=225892.453125
	
	##get the num val out of the hist and sum ##
	num=0

	#iter thru hist bins
	for i in range(201):
	    #sum product of each bin, with scale
	    num += self.num.GetBinCenter(i)*self.num.GetBinContent(i)*weight
	

	## similarly for den, num_errs, den errs ##

	den=0
	#sum all the positive weights
	den += weight*self.den.GetBinContent(2)
	
	#sum all the negative weights
	den += (-1)*weight*self.den.GetBinContent(1)

	#now we can determine A 
	A=num/den

	num_err=0
	for i in range(201):
	    num_err += self.num_err.GetBinCenter(i)*self.num_err.GetBinContent(i)*(weight**2)

	den_err=(weight**2)*self.den_err.GetBinContent(1)

	#determine err
	#calculate the error from the seperate sigmas
	sigAB=0 #covariance is zero for now
	err=abs(num/den)*math.sqrt((num_err/num**2) + (den_err/den**2) - 2*sigAB/(num*den))

	#update A and err according to the appropriate constants
	if self.i == 0:
	    A= A*20/3 + 2/3
	    err = err*20/3
	elif self.i == 1:
	    A= A*5
	    err=err*5
	elif self.i == 2:
	    A = A*10
	    err=err*10
	elif self.i == 3:
	    A = A*4
	    err=err*4
	elif self.i == 4:
	    A = A*4
	    err=err*4
	elif self.i == 5:
	    A = A*5
	    err=err*5
	elif self.i == 6:
	    A = A*5
	    err=err*5
	elif self.i == 7:
	    A = A*4
	    err=err*4


	#return tuple of A and error
	return A, err

    def  sumw(self):

        weight=225892.453125 #from the ttree
	den=0
	#sum all the positive weights
	den += weight*self.den.GetBinContent(2)
	
	#sum all the negative weights
	den += weight*self.den.GetBinContent(1)*(-1)

	return den
	

    def rms(self):
        '''
        Computes the rms error on a coefficient over the A's calculated from the resampled counterparts of the binn
	    Args:
	    self, binn obj
        Returns: ntuple of A and rms error on A
        '''

        #calculate Ahat with calc method
        Ahat,dumerr=self.calc()

        #get other vals from resamples
        Avals=[]

	#iter thru resampled binns stored in resamples list
        for b in self.resamples:

	    #get this binn's A val with calc method
	    A,dumerr = b.calc()
	    Avals.append(A)

        #compute error from data
        rms_sum=0
        for val in Avals:
   	    rms_sum += (Ahat - val)**2

        rms=math.sqrt(rms_sum / len(self.resamples))

        return Ahat, rms


    def rms_reco(self,var):
	'''
	Takes all the resampled hists of the reco/raw data (phi, cost, lep_pt, lep_eta) and 
        converts to a single hist with avg and rms error
	Args:
	var, tells which param to rms (phi, cost, pt, eta
	'''

	#tell all var methods whether to fix yhat
	fixed=False

	### for phi ###
	if var=='phi':
	
	    #get the hists we want
	    hlist_phi=[]

	    #go thru resampled binns
	    for b in self.resamples:
	        hlist_phi.append(b.phics)

	    #create master hist
	    nbins_phi=hlist_phi[0].GetXaxis().GetNbins()
	    master_phi=R.TH1F(str(self)+'_rms_phi','Bootstrapped \phi data',nbins_phi,0,2*math.pi)

	    #run thru bins
	    for i in range(1,nbins_phi+1):

	        #this is where we rms process

	        #determine average
		y_hat=self.phics.GetBinContent(i)

	        #determine error
	        rms_sum=0
	        for h in hlist_phi:
		    rms_sum += (y_hat - h.GetBinContent(i) )**2
	        err=math.sqrt(rms_sum / len(hlist_phi) )

	        #place these vals in the master hist
	        master_phi.SetBinContent(i,y_hat)
	        master_phi.SetBinError(i,err)

	    return master_phi
	
	### for cost ###
	elif var=='cost':
	    #get the hists we want
	    hlist_cost=[]

	    #go thru resampled binns
	    for b in self.resamples:
	        hlist_cost.append(b.costcs)

	    #create master hist
	    nbins_cost=hlist_cost[0].GetXaxis().GetNbins()
	    master_cost=R.TH1F(str(self)+'_rms_cost','Bootstrapped cos(theta) data',nbins_cost,-1,1)

	    #run thru bins
	    for i in range(1,nbins_cost+1):

	        #this is where we rms process

	        #determine average
		y_hat=self.costcs.GetBinContent(i)

	        #determine error
	        rms_sum=0
	        for h in hlist_cost:
		    rms_sum += (y_hat - h.GetBinContent(i) )**2
	        err=math.sqrt(rms_sum / len(hlist_cost) )

	        #place these vals in the master hist
	        master_cost.SetBinContent(i,y_hat)
	        master_cost.SetBinError(i,err)
	
	    return master_cost	
	
	### for lep pt ###
	elif var=='pt':
	    #get the hists we want
	    hlist_pt=[]

	    #go thru resampled binns
	    for b in self.resamples:
	        hlist_pt.append(b.lep_pt)

	    #create master hist
	    nbins_pt=hlist_pt[0].GetXaxis().GetNbins()
	    master_pt=R.TH1F(str(self)+'_rms_pt','Bootstrapped lepton p_{T} data',nbins_pt,0,120)

	    #run thru bins
	    for i in range(1,nbins_pt+1):

	        #this is where we rms process

	        #determine average
		y_hat=self.lep_pt.GetBinContent(i)

	        #determine error
	        rms_sum=0
	        for h in hlist_pt:
		    rms_sum += (y_hat - h.GetBinContent(i) )**2
	        err=math.sqrt(rms_sum / len(hlist_pt) )

	        #place these vals in the master hist
	        master_pt.SetBinContent(i,y_hat)
	        master_pt.SetBinError(i,err)

	    return master_pt	
	
	### for lep eta ###
	elif var=='eta':
	    #get the hists we want
	    hlist_eta=[]

	    #go thru resampled binns
	    for b in self.resamples:
	        hlist_eta.append(b.lep_eta)

	    #create master hist
	    nbins_eta=hlist_eta[0].GetXaxis().GetNbins()
	    master_eta=R.TH1F(str(self)+'_rms_eta','Bootstrapped lepton \eta data',nbins_eta,-3,3)

	    #run thru bins
	    for i in range(1,nbins_eta+1):

	        #this is where we rms process

	        #determine average
		y_hat=self.lep_eta.GetBinContent(i)

	        #determine error
	        rms_sum=0
	        for h in hlist_eta:
		    rms_sum += (y_hat - h.GetBinContent(i) )**2
	        err=math.sqrt(rms_sum / (len(hlist_eta) - 1) )

	        #place these vals in the master hist
	        master_eta.SetBinContent(i,y_hat)
	        master_eta.SetBinError(i,err)

	    return master_eta

	### for gen2d ###
	elif var=='gen2d':

	    #get the hists we want
	    hlist_gen=[]

	    #go thru resampled binns
	    for b in self.resamples:
	        hlist_gen.append(b.gen2d)

	    #create master hist
	    nbins_x=hlist_gen[0].GetXaxis().GetNbins()
	    nbins_y=hlist_gen[0].GetYaxis().GetNbins()
	    master_gen=R.TH2F(str(self)+'gen2d_rms','Bootstrapped \phi and cos(theta) data',nbins_x,0,2*math.pi,nbins_y,-1,1)

	    #run thru bins in x
	    for i in range(1,nbins_x+1):

	        #run thru bins in y
	        for j in range(1,nbins_y+1):

	            #this is where we rms process
		    numbers=[]
	            #determine average  by fixed or unfixed method
		    if fixed:
		        y_hat=self.gen2d.GetBinContent(i,j)

		    else:
		    	hat_sum=0
		    	for h in hlist_gen:
			    hat_sum += h.GetBinContent(i,j)
		    	y_hat=hat_sum/len(hlist_gen)

	            #determine error
	            rms_sum=0
	            for h in hlist_gen:
		        rms_sum += (y_hat - h.GetBinContent(i,j) )**2

	            err=math.sqrt(rms_sum / (len(hlist_gen) - 1) )

	            #place these vals in the master hist
	            master_gen.SetBinContent(i,j,y_hat)
	            master_gen.SetBinError(i,j,err)

	    return master_gen

	### for reco2d ###
	elif var=='reco2d':

	    print 'rms reco'

	    #get the hists we want
	    hlist_reco=[]

	    #go thru resampled binns
	    for b in self.resamples:
	        hlist_reco.append(b.reco2d)

	    #create master hist
	    nbins_x=hlist_reco[0].GetXaxis().GetNbins()
	    nbins_y=hlist_reco[0].GetYaxis().GetNbins()
	    master_reco=R.TH2F(str(self)+'reco2d_rms','Bootstrapped lepton \eta and p_{T}  data',nbins_x,-2.4,2.4,nbins_y,0,120)

	    #run thru bins in x
	    for i in range(1,nbins_x+1):

	        #run thru bins in y
	        for j in range(1,nbins_y+1):

	            #determine average  by fixed or unfixed method
		    if fixed:
		        y_hat=self.reco2d.GetBinContent(i,j)

		    else:
		    	hat_sum=0
		    	for h in hlist_reco:
			    hat_sum += h.GetBinContent(i,j)
		    	y_hat=hat_sum/len(hlist_reco)

	            #determine error
	            rms_sum=0
	            for h in hlist_reco:
		        rms_sum += (y_hat - h.GetBinContent(i,j) )**2

	            err=math.sqrt(rms_sum / (len(hlist_reco) - 1) )

	            #place these vals in the master hist
	            master_reco.SetBinContent(i,j,y_hat)
	            master_reco.SetBinError(i,j,err)

	    return master_reco

	elif var=='dress2d':
	    #get the hists we want
	    hlist_dress=[]

	    #go thru resampled binns
	    for b in self.resamples:
	        hlist_dress.append(b.dress2d)

	    #create master hist
	    nbins_x=hlist_dress[0].GetXaxis().GetNbins()
	    nbins_y=hlist_dress[0].GetYaxis().GetNbins()
	    master_dress=R.TH2F(str(self)+'dress2d_rms','Bootstrapped dressed lepton \eta and p_{T}  data',nbins_x,-2.4,2.4,nbins_y,0,120)

	    #run thru bins in x
	    for i in range(1,nbins_x+1):

	        #run thru bins in y
	        for j in range(1,nbins_y+1):

	            #this is where we rms process

	            #determine average  by fixed or unfixed method
		    if fixed:
		        y_hat=self.dress2d.GetBinContent(i,j)

		    else:
		    	hat_sum=0
		    	for h in hlist_dress:
			    hat_sum += h.GetBinContent(i,j)
		    	y_hat=hat_sum/len(hlist_dress)

	            #determine error
	            rms_sum=0
	            for h in hlist_dress:
		        rms_sum += (y_hat - h.GetBinContent(i,j) )**2

	            err=math.sqrt(rms_sum / (len(hlist_dress) - 1) )

	            #place these vals in the master hist
	            master_dress.SetBinContent(i,j,y_hat)
	            master_dress.SetBinError(i,j,err)

	    return master_dress



	### for bos2d ###
	elif var=='bos2d':
	    #get the hists we want
	    hlist_bos=[]

	    #go thru resampled binns
	    for b in self.resamples:
	        hlist_bos.append(b.bos2d)

	    #create master hist
	    nbins_x=hlist_bos[0].GetXaxis().GetNbins()
	    nbins_y=hlist_bos[0].GetYaxis().GetNbins()
	    master_bos=R.TH2F(str(self)+'bos2d_rms','Bootstrapped boson y and p_{T}  data',nbins_x,-4,4,nbins_y,0,120)

	    #run thru bins in x
	    for i in range(1,nbins_x+1):

	        #run thru bins in y
	        for j in range(1,nbins_y+1):

	            #this is where we rms process

	            #determine average  by fixed or unfixed method
		    if fixed:
		        y_hat=self.bos2d.GetBinContent(i,j)

		    else:
		    	hat_sum=0
		    	for h in hlist_bos:
			    hat_sum += h.GetBinContent(i,j)
		    	y_hat=hat_sum/len(hlist_bos)

	            #determine error
	            rms_sum=0
	            for h in hlist_bos:
		        rms_sum += (y_hat - h.GetBinContent(i,j) )**2

	            err=math.sqrt(rms_sum / (len(hlist_bos) - 1) )

	            #place these vals in the master hist
	            master_bos.SetBinContent(i,j,y_hat)
	            master_bos.SetBinError(i,j,err)

	    return master_bos

	### for bos2d ###
	elif var=='dress2d':
	    #get the hists we want
	    hlist_dress=[]

	    #go thru resampled binns
	    for b in self.resamples:
	        hlist_dress.append(b.dress2d)

	    #create master hist
	    nbins_x=hlist_dress[0].GetXaxis().GetNbins()
	    nbins_y=hlist_dress[0].GetYaxis().GetNbins()
	    master_dress=R.TH2F(str(self)+'dress2d_rms','Dress lepton \eta and p_{T} data',nbins_x,-2.4,2.4,nbins_y,0,120)

	    #run thru bins in x
	    for i in range(1,nbins_x+1):

	        #run thru bins in y
	        for j in range(1,nbins_y+1):

	            #this is where we rms process

	            #determine average  by fixed or unfixed method
		    if fixed:
		        y_hat=self.dress2d.GetBinContent(i,j)

		    else:
		    	hat_sum=0
		    	for h in hlist_dress:
			    hat_sum += h.GetBinContent(i,j)
		    	y_hat=hat_sum/len(hlist_dress)

	            #determine error
	            rms_sum=0
	            for h in hlist_dress:
		        rms_sum += (y_hat - h.GetBinContent(i,j) )**2

	            err=math.sqrt(rms_sum / (len(hlist_dress) - 1) )

	            #place these vals in the master hist
	            master_dress.SetBinContent(i,j,y_hat)
	            master_dress.SetBinError(i,j,err)

	    return master_dress



class holder(object):
    '''
    Container object for 3d array of binns
	binns need to be arranged by A_i, y and pt, thus storing them ins 3d
	this object makes the storage process easier
	initing this object also creates the relevant binns
    '''

    #################### basic methods ####################

    def __init__(self,ptedges,yedges,empty=False,part=9999):
	'''
	attributes:
	main, 3d array of binns
	Args
	ptedges, y edges, lists that define binning
	empty, tells whether to actually fill with binn objects

	'''
	#init list
    	self.main=[]
	self.ptedges=ptedges
	self.yedges=yedges

	#go through y edges and add to main list
	for j in range(len(yedges)-1):
	
	    ptlist=[]

    	    #go thru pt edges and add to pt list
    	    for k in range(len(ptedges) - 1):

		ilist=[]
	
		#go thru i's and add to i list
		for i in range(8):

		    #check if we want to actually fill in binns
		    if not empty:
	    	        ilist.append(binn(ptedges[k],ptedges[k+1],yedges[j],yedges[j+1],'prefsrw',i,part=part))

		#add to pt list
		ptlist.append(ilist)
	    #add to main list
	    self.main.append(ptlist)

    def forin(self):
	'''
	returns contents of the self.main 3d list as a simple list for traversing with a for loop
	'''
	forlist=[]
	
	#iter thru all contents
	for ptlist in self.main:
	    for bgroup in ptlist:
		for b in bgroup:
		    forlist.append(b)

	return forlist

    def bgroups(self):
	'''
	returns bgroups of the self.main 3d list as a simple list for traversing with a for loop
	'''
	forlist=[]
	
	#iter thru all contents
	for ptlist in self.main:
	    for ilist in ptlist:
		forlist.append(ilist)

	return forlist

    def get_As(self,i=0):
	'''
	Return only binns with binn.i == i
	defaults to getting A0
	'''
	Alist=[]
	for b in self.forin():
	    if b.i == i:
		Alist.append(b)

	return Alist

    def Fill(self,TTree):
        '''
        fills the binn holder with the appropriate events from Ttree, according to edges
        event is a Ttree event
        '''

	#iter thru TTree
	for event in TTree:

            #convert from Ttree event to ev object and place in appropriate binns   
	    #iter thru binns
	    for b in self.forin():
		if b.contains(event):
		    b.evs.append(ev(event))

	return


    def write_ttree(self,filename):
	'''
	Writes all A0 binn events to ttrees stored in the root file filename
	See the binn method write_ttree
	'''

	for b in self.get_As():
	    b.write_ttree(filename)




    def resampled_holders(self):
	'''
	Makes holders out of all the sets of resampled binns
	Returns a list of holders, each containing resampled binns
	'''

	#return variable
	resampled_holders=[]
	
	#how many to make
	n_resamples=len(self.forin()[0].resamples)

	#make them
	for n in range(n_resamples):

	    #new holder
	    new=holder(self.ptedges,self.yedges,empty=True)

	    #reset the old main list
	    new.main=[]

	    #refill main list from self main list

	    #run through y edges
	    for j in range(len(self.yedges)-1):
	
	        ptlist=[]

    	        #run thru pt edges and add to pt list
    	        for k in range(len(self.ptedges) - 1):

		    ilist=[]
	
		    #go thru binns of A_i and add to i list
		    for i in range(8):
			#grab appropriate resample from appropriate binn
	    	        ilist.append( self.main[j][k][i].resamples[n] )

		    #add to pt list
		    ptlist.append(ilist)

	        #add to main list
	        new.main.append(ptlist)

	    #add the new holder to the return object
	    resampled_holders.append(new)

	return resampled_holders


    def fill_smooth(self,n_batches,n_resamples,ptedges=[],yedges=[]):
	'''
	fills the associated reweightor with events from the trimmed ttrees
	reimplemented from code previously in bholder.comp_all()

	Args:
	parts, part #'s of trees to fill from
	n_files, #of files to get trees from
	edges, tells what binns to get trees from, defaults to 1 binn for one binn holders
	Returns
	A filled reweightor object
	
	'''

	#if empty lists are passed to edges (default option) we actually get them from self
	if len(ptedges)==0 and len(yedges)==0:
	    ptedges=self.ptedges
	    yedges=self.yedges

	#get reweights
	reweights=self.smooth()


	#batch is 0 for the raw data
	#iter over all trees to fill reweight binns
	batch=0
	if True:
	    if True:

		filename='/eos/user/c/cbunker/WJet_batch'+str(batch)+'.root'
		#filename='/afs/cern.ch/work/c/cbunker/WJet_batch'+str(batch)+'.root'

	        #iter over A0 binns in bholder
	        n_bins=(len(ptedges) - 1)*(len(yedges) - 1)
	        for j in range( n_bins ):

		    #get edges corresponding to this specific binn
		    new_pt, new_y = change_edges( ptedges, yedges, j )
		    #this step should default to just picking out the one A0 binn in the holder

	            #treename varies by part
		    #treename has form 0A_0_1_0_1_r1
	            treename='A0_'+str(new_pt[0])+'_'+str(new_pt[1])+'_'
		    treename += str(new_y[0])+'_'+str(new_y[1])
		    print 'filename is ', filename
		    print 'tree name is ', treename
	            reweights.make_smoothed(filename,treename) #fills the reweightor hists

	#get resample reweights, and fill in to resamples list of reweightor binns

	#fill each reweightor corresponding to each resample
	for i in range(len( self.resampled_holders() )): #iters over resamples

	    #make the reweightor
	    my_reweight = self.resampled_holders()[i].smooth() 

	    #for debugging, find the average reweight
	    print 'resample '+str(i+1)+' has average '+str(my_reweight.avg() )  


	    #determine which batch this resamples is in
	    batch= int( math.floor( i/n_resamples ) )
	    if True:
	        if True:

		    filename='/eos/user/c/cbunker/WJet_batch'+str(batch)+'.root'
		    #filename='/afs/cern.ch/work/c/cbunker/WJet_batch'+str(batch)+'.root'

		    #getresi = 100*batch + (i-n_resamples*batch)

	            #iter over A0 binns in bholder
		    n_bins=(len(ptedges) - 1)*(len(yedges) - 1)
	    	    for j in range( n_bins ):

		        #get edges corresponding to this specific binn
		        new_pt, new_y = change_edges( ptedges, yedges, j )

		        #this step should default to just picking out the one A0 binn in the holder
		        #print new_pt, new_y
  
	                #treename varies by part
	                treename='A0_'+str(new_pt[0])+'_'+str(new_pt[1])+'_'
		        treename += str(new_y[0])+'_'+str(new_y[1])+'_r'+str(i+1)
		        print 'filename is ', filename
		        print 'tree name is ', treename
	                my_reweight.make_smoothed(filename,treename) #fills the reweightor hists for this binn

	    #now add all these filled, resample binns to the 
	    #resamples list of the appropriate smoothed binn
	    for j in range (len ( reweights.forin()) ):
		
		#puts the jth filled resampled reweight binn in the resamples list of the jth filled reweight binn
		reweights.forin()[j].resamples.append( my_reweight.forin()[j] )

	return reweights


    def make_inclusive(self):
	'''
	Takes the data spread out between binns in the A0 hists and combines, as well as combines all the resamples across bins
	'''

	#combine all the A0 hists
	#they all go in the A0 hists of the first bin
	for i in range(1,len(self.get_As() ) ):

	    #add gen2d param hists
	    self.get_As()[0].costcs.Add( self.get_As()[i].costcs )
	    self.get_As()[0].phics.Add( self.get_As()[i].phics )
	    self.get_As()[0].gen2d.Add( self.get_As()[i].gen2d )

	    #add reco2d param hists
	    self.get_As()[0].reco2d.Add( self.get_As()[i].reco2d )

	    #add bos2d param hists
	    self.get_As()[0].bos2d.Add( self.get_As()[i].bos2d )

	    #add dress2d param hists
	    self.get_As()[0].dress2d.Add( self.get_As()[i].dress2d )

	#combine all the A_j hists
	for j in range(8 ):

	    #again iter thru all pT, y binns
	    for k in range(1,len(self.get_As() ) ):

	    	self.get_As(i=j)[0].num.Add( self.get_As(i=j)[k].num )
	    	self.get_As(i=j)[0].num_err.Add( self.get_As(i=j)[k].num_err )
	    	self.get_As(i=j)[0].den.Add( self.get_As(i=j)[k].den )
	    	self.get_As(i=j)[0].den_err.Add( self.get_As(i=j)[k].den_err )

	#now do the resamples
	#we want to put all the resamples across parts into the resamples list of the first pt, y binn
	for j in range(8 ):

	    #again iter thru all pT, y binns
	    for k in range(1,len(self.get_As() ) ):

		#grab the resamples list of the kth binn
	    	for resample in self.get_As(i=j)[k].resamples :

		    #place each resampled binn in the 0th binn list
		    self.get_As(i=j)[0].resamples.append( resample ) 

	return 




    #################### plotting methods ####################


    def plot_all(self, i, unroll=False):  
        '''
        Goal is graphing Ai versus pT for all y binns, except now we find RMS errors also via resampling
	    Args:
            bholder, holder obj which contains 3d list of all the binns needed
	    i, int, which coefficient A_i to plot
	    both, bool, if true plots bootstrapped errors as well

        '''
        #create filename
        stringA='A_{'+str(i)+'}'
        filename='A'+str(i)+'_plots_rms.pdf'

	#if we want to unroll, create the arrays
	if unroll:
	    xvals=array.array('f')
	    xerrs=array.array('f')
	    normvals=array.array('f')
	    rmsvals=array.array('f')
	    normerrs=array.array('f')
	    rmserrs=array.array('f')
	    normrels=array.array('f')
	    rmsrels=array.array('f')

	#otherwise we need canvases:
	else:
            #open file, keeping open
            dum_canvas=R.TCanvas('dum','title',550,800)
            dum_canvas.SetFillColor(R.kGray)
            dum_canvas.Print(filename+'(')

        #iter thru y bins, making plots for each 
	pt_list_counter = -1  
        for ptlist in self.main: #grab list of pt binns from 3d array

	    pt_list_counter += 1
		
            #create title for this plot
            title=stringA+' vs p_{T} for prefsrw data, '+str(ptlist[0][0].ymin)+' < |y| < '+str(ptlist[0][0].ymax)+'; p_{T} [GeV] ;'+stringA

	    #canvas
	    c1=R.TCanvas('c1','title',550,800)

	    #draw top plot: coefs with error bars

	    #upper tpad
	    p1=R.TPad('p1','p1',0,0.35,1,1)
	    p1.SetGrid()
	    p1.Draw()
	    p1.cd()

	    #genw with normal errors (is an Agraph object now!)
            norm=plot(ptlist,i,title)

	    #genw with resampled rms errors
	    rms=plot_resampled(ptlist,i,title)

	    #draw graphs
	    rms.graph.Draw('ap2')
	    norm.graph.Draw('psame')

   	    #make and format legend
	    leg1=R.TLegend(0.7,0.2,1,0.4) #position legend
	    leg1.AddEntry(norm.graph,'Original error')
	    leg1.AddEntry(rms.graph,'RMS error, '+str(len(ptlist[0][0].resamples))+' resamples')
	    leg1.SetTextSize(0.02)
	    leg1.Draw()

	    #lower plot: percentage errors
	    c1.cd()

	    #lower tpad
	    p2=R.TPad('p2','p2',0,0,1,0.35)
	    p2.SetGrid()
	    p2.Draw()
   	    p2.cd()   
	
	    #make relative error graphs
	    norm_g2=R.TGraph(len(norm.xvals),norm.xvals,norm.yerrs)
	    rms_g2=R.TGraph(len(rms.xvals),rms.xvals,rms.yerrs)

	    #make legend
	    leg2=R.TLegend(0.7,0.7,1,1)
	    leg2.AddEntry(norm_g2,'Original percent error')
	    leg2.AddEntry(rms_g2,'RMS percent error, '+str(len(ptlist[0][0].resamples))+' resamples')

	    #make multigraph
	    mg2=R.TMultiGraph()
	    mg2.Add(norm_g2)
	    mg2.Add(rms_g2)
	    mg2.SetTitle('Percent error on '+stringA+'; p_{T} [GeV]; (\sigma_{'+stringA+'} / '+stringA+') x 100')
	    mg2.Draw('al')
	    mg2.GetXaxis().SetLimits(rms.xvals[0]-rms.xerrs[0],rms.xvals[-1]+rms.xerrs[-1])
	    leg2.Draw()

	    #print canvas or fill unrolled arrays as required
	    if unroll:
		
		#iter thru both Agraph objects to fill unrolled arrays
		for i in range( len( norm.xvals ) ):
		    #xvals.append(norm.xvals[i])
		    #update xval for each slice
		    xvals.append( norm.xvals[i]+50*pt_list_counter )
		    xerrs.append( 0 )
		    normvals.append( norm.yvals[i] )
		    rmsvals.append( rms.yvals[i] )
		    normerrs.append( norm.yerrs[i] )
		    rmserrs.append( rms.yerrs[i] )
		    normrels.append( norm.get_rels()[i] )
		    rmsrels.append( rms.get_rels()[i] )
		
	    else:
	        #save canvas
	        c1.Print(filename)

	#make and return graphs as desired
	if unroll:
	    
	    #create the A_i graphs
	    normgraph=R.TGraphErrors(len(xvals),xvals,normvals,xerrs,normerrs)
	    rmsgraph=R.TGraphErrors(len(xvals),xvals,rmsvals,xerrs,rmserrs)

	    #create ratio graphs
	    normratio=R.TGraph(len(xvals),xvals,normrels)
	    rmsratio=R.TGraph(len(xvals),xvals,rmsrels)

	    return normgraph, rmsgraph, normratio, rmsratio

	else:
        
            #close file
            dum_canvas.Print(filename+')')

            return None


    def plot_norm(self, i):  
        '''
        Goal is graphing Ai versus pT for all y binns, with only normal errors
	    Args:
            bholder, holder obj which contains 3d list of all the binns needed
	    i, int, which coefficient A_i to plot
	    both, bool, if true plots bootstrapped errors as well

        '''
        #create filename
        stringA='A_{'+str(i)+'}'
        filename='A'+str(i)+'_plots_norm.pdf'

        #open file, keeping open
        dum_canvas=R.TCanvas('dum','title',550,800)
        dum_canvas.SetFillColor(R.kGray)
        dum_canvas.Print(filename+'(')

        #iter thru y bins, making plots for each   
        for ptlist in self.main: #grab list of pt binns from 3d array
		
            #create title for this plot
            title=stringA+' vs p_{T} for prefsrw data, '+str(ptlist[0][0].ymin)+' < |y| < '+str(ptlist[0][0].ymax)+'; p_{T} [GeV] ;'+stringA

	    #canvas
	    c1=R.TCanvas('c1','title',550,800)

	    #draw top plot: coefs with error bars

	    #upper tpad
	    p1=R.TPad('p1','p1',0,0.35,1,1)
	    p1.SetGrid()
	    p1.Draw()
	    p1.cd()

	    #genw with normal errors (is an Agraph object now!)
            norm=plot(ptlist,i,title)
	    norm.graph.GetXaxis().SetLimits(norm.xvals[0]-5,norm.xvals[-1]+5)

	    #draw graphs
	    norm.graph.Draw('ap')

   	    #make and format legend
	    leg1=R.TLegend(0.7,0.2,1,0.4) #position legend
	    leg1.AddEntry(norm.graph,'Original error')
	    leg1.SetTextSize(0.02)
	    leg1.Draw()

	    #lower plot: percentage errors
	    c1.cd()

	    #lower tpad
	    p2=R.TPad('p2','p2',0,0,1,0.35)
	    p2.SetGrid()
	    p2.Draw()
   	    p2.cd()   
	
	    #make relative error graphs
	    norm_g2=R.TGraph(len(norm.xvals),norm.xvals,norm.get_rels())
	    norm_g2.SetLineColor(R.kRed)
	    norm_g2.SetFillColor(R.kWhite)

	    #make legend
	    leg2=R.TLegend(0.7,0.7,1,1)
	    leg2.AddEntry(norm_g2,'Original percent error')

	    #make multigraph
	    mg2=R.TMultiGraph()
	    mg2.Add(norm_g2)
	    mg2.SetTitle('Percent error on '+stringA+'; p_{T} [GeV]; (\sigma_{'+stringA+'} / '+stringA+') x 100')
	    mg2.Draw('al')
	    mg2.GetXaxis().SetLimits(norm.xvals[0]-5,norm.xvals[-1]+5)
	    leg2.Draw()


	    #save canvas
	    c1.Print(filename)
        
        #close file
        dum_canvas.Print(filename+')')

        return None

    #################### Plotting or smoothing cost or phi seperately ####################


    def smooth_cost(self, plot=False):
	'''
	A function to compute the appropriate reweights, for smoothing the raw cost data to the smooth curve {}

	This is done for each numerical bin of y and pT, across the bgroup 
	( where a bgroup consists of 8 binns of same numerical edges, each corresponding to a different coefficient)
	Within the bgroup a different reweight is determined for each bin of cost, these are returned together as an array

	The args are the bgroups within self.main, and these should already be filled
	The arg plot, boolean, tells whether to make hist/smooth function comparison plots as well

	The final return product is a matrix of reweight arrays, one for each binning of y, pt. The reweights can then be reapplied to the raw data.
	if plot=True, simply returns multigraphs for graphing instead 
	'''
	
	#container variable for return
	all_reweights=[] #will be filled with reweight arrays, one for each bgroup
	
	#for case plot=True
	all_graphs=[] #will be filled with raw and smooth graphs for each bgroup

	#iter through all the bgroups
	for ptlist in self.main:

	    #first round of sub lists for return objects
	    pt_reweights=[]
	    pt_graphs=[]

	    #continue to get bgroups themselves, and iter thru
	    for bgroup in ptlist:

		#determine the angular coefficients for this bgroup
		A=[]
		for b in bgroup: #they are order b.i=0 .. b.i=7
		    A.append(b.calc()[0]) #the 0th element is the ang coefficient

		#get the raw data histogram and its characteristics
		h1=bgroup[0].costcs
		nbins=h1.GetNbinsX()
		sumw=bgroup[0].sumw()

		#form the TGraph for comparison

		#value holder arrays
		costvals=array.array('f')
		fvals=array.array('f') #corresponds to smooth function
		reweights=[]

		#iter over cost vals
		for i in range(1,nbins+1): #n x vals to correspond to n hist bins, hist binning starts at 1 so do likewise

		    #determine cost from x axis of hist
    		    cost= h1.GetBinCenter(i)

		    #determine f from {} definition
    		    f= 1 + cost**2 + (1/2)*A[0]*(1-3*cost**2) + A[4]*cost

		    #get h val from hist
		    h=h1.GetBinContent(i)

    		    #scale based on the total cross section (sum of weights)
		    f=f*(3/8)*sumw*(2/nbins) #last term is bin width

		    #update arrays
    		    costvals.append(cost)
    		    fvals.append(f)

		    #determine the reweights
		    reweights.append(f/h)

		#put this reweights array in the sub list
		pt_reweights.append(reweights)

		#make mg, if required
		if plot:
		    
		    #make graphs
		    g1=R.TGraph(len(costvals),costvals,fvals)

		    #format graph and hist
		    g1.SetLineColor(R.kRed)
		    g1.SetTitle('Smooth function')
		    g1.SetName('Smooth function')
		    h1.SetLineColor(R.kBlue)
		    h1.SetTitle('Raw cos(theta) data')
		    h1.SetName('Raw cos(theta) data')

		    #add graphs to list as sublist
		    pt_graphs.append([g1,h1])


	    #add the sub lists into the main lists
	    all_reweights.append(pt_reweights)
	    all_graphs.append(pt_graphs)


	#return as requested
	if plot:
	    return all_graphs

	return all_reweights



    def smooth_phi(self, plot=False):
	'''
	See smooth_cost above. Same thing, but with phi

	A function to compute the appropriate reweights, for smoothing the raw phi data to the smooth curve {}

	This is done for each numerical bin of y and pT, across the bgroup 
	( where a bgroup consists of 8 binns of same numerical edges, each corresponding to a different coefficient)
	Within the bgroup a different reweight is determined for each bin of cost, these are returned together as an array

	The args are the bgroups within self.main, and these should already be filled
	The arg plot, boolean, tells whether to make hist/smooth function comparison plots as well

	The final return product is a matrix of reweight arrays, one for each binning of y, pt. The reweights can then be reapplied to the raw data.
	if plot=True, simply returns multigraphs for graphing instead 
	'''
	
	#container variable for return
	all_reweights=[] #will be filled with reweight arrays, one for each bgroup
	
	#for case plot=True
	all_graphs=[] #will be filled with comparison multigraphs, one for each bgroup

	#iter through all the bgroups
	for ptlist in self.main:

	    #first round of sub lists for return objects
	    pt_reweights=[]
	    pt_graphs=[]

	    #continue to get bgroups themselves, and iter thru
	    for bgroup in ptlist:

		#determine the angular coefficients for this bgroup
		A=[]
		for b in bgroup: #they are order b.i=0 .. b.i=7
		    A.append(b.calc()[0]) #the 0th element is the ang coefficient

		#get the raw data histogram and its characteristics
		h1=bgroup[0].phics
		nbins=h1.GetNbinsX()
		sumw=bgroup[0].sumw()

		#form the TGraph for comparison

		#value holder arrays
		phivals=array.array('f')
		fvals=array.array('f') #corresponds to smooth function
		reweights=[]

		#iter over cost vals
		for i in range(1,nbins+1): 
		#n x vals to correspond to n hist bins, hist binning starts at 1 so do likewise

		    #determine cost from x axis of hist
    		    phi= h1.GetBinCenter(i)

		    #determine f from {} definition
    		    f=1+(1/4)*math.cos(2*phi)*A[2]+(3*math.pi/16)*math.cos(phi)*A[3]+(1/2)*math.sin(2*phi)*A[5]+(3*math.pi/16)*math.sin(phi)*A[7]


		    #get h val from hist
		    h=h1.GetBinContent(i)

    		    #scale by unpolarised cross section (sum w)
		    f=f*(1/(2*math.pi))*sumw*(2*math.pi/nbins) #last term is bin width

		    #update arrays
    		    phivals.append(phi)
    		    fvals.append(f)

		    #determine the reweights
		    reweights.append(f/h)

		#put this reweights array in the sub list
		pt_reweights.append(reweights)

		#make mg, if required
		if plot:
		    
		    #make graphs
		    g1=R.TGraph(len(phivals),phivals,fvals)

		    #format graphs
		    g1.SetLineColor(R.kRed)
		    g1.SetTitle('Smooth function')
		    g1.SetName('Smooth function')
		    h1.SetLineColor(R.kBlue)
		    h1.SetTitle('Raw \phi data')
		    h1.SetName('Raw \phi data')

		    #add graphs
		    pt_graphs.append([g1,h1])


	    #add the sub lists into the main lists
	    all_reweights.append(pt_reweights)
	    all_graphs.append(pt_graphs)


	#return as requested
	if plot:
	    return all_graphs

	return all_reweights


    def plot_phi(self):
	'''
	Makes comparison plots of raw data vs smooth function

	Gets multigraphs from self.smooth_phi(plot=True) and then makes and stores plots
	'''

	#get mgs
	all_graphs=self.smooth_phi(plot=True)

	#determine filename
	filename='phi_comp_plots.pdf'

	#begin file with dummy canvas
	dum_c=R.TCanvas()
	dum_c.SetFillColor(R.kGray)
	dum_c.Print(filename+'(') #keep file open

	#iter thru mgs
	for i in range(len(all_graphs)):
	
	    #these mgs all correspond to bgroups of constant y
	    stringy=str(self.main[i][0][0].ymin)+' > |y| > '+str(self.main[i][0][0].ymax) #grabs the first pt binn, A0 binn, of this pt column

	    #keep itering through
	    for j in range(len(all_graphs[i])):

		#these mgs correspond to bgroups with constant pt
		stringpt=str(self.main[i][j][0].ptmin)+' > p_{T} > '+str(self.main[i][j][0].ptmax)

		#get graph and hist
		g1=all_graphs[i][j][0]
		h1=all_graphs[i][j][1]

		#make legend
		leg1=R.TLegend()
		leg1.AddEntry(g1,g1.GetTitle())
		leg1.AddEntry(h1,h1.GetTitle())
		
		#draw
		c1=R.TCanvas()
		h1.Draw('e')
		g1.Draw('same')
		leg1.Draw()
		c1.Print(filename)


	#close the file with dummy canvas
	dum_c.Print(filename+')')
		   	
	return None




    def plot_cost(self):
	'''
	Makes comparison plots of raw data vs smooth function

	Gets graphs and hists from self.smooth_cost(plot=True) and then makes and stores plots
	'''

	#get mgs
	all_graphs=self.smooth_cost(plot=True)

	#determine filename
	filename='cost_comp_plots.pdf'

	#begin file with dummy canvas
	dum_c=R.TCanvas()
	dum_c.SetFillColor(R.kGray)
	dum_c.Print(filename+'(') #keep file open

	#iter thru mgs
	for i in range(len(all_graphs)):
	
	    #these all correspond to bgroups of constant y
	    stringy=str(self.main[i][0][0].ymin)+' > |y| > '+str(self.main[i][0][0].ymax) 
	    #grabs the first pt binn, A0 binn, of this pt column

	    #keep itering through
	    for j in range(len(all_graphs[i])):

		#these graphs correspond to bgroups with constant pt
		stringpt=str(self.main[i][j][0].ptmin)+' > p_{T} > '+str(self.main[i][j][0].ptmax)

		#get graph and hist
		g1=all_graphs[i][j][0]
		h1=all_graphs[i][j][1]

		#make legend
		leg1=R.TLegend()
		leg1.AddEntry(g1,g1.GetTitle())
		leg1.AddEntry(h1,h1.GetTitle())
		
		#draw
		c1=R.TCanvas()
		h1.Draw('e')
		g1.Draw('same')
		leg1.Draw()
		c1.Print(filename)


	#close the file with dummy canvas
	dum_c.Print(filename+')')
		   
	return None




    def get_distribution(self,var,bgroup_index,n_resamples,nbins ):
	'''
	Looks at the distribution of resampled cross section w.r.t a specified parameter.
	Ideally the distribution should be Gaussian
	Args: var, tells parameter to choose bins for
	'''

	#pick the binn we want
	mybinn = self.get_As()[bgroup_index]

	#get list of 1d histograms w.r.t var
	hlist=[]
	for b in mybinn.resamples: #iter thru resamples

	    #grab hist according to var
	    #also define other var-based variables here
	    if var=='phi':
		hlist.append( b.gen2d.ProjectionX(str(b),1,20,'e') )
		low,hi = 0,2*math.pi
		binsize=(hi-low)/20
	    elif var=='cost':
		hlist.append( b.gen2d.ProjectionY(str(b),1,20,'e') )
		low,hi = -1,1
		binsize=(hi-low)/20
	    elif var=='eta':
		hlist.append( b.reco2d.ProjectionX(str(b),1,20,'e') )
		low,hi = -2.4,2.4
		binsize=(hi-low)/20
	    elif var=='pt':
		hlist.append( b.reco2d.ProjectionY(str(b),1,20,'e') )
		low,hi = 0,60
		binsize=(hi-low)/20

	#make a TH2 of the resample distribution
	# x axis is resample #, y axis is var bin
	distribution=R.TH2F(var,var,n_resamples,1,n_resamples+1,nbins,1,nbins+1)

	#format
	distribution.GetXaxis().SetTitle(var)
	distribution.GetYaxis().SetTitle('resample number')
	distribution.GetYaxis().SetNdivisions(n_resamples)

	#iter thru var bins
	for resample in range(1,1+n_resamples):

	    #iter thru resamples
	    for bin in range(1,1+nbins):
	    
		#now fill the distribution th2
		distribution.SetBinContent( resample,bin, hlist[resample-1].GetBinContent(bin) )

	return distribution


    def plot_distribution(self,var,bgroup_index,n_resamples,nbins=10, choosebin=False ):
	'''
	After the resample distribution for a given var is stored in a TH2F by get_distribution(),
	makes an unrolled plot of the xsec of the distribution unrolled along the var distribution
	'''

	#get the distribution
	distribution = self.get_distribution(var,bgroup_index,n_resamples,nbins)

	#c1=R.TCanvas()
	#c1.SetGrid()
	#distribution.Draw()
	#c1.Print('gauss.pdf')

	#find the maximum xsec for each binn so that we can plot in terms of relative xsec
	maxes=[]
	for bin in range(1,1+nbins): #iter thru bins
	    xsec_vals=[]
	    for resample in range(1,1+n_resamples): #iter thru resamples
		xsec_vals.append( distribution.GetBinContent(resample,bin) )

	    #find the max of the resample xsecs in this bin
	    maxes.append( max( xsec_vals ) )
	

	#new hist with xsec distribution
	#x axis is xsec, y axis is var bins


	if not choosebin:

	    xsec_hist=R.TH2F(var+'xsec','title',20,0.85,1.01,nbins,1,nbins+1)

	    #iter thru var bins
	    for resample in range(1,1+n_resamples):

	        #iter thru resamples
	        for bin in range(1,1+nbins):

		    #get xsec out of distribution and fill into xsec_hist
		    #divide by xsec_max to get relative height
		    xsec_hist.Fill( distribution.GetBinContent(resample,bin)/maxes[bin-1] , bin )

	    return xsec_hist


	elif choosebin:

	    #set limits according to max
	    low, hi = maxes[nbins-1]*0.7,maxes[nbins-1]*1.2


	    #new hist with xsec distribution
	    #x axis is xsec, y axis is var bins
	    xsec_hist=R.TH1F(var+'xsec','title',200,low,hi)

	    bin_center=self.get_As()[bgroup_index].reco2d.ProjectionY().GetBinCenter(nbins)

	    #iter thru var bins
	    for resample in range(1,1+n_resamples):


		#get xsec out of distribution and fill into xsec_hist
		#divide by xsec_max to get relative height
		xsec_hist.Fill( distribution.GetBinContent(resample, nbins)  )

		print distribution.GetBinContent(resample, nbins)/10**6



	    return xsec_hist, bin_center








    #################### Smoothing and plotting together ####################


    def smooth(self, plot=False):
	'''
	A function to compute the appropriate reweights, for smoothing all the raw data to the smooth curve {}

	This is done for each numerical bin of y and pT, across the bgroup 
	( where a bgroup consists of 8 binns of same numerical edges, each corresponding to a different coefficient)
	Within the bgroup a different reweight is determined for each bin of cost, and phi

	The args are the bgroups within self.main, and these should already be filled
	The arg plot, boolean, tells whether to make hist/smooth function comparison plots as well

	'''

	#if we want to plot, create pdf
	if plot:
	    c_dum=R.TCanvas('c_dum','c_dum',2000,800)
	    c_dum.SetFillColor(R.kGray)
	    c_dum.Print('smooth_data.pdf(') #leave open
	
	#init reweightor
	reweights=reweightor(self.ptedges,self.yedges,self.forin()[0].resamp) 

	#index to keep track of which bgroup we are on
	b_index= -1
	#iter through all the bgroups
	for ptlist in self.main:

	    #continue to get bgroups themselves, and iter thru
	    for bgroup in ptlist:

		b_index += 1

		#determine the angular coefficients for this bgroup
		A=[]
		for b in bgroup: #they are order b.i=0 .. b.i=7
		    A.append(b.calc()[0]) #the 0th element of the tuple is the ang coefficient

		#get the raw data histogram and its characteristics
		h_phi=bgroup[0].phics
		h_cost=bgroup[0].costcs
		h_both=bgroup[0].gen2d

		#get sum of weights
		sumw=bgroup[0].sumw()

		#get binning for phi and cost
		nbins_phi=bgroup[0].nbins_phi
		nbins_cost=bgroup[0].nbins_cost

		#create a hist to hold f vals
		h_f=R.TH2F('h_f'+str(bgroup[0]),'h_f; \phi ; cos(theta) ; cross section',nbins_phi,0,2*math.pi,nbins_cost,-1,1)


		#iter over phi vals
		for i in range(1,nbins_phi+1): #no of cost vals to correspond to n hist bins in x, hist binning starts at 1 so do likewise

		    #determine phi from the phi hist
    		    phi= h_phi.GetBinCenter(i)

		    #iter over cost vals
		    for j in range(1,nbins_cost+1): 
	
			#determine cost from the cost hist
			cost=h_cost.GetBinCenter(j)

		        #determine f from full {} definition
			#helpful vars
			sint=math.sqrt(1-cost**2)
			sin2t=2*sint*cost
    		        f= 1 + cost**2 + (1/2)*A[0]*(1-3*cost**2)+A[1]*sin2t*math.cos(phi)
		        f += (1/2)*A[2]*(sint**2)*math.cos(2*phi)+A[3]*sint*math.cos(phi)
			f += A[4]*cost +A[5]*(sint**2)*math.sin(2*phi)+A[6]*sin2t*math.sin(phi)
			f += A[7]*sint*math.sin(phi)

			#scale f according to the unpolarised cross section
			#sum of weights is total unpolarised cross section
			#divide total across bins
			f= f*(3/(16*math.pi))*sumw*(2*2*math.pi)/(nbins_phi*nbins_cost)

			#fill f val into its hist
			h_f.Fill(phi,cost,f)

		        #get h val from hist
		        h=h_both.GetBinContent(i,j)

			#determine reweight
			#avoid zero division / nearly empty bins
			if h<=10:
			    re=1
			    print 'zero division avoided','phi is',phi,'cost is',cost
			    print 'binn is',bgroup[0], 
			else:
			    re=f/h


			#give the reweight to the reweightor
			#first we need the appropriate binn
			#then the hist bin of the reweights histogram corresponding to phi, cost
			reweights.forin()[b_index].reweights.SetBinContent(i,j,re)

		#make unrolled histograms, if required
		if plot:

		    n_slices=7

		    #unroll hists (converting TH2F to TH1F
		    h_both_1d=unroll(h_both,n_slices,histid=str(bgroup[0])) #add str representation of binn
		    h_f_1d=unroll(h_f,n_slices,histid=str(bgroup[0]))	   #to avoid memory leaks

		    #format hists
		    h_both_1d.SetLineColor(R.kBlue)
		    h_f_1d.SetLineColor(R.kRed)

		    #create legend
		    leg_unroll=R.TLegend()
		    leg_unroll.AddEntry(h_both_1d, "Raw data")
		    leg_unroll.AddEntry(h_f_1d,"Smooth function")

		    c_unroll=R.TCanvas('c_unroll','c'+str(bgroup[0]),2000,800)
		    h_f_1d.Draw()
		    h_both_1d.Draw('esame')
		    leg_unroll.Draw()
		    c_unroll.Print('smooth_data.pdf')

		    
		
	#if we need to plot, finish pdf
	if plot:
	    c_dum.Print('smooth_data.pdf)') #close file

	return reweights

    def comp_bootstrapped(self,var,inclusive=False):
	'''
	Compare the raw data stored in self binns with the bootsrapped data in binn.resamples[]
	var, tells which var to compare for, phi, cost, pt, or eta
	'''


	#lists of all the hists
	raw_hists=[]
	raw_rms=[]

	
	#get hist types out of bins
	for i in range(len(self.get_As())):
	
	    #get the binns we want
	    b_raw=self.get_As()[i] #raw binn comes from holder

	    #choose hist based on var option
	    if var=='phi':
		raw_hists.append(b_raw.phics)
		raw_rms.append(b_raw.rms_reco(var))

	    if var=='cost':
		raw_hists.append(b_raw.costcs)
		raw_rms.append(b_raw.rms_reco(var))		

	    if var=='pt':
		raw_hists.append(b_raw.lep_pt)
		raw_rms.append(b_raw.rms_reco(var))

	    if var=='eta':
		raw_hists.append(b_raw.lep_eta)
		raw_rms.append(b_raw.rms_reco(var))

	
	#if we want inclusive hists instead
	if inclusive:
	    
	    #make inclusive hists
	    raw_inc=raw_hists[0]
	    raw_rms_inc=raw_rms[0]

	    #fill inclusive hists
	    for i in range(1,len(raw_hists) ):
	        raw_inc.Add(raw_hists[i])
	        raw_rms_inc.Add(raw_rms[i])

	    #convert to graphs
	    glist=hists_to_graph([raw_inc,raw_rms_inc])

	    #convert to multigraph
	    c1=R.TCanvas()
	    mg1=R.TMultiGraph()
	    for g in glist:
		mg1.Add(g)
	    mg1.Draw('ape')

	    #make the strings nicer
	    if var=='pt':
		nicevar='lepton p_{T}'
	    if var=='eta':
		nicevar='lepton \eta'
	    if var=='cost':
		nicevar='cos(theta)'
	    if var=='phi':
		nicevar='\phi'

	    #format the multigraph
	    mg1.SetTitle('Different sources of error on '+nicevar+', inclusive over y and p_{T};'+nicevar+'; cross section')
	    
	    #make the legend
	    leg1=R.TLegend(0.7,0.7,0.9,0.9)
	    leg1.AddEntry(glist[0],'Raw '+nicevar+' data' )
	    leg1.AddEntry(glist[1],'Bootstrapped '+nicevar+' data' )
	    leg1.Draw()

	    c1.Print(var+'_boot_inclusive.pdf')

	    return


	#otherwise

	#convert from hists to graphs
	graphs=[]
	for i in range(len(raw_hists)):
	    graphs.append(hists_to_graph([raw_hists[i],raw_rms[i] ]))

	#graph everything and print
	filename=var+'_comp_boot.pdf'
	c_dum=R.TCanvas()
	c_dum.SetFillColor(R.kGray)
	c_dum.Print(filename+'(')

	#get the list of graphs we made and draw each
	for i in range(len(graphs)): #each is a list of 2 graphs
	
	    #make and format canvas
	    c1=R.TCanvas()

	    #make multigraph and draw
	    mg1=R.TMultiGraph()
	    for g in graphs[i]:
		mg1.Add(g)
	    mg1.Draw('ape')

	    #make the strings nicer
	    if var=='pt':
		nicevar='lepton p_{T}'
	    if var=='eta':
		nicevar='lepton \eta'
	    if var=='cost':
		nicevar='cos(theta)'
	    if var=='phi':
		nicevar='\phi'

	    #format the multigraph
	    mg1.SetTitle('Different sources of error on '+nicevar+' '+self.get_As()[i].nice_str()+';'+nicevar+'; cross section')
	    
	    #make the legend
	    leg1=R.TLegend(0.7,0.7,0.9,0.9)
	    leg1.AddEntry(graphs[i][0],'Raw '+nicevar+' data' )
	    leg1.AddEntry(graphs[i][1],'Bootstrapped '+nicevar+' data' )
	    leg1.Draw()

	    c1.Print(filename)

	#close file
	c_dum.Print(filename+')')

    def comp_smoothed(self,files,var,hists=True):
	'''
	Compare the raw data stored in self binns with the smoothed data gotten from ttree
	ttree, ttree obj, has events we can smooth
	var, tells which var to compare for, phi, cost, pt, or eta
	hists, bool, tells whether to return TH1s instead of writing it to pdf
	'''

	#lists of all the hists
	raw_hists=[]
	smooth_hists=[]

	#get reweights
	reweights=self.smooth()

	#get trimmed tree name
	treename=str(self.forin()[0])

	#smooth trimmed tree, once for each part
	for filename in files:
	    reweights.make_smoothed(filename,treename)

	#get hists out of bins
	for i in range(len(self.get_As())):
	
	    #get the binns we want
	    b_raw=self.get_As()[i]
	    b_smooth=reweights.forin()[i]

	    #choose based on var option
	    if var=='phi':
		raw_hists.append(b_raw.phics)
		smooth_hists.append(b_smooth.phics)

	    if var=='cost':
		raw_hists.append(b_raw.costcs)
		smooth_hists.append(b_smooth.costcs)


	    if var=='pt':
		print b_raw
		print b_raw.lep_pt
		raw_hists.append(b_raw.lep_pt)
		smooth_hists.append(b_smooth.lep_pt)

	    if var=='eta':
		raw_hists.append(b_raw.lep_eta)
		smooth_hists.append(b_smooth.lep_eta)

	print raw_hists
	print smooth_hists

	#return hists here if asked
	if hists:
	    return raw_hists[0], smooth_hists[0]

	#convert from hists to graphs
	graphs=[]
	for i in range(len(raw_hists)):
	    graphs.append(hists_to_graph([raw_hists[i],smooth_hists[i] ]))

	#graph everything and print
	filename=var+'_comp_smoothed.pdf'
	c_dum=R.TCanvas()
	c_dum.SetFillColor(R.kGray)
	c_dum.Print(filename+'(')

	#make the strings nicer
	if var=='pt':
	    nicevar='lepton p_{T}'
	if var=='eta':
	    nicevar='lepton \eta'
	if var=='cost':
	    nicevar='cos(theta)'
	if var=='phi':
	    nicevar='\phi'


	#get the list of graphs we made and draw each
	for i in range(len(graphs)): #each is a list of graphs
	
	    #make and draw mg
	    c1=R.TCanvas()
	    mg1=R.TMultiGraph()
	    mg1.SetTitle(nicevar+' data before and after smoothing, '+self.get_As()[i].nice_str()+';'+nicevar+'; cross section')
	    mg1.Add(graphs[i][0])
	    mg1.Add(graphs[i][1])
	    mg1.Draw('ape')

	    #make the legend
	    leg1=R.TLegend(0.7,0.7,0.9,0.9)
	    leg1.AddEntry(graphs[i][0],'Raw '+nicevar+' data' )
	    leg1.AddEntry(graphs[i][1],'Smooth '+nicevar+' data' )
	    leg1.Draw()

	    c1.Print(filename)

	#close file
	c_dum.Print(filename+')')


    def comp_all(self, reweights, var, hist=False):
	'''
	Compare:
	raw data with w^2 error
	bootstrapped data with rms error
	smoothed raw data with w^2 error
	smoothed bootstrapped data with rms error
	Args:
	ttreelist, list of ttree obj, each has events we can smooth, first is real data, rest is resampled
	var, tells which var to compare for, phi, cost, pt, or eta
	hist, if true tells to return hists to comp_inclusive before making graphs
	
	'''

	#lists of all the hists
	raw_hists=[]
	smooth_hists=[]
	rms_hists=[]
	boot_hists=[]

	
	#get first 3 hist types out of bins
	#for i in range(len(self.get_As())):
	if True:
	
	    #get the binns we want
	    b_raw=self.get_As()[0] #raw binn comes from holder
	    b_smooth=reweights.forin()[0] #smoothed binn comes from reweightor

	    print 'b raw',b_raw
	    print 'b smooth', b_smooth


	    #choose hist based on var option
	    if var=='phi':
		raw_hists.append(b_raw.phics)
		smooth_hists.append(b_smooth.phics)
		rms_hists.append(b_smooth.rms_reco(var))
		boot_hists.append(b_raw.rms_reco(var))

	    elif var=='cost':
		raw_hists.append(b_raw.costcs)		
		smooth_hists.append(b_smooth.costcs)
		rms_hists.append(b_smooth.rms_reco(var))
		boot_hists.append(b_raw.rms_reco(var))

	    elif var=='pt':
		raw_hists.append(b_raw.lep_pt)
		smooth_hists.append(b_smooth.lep_pt)
		rms_hists.append(b_smooth.rms_reco(var))
		boot_hists.append(b_raw.rms_reco(var))

	    elif var=='eta':
		raw_hists.append(b_raw.lep_eta)
		smooth_hists.append(b_smooth.lep_eta)
		rms_hists.append(b_smooth.rms_reco(var))
		boot_hists.append(b_raw.rms_reco(var))

	    elif var=='gen2d':
		raw_hists.append(b_raw.gen2d)
		smooth_hists.append(b_smooth.gen2d)
		rms_hists.append(b_smooth.rms_reco(var))
		boot_hists.append(b_raw.rms_reco(var))

	    elif var=='reco2d':
		raw_hists.append(b_raw.reco2d)
		smooth_hists.append(b_smooth.reco2d)
		rms_hists.append(b_smooth.rms_reco(var))
		boot_hists.append(b_raw.rms_reco(var))

	    elif var=='bos2d':
		raw_hists.append(b_raw.bos2d)
		smooth_hists.append(b_smooth.bos2d)
		rms_hists.append(b_smooth.rms_reco(var))
		boot_hists.append(b_raw.rms_reco(var))

	    elif var=='dress2d':
		raw_hists.append(b_raw.dress2d)
		smooth_hists.append(b_smooth.dress2d)
		rms_hists.append(b_smooth.rms_reco(var))
		boot_hists.append(b_raw.rms_reco(var))


	#we might just want the hists (so that we can add them inclusively before conversion)
	if hist:
	    return raw_hists, smooth_hists, rms_hists, boot_hists

	#convert from hists to graphs
	graphs=[]
	for i in range(len(raw_hists)):
	    graphs.append(hists_to_graph([raw_hists[i],smooth_hists[i],raw_rms[i],smooth_rms[i] ]))

	#graph everything and print
	filename=var+'_comp_all.pdf'
	c_dum=R.TCanvas()
	c_dum.SetFillColor(R.kGray)
	c_dum.Print(filename+'(')

	#get the list of graphs we made and draw each
	for i in range(len(graphs)): #each is a list of 4 graphs
	
	    #make and format canvas
	    c1=R.TCanvas()
	    #pt plots should be logarithmic

	    #make multigraph and draw
	    mg1=R.TMultiGraph()
	    for g in graphs[i]:
		mg1.Add(g)
	    mg1.Draw('ape')

	    #make the strings nicer
	    if var=='pt':
	        nicevar='lepton p_{T}'
	    if var=='eta':
	        nicevar='lepton \eta'
	    if var=='cost':
	        nicevar='cos(theta)'
	    if var=='phi':
	        nicevar='\phi'

	    #format the multigraph
	    mg1.SetTitle('Different sources of error on '+nicevar+' '+self.get_As()[i].nice_str()+';'+nicevar+'; cross section')
	    
	    #make the legend
	    leg1=R.TLegend(0.7,0.7,0.9,0.9)
	    leg1.AddEntry(graphs[i][0],'Raw '+nicevar+' data' )
	    leg1.AddEntry(graphs[i][1],'Smooth '+nicevar+' data' )
	    leg1.AddEntry(graphs[i][2],'Bootstrapped '+nicevar+' data' )
	    leg1.AddEntry(graphs[i][3],'Smoothed bootstrapped '+nicevar+' data' )
	    leg1.Draw()

	    c1.Print(filename)

	#close file
	c_dum.Print(filename+')')


    def comp_inclusive(self,ttree,var):
	'''
	Obtains the results of comp_all above, which is a bunch of hists for different pt, y binns. 
	Combines them into inclusive hists then converts to graphs
	'''
	
	filename=var+'_comp_inclusive.pdf'

	#get lists of histograms
	raw,smooth,raw_rms,smooth_rms=self.comp_all(ttree,var,hist=True)

	#make inclusive hists
	raw_inc=raw[0]
	smooth_inc=smooth[0]
	raw_rms_inc=raw_rms[0]
	smooth_rms_inc=smooth_rms[0]

	#fill inclusive hists
	for i in range(1,len(raw) ):
	    raw_inc.Add(raw[i])
	    smooth_inc.Add(smooth[i])
	    raw_rms_inc.Add(raw_rms[i])
	    smooth_rms_inc.Add(smooth_rms[i])

	#convert to graphs
	glist=hists_to_graph([raw_inc,smooth_inc,raw_rms_inc,smooth_rms_inc])

	#make into multigraph and draw
	c1=R.TCanvas()
	mg1=R.TMultiGraph()
	for g in glist:
	    mg1.Add(g)
	mg1.Draw('ape')

	#make the strings nicer
	if var=='pt':
	    nicevar='lepton p_{T}'
	if var=='eta':
	    nicevar='lepton \eta'
	if var=='cost':
	    nicevar='cos(theta)'
	if var=='phi':
	    nicevar='\phi'

	#format the multigraph
	mg1.SetTitle(nicevar+' error formulations, inclusive over y and p_{T};'+nicevar+'; cross section')
	    
	#make the legend
	leg1=R.TLegend(0.7,0.7,0.9,0.9)
	leg1.AddEntry(glist[0],'Raw '+nicevar+' data' )
	leg1.AddEntry(glist[1],'Smooth '+nicevar+' data' )
	leg1.AddEntry(glist[2],'Bootstrapped '+nicevar+' data' )
	leg1.AddEntry(glist[3],'Smoothed bootstrapped '+nicevar+' data' )
	leg1.Draw()

	#print canvas
	c1.Print(filename)

		
		
	
class reweightor(object):
    '''
    Object for utilizing the smoothing reweights found above

    First, is filled with reweights 
    Then, iterates through tree events, determining the corresponding reweight and filling
    '''

    def __init__(self,ptedges,yedges,resamp):
	'''
	Args
	ptedges, y edges, lists that define binning
	nbins_phi,nbins_cost, define minibinning
	resample, number of resample

	'''
	
	#divide space into binns like in holder

	#init list
    	self.main=[]

	#init attributes
	self.ptedges=ptedges
	self.yedges=yedges
	self.resamp=resamp

	#go through y edges and add to main list
	for j in range(len(yedges)-1):
	
	    ptlist=[]

    	    #go thru pt edges and add to pt list
    	    for k in range(len(ptedges) - 1):

	
		#make binn an A8 to distinguish it from bholder binns
	    	ptlist.append(binn(ptedges[k],ptedges[k+1],yedges[j],yedges[j+1],'prefsrw',8,self.resamp))

	    #add to main list
	    self.main.append(ptlist)



    def forin(self):
	'''
	returns contents of the self.main 3d list as a simple list for traversing with a for loop
	'''
	forlist=[]
	
	#iter thru all contents
	for ptlist in self.main:
	    for b in ptlist:
		forlist.append(b)

	return forlist


    def avg(self):

	#iter over bins
	for b in self.forin():
	    
	    #init count vars
	    mysum=0
	    count=0

	    #loop thru reweights
	    for i in range(1,b.nbins_phi+1):
		for j in range(1,b.nbins_cost+1):
		    mysum += b.reweights.GetBinContent(i,j)
		    count += 1

	#average and print
	return mysum/count


    def get_reweight(self,event):
	'''
	Determines the appropriate reweight that should modify the event's weight
	e is a ttree event
	returns a float, reweight
	'''

	#determine which binn the event falls in
	for b in self.forin():

	    if b.contains_trimmed(event):

		#find bin coordinates for reweight
		i=int( math.ceil(event.phi/(2*math.pi/20) ) ) #20 is nbins_phi, 2pi is width of hist
		j=int( math.ceil((event.cost+1)/(2/20) ) ) #20 is nbins_cost, 2 is width of hist

		#print 'phi is ', event.phi, ' and should be > ',b.reweights.GetXaxis().GetBinLowEdge(i),' and < ',b.reweights.GetXaxis().GetBinLowEdge(i+1)
		if not ((event.phi > b.reweights.GetXaxis().GetBinLowEdge(i)) and (event.phi < b.reweights.GetXaxis().GetBinLowEdge(i+1))):
		    print 'False for phi'


	
		if not ((event.cost > b.reweights.GetYaxis().GetBinLowEdge(j)) and (event.cost < b.reweights.GetYaxis().GetBinLowEdge(j+1))):

		    print 'cost is ', event.cost, ' and should be > ',b.reweights.GetYaxis().GetBinLowEdge(j),' and < ',b.reweights.GetYaxis().GetBinLowEdge(j+1)
		    print 'False for cost'

		#get reweight from TH2
		return b.reweights.GetBinContent(i,j)



    def make_smoothed(self,filename,treename):
	'''
	Makes params in ttree smooth with reweights applied
	For efficiency does only one event at a time (thus only 1 ttree loop required)
	Args:
	filename, tells where the ttree object (should be a main tree) containing required data is
	var, tells which variable hist to fill
	'''

	file1=R.TFile(filename)
	tree1=file1.Get(treename)

	print 'start ttree loop'
	for event in tree1:

	    #find the appropriate binn(s)
	    for b in self.forin():

		if b.contains_trimmed(event):

	    	    #get the reweight
	    	    reweight=self.get_reweight(event)

	    	    
		    #for pt in event.lep_pt:
			#b.lep_pt.Fill(pt,event.w*reweight)

		    #for eta in event.lep_eta:
			#b.lep_eta.Fill(eta,event.w*reweight)

		    b.phics.Fill(event.phi,event.w*reweight)

		    b.costcs.Fill(event.cost,event.w*reweight)

		    b.gen2d.Fill(event.phi, event.cost,event.w*reweight)

		    lep_eta=tuple(event.lep_eta)
		    lep_pt=tuple(event.lep_pt)

		    for i in range( len( lep_eta) ):

	    		b.reco2d.Fill(lep_eta[i],lep_pt[i],event.w*reweight )

		    dress_eta=tuple(event.dress_eta)
		    dress_pt=tuple(event.dress_pt)

		    for i in range( len( dress_eta) ):
			b.dress2d.Fill(dress_eta[i],dress_pt[i],event.w*reweight )

		    b.bos2d.Fill(event.y, event.pt,event.w*reweight)


	return None



    def get_distribution(self,var,bgroup_index,n_resamples,nbins ):
	'''
	Looks at the distribution of resampled cross section w.r.t a specified parameter.
	Ideally the distribution should be Gaussian
	Args: var, tells parameter to choose bins for
	'''

	#pick the binn we want
	mybinn = self.forin()[bgroup_index]

	#get list of 1d histograms w.r.t var
	hlist=[]
	for b in mybinn.resamples: #iter thru resamples

	    #grab hist according to var
	    #also define other var-based variables here
	    if var=='phi':
		hlist.append( b.gen2d.ProjectionX(str(b),1,20,'e') )
		low,hi = 0,2*math.pi
		binsize=(hi-low)/20
	    elif var=='cost':
		hlist.append( b.gen2d.ProjectionY(str(b),1,20,'e') )
		low,hi = -1,1
		binsize=(hi-low)/20
	    elif var=='eta':
		hlist.append( b.reco2d.ProjectionX(str(b),1,20,'e') )
		low,hi = -2.4,2.4
		binsize=(hi-low)/20
	    elif var=='pt':
		hlist.append( b.reco2d.ProjectionY(str(b),1,20,'e') )
		low,hi = 0,60
		binsize=(hi-low)/20

	#make a TH2 of the resample distribution
	# x axis is resample #, y axis is var bin
	distribution=R.TH2F(var,var,n_resamples,1,n_resamples+1,nbins,1,nbins+1)

	#format
	distribution.GetXaxis().SetTitle(var)
	distribution.GetYaxis().SetTitle('resample number')
	distribution.GetYaxis().SetNdivisions(n_resamples)

	#iter thru var bins
	for resample in range(1,1+n_resamples):

	    #iter thru resamples
	    for bin in range(1,1+nbins):
	    
		#now fill the distribution th2
		distribution.SetBinContent( resample,bin, hlist[resample-1].GetBinContent(bin) )

	return distribution


    def plot_distribution(self,var,bgroup_index,n_resamples,nbins=10 , choosebin=False):
	'''
	After the resample distribution for a given var is stored in a TH2F by get_distribution(),
	makes an unrolled plot of the xsec of the distribution unrolled along the var distribution
	'''

	#get the distribution
	distribution = self.get_distribution(var,bgroup_index,n_resamples,nbins)

	#c1=R.TCanvas()
	#c1.SetGrid()
	#distribution.Draw()
	#c1.Print('gauss.pdf')

	#find the maximum xsec for each binn so that we can plot in terms of relative xsec
	maxes=[]
	for bin in range(1,1+nbins): #iter thru bins
	    xsec_vals=[]
	    for resample in range(1,1+n_resamples): #iter thru resamples
		xsec_vals.append( distribution.GetBinContent(resample,bin) )

	    #find the max of the resample xsecs in this bin
	    maxes.append( max( xsec_vals ) )
	

	#new hist with xsec distribution
	#x axis is xsec, y axis is var bins


	if not choosebin:

	    xsec_hist=R.TH2F(var+'xsecsmooth','title',20,0.9,1,nbins,1,nbins+1)

	    #iter thru var bins
	    for resample in range(1,1+n_resamples):

	        #iter thru resamples
	        for bin in range(1,1+nbins):

		    #get xsec out of distribution and fill into xsec_hist
		    #divide by xsec_max to get relative height
		    xsec_hist.Fill( distribution.GetBinContent(resample,bin)/maxes[bin-1] , bin )

	    return xsec_hist

	elif choosebin:

	    #set limits according to max
	    low, hi = maxes[nbins-1]*0.7,maxes[nbins-1]*1.2

	    xsec_hist=R.TH1F(var+'xsecsmooth','title',200,low,hi)

	    bin_center=self.forin()[bgroup_index].reco2d.ProjectionY().GetBinCenter(nbins)

	    #iter thru var bins
	    for resample in range(1,1+n_resamples):


		#get xsec out of distribution and fill into xsec_hist
		#divide by xsec_max to get relative height
		xsec_hist.Fill( distribution.GetBinContent(resample,nbins)  )

		print distribution.GetBinContent(resample, nbins)/10**6


	    return xsec_hist, bin_center
	    
    def make_inclusive(self):
	'''
	Takes the data spread out between binns in the A0 hists and combines, as well as combines all the resamples across bins
	'''

	#combine all the A8 hists
	#they all go in the A8 hists of the first bin
	for i in range(1,len(self.forin() ) ): #iters thru binns

	    #add gen2d param hists
	    self.forin()[0].costcs.Add( self.forin()[i].costcs )
	    self.forin()[0].phics.Add( self.forin()[i].phics )
	    self.forin()[0].gen2d.Add( self.forin()[i].gen2d )

	    #add reco2d param hists
	    self.forin()[0].reco2d.Add( self.forin()[i].reco2d )

	    #add bos2d param hists
	    self.forin()[0].bos2d.Add( self.forin()[i].bos2d )

	    #add dress2d param hists
	    self.forin()[0].dress2d.Add( self.forin()[i].dress2d )


	#now do the resamples
	#we want to put all the resamples across parts into the resamples list of the first pt, y binn
	#again iter thru all pT, y binns
	for k in range(1,len(self.forin() ) ):

	    #grab the resamples list of the kth binn
	    for resample in self.forin()[k].resamples :

		#place each resampled binn in the 0th binn list
		self.forin()[0].resamples.append( resample ) 

	return 






    def convert_to_holder(self):
	newholder=holder(self.ptedges,self.yedges)
	newholder.main=self.main

	return newholder

    def write_hists(self, filename):

	for b in self.forin():
	    b.write(filename)

    def plot_hists(self):
	
	histlist=[]
	for b in self.forin():
	    histlist.append(b.phics)
	    histlist.append(b.costcs)
	    histlist.append(b.lep_pt)
	    histlist.append(b.lep_eta)

	c_dum=R.TCanvas()
	c_dum.SetFillColor(R.kGray)
	filename='reweightor_hists.pdf'
	c_dum.Print(filename+'(')

	for h in histlist:
	    c1=R.TCanvas()
	    h.Draw()
	    c1.Print(filename)

	c_dum.Print(filename+')')

	

class Agraph(object):
    '''
    Essentially a tgrapherrors object specifically for plotting angular coefficients
    '''
    def __init__(self,xvals,yvals,xerrs,yerrs,dtype, etype,title):
	'''
	Attributes:
	graph: tgrapherrors obj
	xvals, yvals: arrays of the graph points
	xerrs: error to include on x, usually 0
	yerrs: error on ang coefficient values
	dtype: string, whether genw or prefsrw
	etype: string, whether norm or rms 
	'''
	#create graph
	self.graph=R.TGraphErrors(len(xvals),xvals,yvals,xerrs,yerrs)

	#add array attributes
	self.xvals=xvals
	self.yvals=yvals
	self.xerrs=xerrs
	self.yerrs=yerrs

	#format graph

	#type specific formatting
	if dtype == 'genw':
	    if etype == 'norm':
		self.graph.SetLineColor(R.kBlack) 
		self.graph.SetMarkerStyle(21)
	    elif etype == 'rms':
		self.graph.SetLineColor(R.kGreen) 
		self.graph.SetMarkerStyle(21)
	elif dtype == 'prefsrw':
	    if etype == 'norm':
		self.graph.SetLineColor(R.kRed) 
	        self.graph.SetFillColor(R.kWhite)
		self.graph.SetMarkerStyle(20)
	    elif etype == 'rms':
		self.graph.SetLineColor(R.kBlue) 
		self.graph.SetMarkerStyle(20)
		self.graph.SetFillColorAlpha(R.kBlue,0.5)

	#general formatting
        self.graph.SetTitle(title)
        #self.graph.GetXaxis().SetMaxDigits(3)
        #self.graph.GetYaxis().SetMaxDigits(3)
        self.graph.SetLineWidth(1)
	self.graph.GetXaxis().SetLimits(xvals[0]-xerrs[0],xvals[-1]+xerrs[-1])

    def get_rels(self):
	'''
	returns the yerrs, converted to perentage relative errors calculated against yvals
	'''

	rels=array.array('f')

	#iter thru yvals
	for i in range(len(self.yvals)):
	    rels.append(abs(100*self.yerrs[i]/self.yvals[i]))

	return rels



#####################################################################

####################### functions below here ########################	

#####################################################################


def makehists(name):
    '''
    creates and outputs four TH1Fs of the kind we need for storing data
    Args:
	name, string which forms first part of hist name (rest is generic endings)
    Returns, tupe of four hists		
    '''

    #create hists

    num=R.TH1F(name+'_num',name+'_num',200,-1.5,1.5)

    den=R.TH1F(name+'_den',name+'_den',2,-2,2) #will be filled with +/- 1

    num_err=R.TH1F(name+'_num_err',name+'_num_err',200,0,2.5)

    den_err=R.TH1F(name+'_den_err',name+'_den_err',1,0,2) #will be filled with 1

    return num, den, num_err, den_err


def clone_group(bgroup):
    '''
    After an A0 binn is filled and resampled, passes its resamples to the other binns in the bgroup
    '''
    b0=bgroup[0]

    for i in range(1,8):

	#pass resamples and rename
	for resamp in b0.resamples:

	    #copy resample evs to new binn
	    new=binn(b0.ptmin,b0.ptmax,b0.ymin,b0.ymax,b0.type,i,resamp.resamp,part=b0.part)

	    #pass evs
	    for e in resamp.evs:
		new.evs.append(e)

	    #add to list
	    bgroup[i].resamples.append(new)

    return None




def plot(blist, i, title):
    """
    Function that graphs values of An over bins of pT for a given y
    Args: An: the nth angular coefficient (actually a function)
        blist: list of pt bins that we are plotting A_i for, each element is a list of binns for all of
    	the A_i's
	i: tells which coefficient
        title: string that will be title of graph created
    Returns: tuple of TGraphErrors obj (and array of errs if error=True)
    """
    
    #init graph value arrays
    xvals=array.array('f')
    yvals=array.array('f')
    xerrs=array.array('f')
    yerrs=array.array('f')
    
    
    #find array vals for graph from binns
    for ilist in blist: #blist is a list of lists, ilist gives a list of binns for each Ai

	#correct bin is the ith one in the list
	b=ilist[i]

        xvals.append(0.5*(b.ptmax+b.ptmin))
	xerrs.append(0)

	#use calc method to get A and error
	A, err =b.calc()

        yvals.append(A)
	yerrs.append(err)
    
    #return Agraph object (see treesuite)
    return Agraph(xvals,yvals,xerrs,yerrs,blist[0][0].type,'norm',title)	



def plot_resampled(blist, i, title):
    '''
    Same goal as plot, but with errors computed by resampling, and taking rms error
	Args:
        blist: list of pt bins that we are plotting A_i for, each element is a list of binns for all of
    	the A_i's
	i: tells which coefficient
        title: string that will be title of graph created
    Returns: tuple of TGraphErrors obj (and array of errs if error=True)
    '''

    #init graph value arrays
    xvals=array.array('f')
    yvals=array.array('f')
    xerrs=array.array('f')
    yerrs=array.array('f')

    #find array vals for graph from binns
    for ilist in blist: #blist is a list of lists, ilist gives a list of binns for each Ai

	#correct bin is the ith one in the list
	b=ilist[i] 

        xvals.append(0.5*(b.ptmax+b.ptmin))
	xerrs.append(0.5*(b.ptmax-b.ptmin)) #here we want xerr to be pt width of bin

	#fill with rms error
	Ahat, err = b.rms()
        yvals.append(Ahat)
        yerrs.append(err)
    
    #returns Agraph object (see treesuite)
    return Agraph(xvals,yvals,xerrs,yerrs,blist[0][0].type,'rms',title)


def unroll(h2d,n_slices,histid=str(0)):
    '''
    Function to unroll a 2d histogram into a 1d one
    This PyRoot version translated from the cpp form in TGHistogramUtils::unrollHistogram

    Args
    h2d, a TH2F to be unrolled
    n_slices,int,tels how many y slices to pick out and display (instead of them all)
    histid, string for naming the new hist and avoiding memory leaks

    Returns, h2d, a TH2F
    '''

    #get all the information about this 2d hist
    title=h2d.GetTitle()
    nx=h2d.GetXaxis().GetNbins()
    ny=h2d.GetYaxis().GetNbins()
    label=h2d.GetYaxis().GetTitle()

    #make a new 1d hist
    h1d=R.TH1F(title+'1d'+histid,title+'1d; For\,each\,cos(theta),\,0 < \phi < 2\pi; cross section',nx*n_slices,0,1) #can go from 0 to 1 since we are relabeling later anyway

    #pick out certain slices of y to display more sanely
    ylist=[]
    for i in range(n_slices):
	ylist.append(int((i+0.5)*ny/n_slices))

    #fill and label new bins
    bincount=1
    binsize=h2d.GetYaxis().GetBinCenter(2)-h2d.GetYaxis().GetBinCenter(1)

    #iter over y
    for j in ylist:

	binlabel=str(h2d.GetYaxis().GetBinCenter(j-1)+binsize/2)+'<'+label+'<'+str(h2d.GetYaxis().GetBinCenter(j+1)-binsize/2)
	
	#iter over x
	for i in range(1,nx):

	    #fill the 1d bin with contents of corresponding 2d binn
	    h1d.SetBinContent(bincount,h2d.GetBinContent(i,j))

	    #label only the first bin in each j series
	    if i==1:
		h1d.GetXaxis().SetBinLabel(bincount,binlabel)

	    bincount += 1

    #format hist
    h1d.GetXaxis().LabelsOption('h')

    return h1d



def hists_to_graph(histlist):
    '''
    Given a list of hists with errors, converts them to tgrapherrors and plots them offset 
    so that all the error bars can be easily seen
    '''

    #colors for plotting
    colors=[R.kBlue+1,R.kRed+1,R.kViolet,R.kBlack,R.kRed-2,R.kBlue-2,R.kGray,R.kGreen-3,R.kGreen+4,R.kOrange]

    #get data out of first hist
    h0=histlist[0]
    nbins=h0.GetXaxis().GetNbins()
    binstep=h0.GetBinCenter(2)-h0.GetBinCenter(1)
    offset=binstep/(len(histlist)+2) #leave some space at edges

    #check to ensure the hists align properly
    for h in histlist:
	if (h.GetBinCenter(1)-h0.GetBinCenter(1)) > 10**(-6):
	    print 'histograms not aligned'
	if (h.GetBinCenter(2)-h.GetBinCenter(1))-(binstep) > 10**(-6):
	    print 'histograms not aligned'
	if h.GetXaxis().GetNbins() != nbins:
	    print 'histograms not aligned'

    #create holders arrays
    xvalarrays=[]
    yvalarrays=[]
    xerrarrays=[]
    yerrarrays=[]

    #fill in the arrays for each hist
    for h in histlist:
	xvalarrays.append(array.array('f'))
	yvalarrays.append(array.array('f'))
	xerrarrays.append(array.array('f'))
	yerrarrays.append(array.array('f'))

    #now iterate thru the hists and get the vals out
    for i in range(len(histlist)):

	#iter thru hist bins
	for j in range(1,nbins+1):

	    #yval and yerr are simple
            yvalarrays[i].append(histlist[i].GetBinContent(j))
	    yerrarrays[i].append(histlist[i].GetBinError(j))

	    #xvals have to be offset
	    xval=histlist[i].GetBinCenter(j) #start at center
	    xoffset=binstep/(2*len(histlist) )
	    xvalarrays[i].append( xval + (i-1)*xoffset ) #second one will be in center

	    #fill in xerrs
	    if i==0: #the first has to show the bin width
		xerrarrays[i].append(binstep/2)
	    elif i==2: #the botstrapped data is an opaque rectangle
		xerrarrays[i].append(binstep/2)
	    else:
		xerrarrays[i].append(0)

    #make into tgrapherrors
    tgraphs=[]
    for i in range(len(histlist)):
	tgraphs.append(R.TGraphErrors(len(xvalarrays[i]),xvalarrays[i],yvalarrays[i],xerrarrays[i],yerrarrays[i]))

	#format tgraph
	tgraphs[i].SetTitle(histlist[i].GetTitle())
	tgraphs[i].SetName(histlist[i].GetName()+'1') #to avoid memory leak
	tgraphs[i].SetLineColor(colors[i])

	if i ==2: #this will be bootstrapped data
	    tgraphs[i].SetLineColorAlpha(colors[i],0.4)

    return tgraphs
	


	

def make_binns(ttree,n_pt,n_y):
    '''
    Given a ttree of events, break up into bins of pt, y with roughly equal number of events per bin
    n_pt, n_y, ints, tell how many binns along each
    '''

    #get representative sample
    ptlist=[]
    ylist=[]
    for event in ttree:
	
	#convert to ev
	ptlist.append(event.prefsrw_pt)
	ylist.append(abs(event.prefsrw_y))

    #sort the lists
    ptlist.sort()
    ylist.sort()

    #find pt edges
    ptedges=[0]
    thresh=len(ptlist)/n_pt

    for i in range(n_pt):
	#find the cutoff index
	j=0
	while j+1 < (i+1)*thresh:
	    j += 1

	#while loop exits as thresh is reached
	if j > len(ptlist): #end of the list
	    ptedges.append( math.ceil( ptlist[-1] ) ) #just append the last one
	else:
	    ptedges.append( math.ceil( ptlist[j] ) )

    #same for y
    yedges=[0]
    thresh=len(ylist)/n_y

    for i in range(n_y):
	#find the cutoff index
	j=0
	while j+1 < (i+1)*thresh:
	    j += 1

	#while loop exits as thresh is reached
	if j > len(ylist): #end of the list
	    yedges.append( math.ceil( ylist[-1]*100 )/100 ) #just append the last one
	else:
	    yedges.append( math.ceil( ylist[j]*100 ) /100 ) #factors of 100 round to 100th place



    return ptedges, yedges


def change_edges(ptedges,yedges,index):
    '''
    Takes a full list og ptedges, yedges and returns a specific bin of each based on the index of the bgroup
    used for parallelizing smooth_par.py
    '''

    if index >= ( len(ptedges) -1)*( len(yedges) -1):
	index = index - ( len(ptedges) -1)*( len(yedges) -1)

    #divide index into pt and y components
    y_comp=int(math.floor(index / (len(ptedges) - 1) ) )
    pt_comp= int( index - y_comp*(len(ptedges) - 1) )

    #grab edges from the original list
    new_y=[ yedges[y_comp], yedges[y_comp+1] ]
    new_pt=[ ptedges[pt_comp], ptedges[pt_comp+1] ]

    return new_pt, new_y
    


    





