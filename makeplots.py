
from __future__ import division
import ROOT as R
import math
import array
import random
import time
import sys
import treesuite as t

'''
Module for making plots from the th1 .root file output of smooth_par
'''

################################    Functions     ################################ 


def get(var,index,histlist,ptedges,yedges):
    '''
    Gets a hist as specified by var, index and type from the file as specified by index. 
    To keep everything in scope, also adds the hist to the list
    '''

    #get file name
    #filestring='/afs/cern.ch/user/c/cbunker/Christian/CMSSW_10_3_0/ang_job/rootfiles/'
    filestring='/eos/user/c/cbunker/'
    filestring += var
    file0=R.TFile(filestring+str(index)+'.root')

    #we need the edges corresponding to file number
    new_pt, new_y = t.change_edges(ptedges, yedges, index) #outputs edges as list

    #convert edges to string representation of the corresponding binn
    raw_string='A0_'+str(new_pt[0])+'_'+str(new_pt[1])+'_'+str(new_y[0])+'_'+str(new_y[1])
    smooth_string='A8_'+str(new_pt[0])+'_'+str(new_pt[1])+'_'+str(new_y[0])+'_'+str(new_y[1])

    #update according to var
    if var=='eta':
	raw_string += 'lep_eta'
	smooth_string += 'lep_eta'
    elif var=='pt':
	raw_string += 'lep_pt'
	smooth_string += 'lep_pt'
    elif var=='phi':
	raw_string += 'phics'
	smooth_string += 'phics'
    elif var=='cost':
	raw_string += 'costcs'
	smooth_string += 'costcs'
    elif var=='gen2d':
        raw_string += 'gen2d'
        smooth_string += 'gen2d'     
    elif var=='reco2d':
        raw_string += 'reco2d'
        smooth_string += 'reco2d'

    #get the hists
    raw, smooth = file0.Get(raw_string) , file0.Get(smooth_string)
    #histlist.append(raw)
    #histlist.append(smooth)
    print raw, smooth

    return raw, smooth





def get4(var,index,ptedges,yedges):
    '''
    Gets a hist as specified by var, index and type from the file as specified by index. 
    To keep everything in scope, also adds the hist to the list
    '''

    #get file name
    #filestring='/afs/cern.ch/user/c/cbunker/Christian/CMSSW_10_3_0/ang_job/rootfiles/'
    filestring='/eos/user/c/cbunker/'
    filestring += var
    file0=R.TFile(filestring+str(index)+'.root')

    #we need the edges corresponding to file number
    new_pt, new_y = t.change_edges(ptedges, yedges, index) #outputs edges as list

    #convert edges to string representation of the corresponding binn
    raw_string='A0_'+str(ptedges[0])+'_'+str(ptedges[1])+'_'+str(yedges[0])+'_'+str(yedges[1])
    smooth_string='A8_'+str(ptedges[0])+'_'+str(ptedges[1])+'_'+str(yedges[0])+'_'+str(yedges[1])
    boot_string='A0_'+str(ptedges[0])+'_'+str(ptedges[1])+'_'+str(yedges[0])+'_'+str(yedges[1])
    smooth_rms_string='A8_'+str(ptedges[0])+'_'+str(ptedges[1])+'_'+str(yedges[0])+'_'+str(yedges[1])

    #update according to var
    if var=='eta':
	raw_string += 'lep_eta'
	smooth_string += 'lep_eta'
	boot_string += 'eta'
	smooth_rms_string += 'eta'
    elif var=='pt':
	raw_string += 'lep_pt'
	smooth_string += 'lep_pt'
	boot_string += 'pt'
	smooth_rms_string += 'pt'
    elif var=='phi':
	raw_string += 'phics'
	smooth_string += 'phics'
	boot_string += 'phi'
	smooth_rms_string += 'phi'
    elif var=='cost':
	raw_string += 'costcs'
	smooth_string += 'costcs'
	boot_string += 'cost'
	smooth_rms_string += 'cost'
    elif var=='gen2d':
        raw_string += 'gen2d'
        smooth_string += 'gen2d'
        boot_string += 'gen2d_rms'
        smooth_rms_string += 'gen2d_rms'      
    elif var=='reco2d':
        raw_string += 'reco2d'
        smooth_string += 'reco2d'
        boot_string += 'reco2d_rms'
        smooth_rms_string += 'reco2d_rms'
    elif var=='bos2d':
        raw_string += 'bos2d'
        smooth_string += 'bos2d'
        boot_string += 'bos2d_rms'
        smooth_rms_string += 'bos2d_rms'

    #print raw_string, smooth_string, boot_string, smooth_rms_string


    #get the hists
    raw= file0.Get(raw_string) 
    smooth=file0.Get(smooth_string)
    rms=file0.Get(smooth_rms_string)
    boot=file0.Get(boot_string)

    print raw, smooth, rms, boot

    new=raw.ProjectionX('name',1,5)

    print new

    return [raw, smooth, rms, boot]


def project(h2d, xaxis=True, out=False):
    '''
    Return a projection on x or y of the TH2F obj h2d
    '''

    #get projection based on axis requested
    if xaxis:
	h1d=h2d.ProjectionX('name',1,5)
    else:
	h1d=h2d.ProjectionY('name',1,5)

    #print to a test file if requested
    if out:
	c1=R.TCanvas()
	h1d.Draw()
	c1.Print('projection.pdf')

    return h1d


def unroll(h2d,var,typestring,choosebin=False):
    '''
    Function to unroll a 2d histogram into a 1d one
    This PyRoot version translated from the cpp form in TGHistogramUtils::unrollHistogram

    Returns, h1d, a TH1F
    '''

    #get all the information about this 2d hist
    nx=h2d.GetXaxis().GetNbins()
    ny=h2d.GetYaxis().GetNbins()


    #determine the hist name and title from var
    if var=='gen2d':
	name=var+typestring+'1d'
	title='Lepton cos(#theta) and \phi; ; #frac{d\sigma}{dp_{T}dydcos(#theta)d\phi}'
	label='cos(#theta) \in '
	binsize=int( 100* 2/20)/100
	low,high = -1, 1

    elif var=='reco2d':
	name=var+typestring+'1d'
	title='Lepton p_{T} and \eta; ; d\sigma/dp_{T}dydcos(#theta)d\phi'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low,high = 0,120

    elif var=='bos2d':
	name=var+typestring+'1d'
	title='Boson p_{T} and y; ; d\sigma/dp_{T}dydcos(#theta)d\phi'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low,high = 0, 120

    elif var=='dress2d':
	name=var+typestring+'1d'
	title='Dressed lepton p_{T} and \eta; ; d\sigma/dp_{T}dydcos(#theta)d\phi'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low,high = 0,120

    elif var=='cost':
	name=var+typestring+'1d'
	title='Cost; ; d\sigma/dp_{T}dydcos(#theta)d\phi'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low,high = -1,1

    elif var=='phi':
	name=var+typestring+'1d'
	title='phi; ; d\sigma/dp_{T}dydcos(#theta)d\phi'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low,high = 0,2*math.pi

    #limits must be changed for distribution histograms
    if var=='cost' and typestring=='distribution':
	name=var+typestring+'1d'
	low,high = 1,ny+1
	title='resample distribution for cos(#theta); unrolled along cos(#theta) bin number; count '

    elif var=='phi' and typestring=='distribution':
	name=var+typestring+'1d'
	low,high = 1,ny+1
	title='resample distribution for \phi; unrolled along \phi bin number; count '

    elif var=='eta' and typestring=='distribution':
	name=var+typestring+'1d'
	low,high = 1,ny+1
	title='resample distribution for lepton \eta; unrolled along \eta bin number; count '

    elif var=='pt' and typestring=='distribution':
	name=var+typestring+'1d'
	low,high = 1,ny+1
	title='resample distribution for lepton p_{T}; unrolled along p_{T} bin number; count '


    if var=='cost' and typestring=='distribution' and choosebin:
	name=var+typestring+'1d'
	low,high = 2*10**8,4*10**8
	title='resample distribution for cos(#theta); unrolled along cos(#theta) bin number; count '

    elif var=='phi' and typestring=='distribution' and choosebin:
	name=var+typestring+'1d'
	low,high = 2*10**8,4*10**8
	title='resample distribution for \phi; unrolled along \phi bin number; count '

    elif var=='eta' and typestring=='distribution' and choosebin:
	name=var+typestring+'1d'
	low,high = 2*10**8,4*10**8
	title='resample distribution for lepton \eta; unrolled along \eta bin number; count '

    elif var=='pt' and typestring=='distribution' and choosebin:
	name=var+typestring+'1d'
	low,high = 0,4*10**8
	title='resample distribution for lepton p_{T}; unrolled along p_{T} bin number; count '


    #make a new 1d hist
    h1d=R.TH1F(name,title,nx*ny,low,high) 

    #fill and label new bins
    bincount=1

    #iter over y
    for j in range(1,ny+1):
	
	#iter over x
	for i in range(1,nx+1):

	    #fill the 1d bin with contents of corresponding 2d binn
	    h1d.SetBinContent(bincount,h2d.GetBinContent(i,j))
	    h1d.SetBinError(bincount,h2d.GetBinError(i,j))
	    bincount += 1

    print bincount


    #format hist based on typestring
    if typestring=='raw':
	h1d.SetFillColorAlpha(R.kGreen, 0.6)
	h1d.SetLineColorAlpha(R.kGreen, 0.6)
    elif typestring=='rms':
	h1d.SetLineColor(R.kBlack)
    elif typestring=='boot':
	h1d.SetFillColorAlpha(R.kAzure +2, 0.8)
	h1d.SetLineColorAlpha(R.kAzure+2, 0.8)
    elif typestring=='distribution':
	h1d.GetXaxis().SetNdivisions(ny)

    return h1d


def condense(oldhist,var,thresh=0):
    '''
    trims an unrolled hist by removing the bins that have zero or negative contents, or contents
    below a certain threshhold. Converts to a tgraph in the process

    Args:
    oldhist, TH1 that we made from unrolling a reco2d, dress2d th2
    var, tells which var the hist is of, to determine threshold

    Returns: an Agraph obj (see treesuite)
    '''

    #determine thresh based on var
    if var=='gen2d':
	thresh=0
    elif var=='reco2d' or var=='dress2d':
	thresh=2*10**6
    elif var=='bos2d':
	thresh=10**6

    #init arrays for the new graph
    xvals=array.array('f')
    yvals=array.array('f')
    xerrs=array.array('f')
    yerrs=array.array('f')

    n_bins=oldhist.GetXaxis().GetNbins()

    for i in range(1,n_bins+1):

	#we exlude contents below the threshhold
	content=oldhist.GetBinContent(i)
	if content > thresh:

	    #copy everything to the tgraph arrays
	    xvals.append(oldhist.GetBinCenter(i) )
	    yvals.append(content)
	    xerrs.append(0)
	    yerrs.append(oldhist.GetBinError(i) )

    #returns an Agraph object
    return t.Agraph( xvals,yvals,xerrs,yerrs, 'none','none','title')



def project3(raw,rms,boot,var,x=True,hist=True):
    '''
    Copare the main 3 hists as projections, either x axis or y axis
    '''

    #get either x or y projection
    if x:

	raw_p=raw.ProjectionX('raw_x',1,20,'e')
	rms_p=rms.ProjectionX('rms_x',1,20,'e')
	boot_p=boot.ProjectionX('boot_x',1,20,'e')

    else:

	raw_p=raw.ProjectionY('raw_y',1,20,'e')
	rms_p=rms.ProjectionY('rms_y',1,20,'e')
	boot_p=boot.ProjectionY('boot_y',1,20,'e')

    #format
    raw_p.SetLineColorAlpha(R.kGreen,0.6)
    rms_p.SetLineColor(R.kBlack)
    boot_p.SetLineColorAlpha(R.kAzure+2,0.8)
    raw_p.SetFillColorAlpha(R.kGreen,0.6)
    rms_p.SetFillColor(R.kBlack)
    boot_p.SetFillColorAlpha(R.kAzure+2,0.8)

    if hist:
	return raw_p,rms_p,boot_p

    #make canvas and draw
    c1=R.TCanvas()
    raw_p.Draw('e3')
    rms_p.Draw('same')
    boot_p.Draw('e3same')
    c1.Print(var+'1d.pdf')



def comp_plot(raw,rms,boot,var,index,ptedges,yedges,total_resamples,inc=False):
    '''
    Makes plot of errors as well as ratio plot comparing errors
    The raw,rms,boot args are Agraph objects now
    '''

    #make canvas and first pad
    c1=R.TCanvas()
    p1=R.TPad('p1','p1',0.05,0.4,1,1)
    p1.SetGrid()
    p1.Draw()
    p1.cd()

    #clean up decimals of pt, y lists
    for mylist in [ptedges,yedges]:
	for i in range(len(mylist )):
	    mylist[i]= str(mylist[i])[:4] 

    #determine the hist name and title from var
    if var=='gen2d':
	xtitle=str(total_resamples)+' resamples       cos(#theta) bin limits shown       cross section unrolled along \phi, \phi \in [0,2\pi]'
	ytitle=' #frac{d\sigma}{dp_{T}dydcos(#theta)d\phi}'
	title='W cos(#theta) and \phi (p_{T}^{W} \in '+str(ptedges)+', y^{W} \in '+str(yedges)+')'
	label='cos(#theta) \in '
	binsize=int( 100* 2/20)/100
	low,high = -1, 1

    if var=='reco2d':
	xtitle=str(total_resamples)+' resamples       Lepton p_{T} bin limits shown       cross section unrolled along lepton \eta, \eta \in [-2.4,2.4]'
	title='Lepton p_{T} and \eta (p_{T}^{W} \in '+str(ptedges)+', y^{W} \in '+str(yedges)+')'
	ytitle=' #frac{d\sigma}{dp_{T}dydcos(#theta)d\phi}'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low,high = 0, 120

    if var=='dress2d':
	xtitle=str(total_resamples)+' resamples       Lepton p_{T} bin limits shown       cross section unrolled along lepton \eta, \eta \in [-2.4,2.4]'
	title='Dressed lepton p_{T} and \eta (p_{T}^{W} \in '+str(ptedges)+', y^{W} \in '+str(yedges)+')'
	ytitle=' #frac{d\sigma}{dp_{T}dydcos(#theta)d\phi}'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low,high = 0, 120

    if var=='bos2d':
	xtitle=str(total_resamples)+' resamples       Lepton p_{T} bin limits shown       cross section unrolled along W y, y \in '+str(yedges)
	ytitle=' #frac{d\sigma}{dp_{T}dydcos(#theta)d\phi}'
	title='W p_{T} and y'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low,high = float(ptedges[0]), float(ptedges[1])

    if var=='cost':
	xtitle='\cost \in [-1,1]'
	title='W cos(#theta) (p_{T}^{W} \in '+str(ptedges)+', y^{W} \in '+str(yedges)+')'
	ytitle=' #frac{d\sigma}{dp_{T}dydcos(#theta)d\phi}'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low,high = -1, 1

    if var=='phi':
	xtitle='\phi \in [0,2\pi]'
	title='W \phi (p_{T}^{W} \in '+str(ptedges)+', y^{W} \in '+str(yedges)+')'
	ytitle=' #frac{d\sigma}{dp_{T}dydcos(#theta)d\phi}'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low,high = 0, 2*math.pi

    #format the main graph
    rms.graph.GetXaxis().SetNdivisions( 20 )
    rms.graph.GetXaxis().SetLimits(0,60)

    #format main graph title
    rms.graph.GetXaxis().SetTitle(xtitle)
    rms.graph.GetYaxis().SetTitle(ytitle)
    rms.graph.SetTitle(title)
    rms.graph.GetXaxis().SetTitleSize(0.05)
    rms.graph.GetXaxis().SetTitleOffset(0.9)
    rms.graph.GetYaxis().SetTitleSize(0.05)
    rms.graph.GetYaxis().SetTitleOffset(0.6)

    #draw hists on upper pad
    rms.graph.Draw('ape')
    raw.graph.Draw('3same')
    boot.graph.Draw('3same')

    #make error graohs
    raw_g=R.TGraph(len(raw.xvals),raw.xvals,raw.get_rels() )
    rms_g=R.TGraph(len(rms.xvals),rms.xvals,rms.get_rels() )
    boot_g=R.TGraph(len(boot.xvals),boot.xvals,boot.get_rels() )

    #color format all graphs
    raw.graph.SetLineColorAlpha(R.kGreen,0.6)
    raw.graph.SetFillColorAlpha(R.kGreen,0.6)
    raw_g.SetLineColorAlpha(R.kGreen,0.6)
    raw_g.SetFillColorAlpha(R.kGreen,0.6)

    rms.graph.SetLineColor(R.kBlack)  
    rms_g.SetLineColor(R.kBlack)

    boot.graph.SetLineColorAlpha(R.kAzure+2,0.6)
    boot.graph.SetFillColorAlpha(R.kAzure+2,0.6)  
    boot_g.SetLineColorAlpha(R.kAzure+2,0.6)
    boot_g.SetFillColorAlpha(R.kAzure+2,0.6)

    #make the strings nicer
    if var=='gen2d':
	nicevar='W \phi and cos(#theta)'
    if var=='reco2d' or 'dress2d':
	nicevar='Lepton \eta and p_{T}'
    if var=='bos2d':
	nicevar='W y and p_{T}'
    if var=='cost':
	nicevar='Lepton cos(#theta)'
    if var=='phi':
	nicevar='Lepton \phi'
	    
    #make the legend
    leg1=R.TLegend(0.7,0.7,0.9,0.9)
    leg1.AddEntry(raw_g,'Monte Carlo '+nicevar,'f' )
    leg1.AddEntry(boot_g,'Bootstrapped '+nicevar ,'f')
    leg1.AddEntry(rms_g,'Smooth bootstrapped '+nicevar,'l')
    leg1.Draw()

    #make the ratio plot
    c1.cd()

    #lower tpad
    p2=R.TPad('p2','p2',0.05,0,1,0.4)
    p2.SetGrid()
    p2.Draw()
    p2.cd()

    #format error graphs
    raw_g.GetXaxis().SetNdivisions(20)
    raw_g.GetXaxis().SetLimits(0,60)

    #format error graph titles
    raw_g.SetTitle()
    raw_g.GetYaxis().SetTitle(' % error ')
    raw_g.GetYaxis().SetTitleSize(0.07)
    raw_g.GetYaxis().SetTitleOffset(0.4)

    #draw error graphs
    raw_g.Draw('alsame')
    rms_g.Draw('lsame')
    boot_g.Draw('lsame')


    #print canvas depending on job index
    if inc:
	c1.Print(var+'_inc.pdf')

    else:
        c1.Print(var+str(index)+'.pdf')




def old_comp_plot(raw,rms,boot,var,index,ptedges,yedges,total_resamples,inc=False):
    '''
    Makes plot of errors as well as ratio plot comparing errors
    '''


    #make canvas and first pad
    c1=R.TCanvas()
    p1=R.TPad('p1','p1',0.05,0.4,1,1)
    p1.SetGrid()
    p1.Draw()
    p1.cd()

    #clean up decimals of pt, y lists
    for mylist in [ptedges,yedges]:
	for i in range(len(mylist )):
	    mylist[i]= str(mylist[i])[:4] 

    #info from the graphs
    n_bins=raw.GetXaxis().GetNbins()

    print n_bins

    #determine the hist name and title from var
    if var=='gen2d':
	xtitle=str(total_resamples)+' resamples       cos(#theta) bin limits shown       cross section unrolled along \phi, \phi \in [0,2\pi]'
	ytitle=' #frac{d\sigma}{dp_{T}dydcos(#theta)d\phi}'
	title='W cos(#theta) and \phi (p_{T}^{W} \in '+str(ptedges)+', y^{W} \in '+str(yedges)+')'
	label='cos(#theta) \in '
	binsize=int( 100* 2/20)/100
	low, high = -1, 1
	err_low, err_high = 0,40

    if var=='lep_eta':
	xtitle='Lepton \eta \in [-2.4,2.4]'
	ytitle=' d\sigma/dp_{T}dydcos(#theta)d\phi'
	title='Lepton \eta, (p_{T}^{W} \in '+str(ptedges)+', y^{W} \in '+str(yedges)+')'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low, high = -2.4,2.4
	err_low, err_high = 0,40

    if var=='lep_pt':
	xtitle='Lepton p_{T} \in [0,60]'
	ytitle=' d\sigma/dp_{T}dydcos(#theta)d\phi'
	title='Lepton p_{T}, (p_{T}^{W} \in '+str(ptedges)+', y^{W} \in '+str(yedges)+')'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low, high = 0,60
	err_low, err_high = 0,40	


    if var=='cost':
	xtitle='\cost \in [-1,1]'
	ytitle=' d\sigma/dp_{T}dydcos(#theta)d\phi'
	title='W cos(#theta) (p_{T}^{W} \in '+str(ptedges)+', y^{W} \in '+str(yedges)+')'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low, high = -1, 1

	#color format the hists as well
        raw.SetFillColorAlpha(R.kGreen,0.6)
        rms.SetLineColor(R.kBlack)
        boot.SetFillColorAlpha(R.kAzure+2,0.6)

    if var=='phi':
	xtitle='\phi \in [0,2\pi]'
	ytitle=' d\sigma/dp_{T}dydcos(#theta)d\phi'
	title='W \phi (p_{T}^{W} \in '+str(ptedges)+', y^{W} \in '+str(yedges)+')'
	label='p_{T} \in '
	binsize=int( 100* 120/20)/100
	low, high = 0, 2*math.pi

	#color format the hists as well
        raw.SetFillColorAlpha(R.kGreen,0.6)
        rms.SetLineColor(R.kBlack)
        boot.SetFillColorAlpha(R.kAzure+2,0.6)

    #format the main graph
    rms.GetXaxis().SetNdivisions( 20 )
    #rms.GetXaxis().SetLimits(low, high)

    #format main graph title
    rms.GetXaxis().SetTitle(xtitle)
    rms.SetTitle(title)
    rms.GetXaxis().SetTitleSize(0.05)
    rms.GetXaxis().SetTitleOffset(0.9)
    rms.GetYaxis().SetTitleSize(0.05)
    rms.GetYaxis().SetTitleOffset(0.6)

    print 'getting rms xvals from old_comp_plot'
    for i in range(1,n_bins+1):
	print rms.GetBinCenter(i)


    #draw hists on upper pad
    rms.Draw()
    raw.Draw('E3same')
    boot.Draw('E3same')

    #get errors from hists and put them in arrays
    xvals=array.array('f')
    raw_y=array.array('f')
    boot_y=array.array('f')
    rms_y=array.array('f')

    #iter thru hist bins
    for i in range(1,1+n_bins):
	xvals.append(raw.GetBinCenter(i))
	raw_y.append(100*raw.GetBinError(i)/(raw.GetBinContent(i)+1))
	boot_y.append(100*boot.GetBinError(i)/(boot.GetBinContent(i)+1))
	rms_y.append(100*rms.GetBinError(i)/(rms.GetBinContent(i)+1))


    #make error graohs
    raw_g=R.TGraph(len(xvals),xvals,raw_y)
    rms_g=R.TGraph(len(xvals),xvals,rms_y)
    boot_g=R.TGraph(len(xvals),xvals,boot_y)

    print xvals

    #format error graphs
    raw_g.SetLineColorAlpha(R.kGreen,0.6)
    raw_g.SetFillColorAlpha(R.kGreen,0.6)
    rms_g.SetLineColor(R.kBlack)
    boot_g.SetLineColorAlpha(R.kAzure+2,0.6)
    boot_g.SetFillColorAlpha(R.kAzure+2,0.6)

    #format error graph titles
    raw_g.GetXaxis().SetNdivisions( 20 )
    #raw_g.GetXaxis().SetLimits(low, high)
    raw_g.GetHistogram().SetMinimum(err_low)
    raw_g.GetHistogram().SetMaximum(err_high)

    #make the strings nicer
    if var=='gen2d':
	nicevar='W \phi and cos(#theta)'
    elif var=='cost':
	nicevar='W cos(#theta)'
    elif var=='phi':
	nicevar='W \phi'
    elif var=='lep_eta':
	nicevar='Lepton \eta'
    elif var=='lep_pt':
	nicevar='Lepton p_{T}'
	    
    #make the legend
    leg1=R.TLegend(0.7,0.7,0.9,0.9)
    leg1.AddEntry(raw_g,'Monte Carlo '+nicevar,'f' )
    leg1.AddEntry(boot_g,'Bootstrapped '+nicevar,'f' )
    leg1.AddEntry(rms_g,'Smooth bootstrapped '+nicevar,'l')
    leg1.Draw()

    #make the ratio plot
    c1.cd()

    #lower tpad
    p2=R.TPad('p2','p2',0.05,0,1,0.4)
    p2.SetGrid()
    p2.Draw()
    p2.cd()

    #format error graph titles
    raw_g.SetTitle()
    raw_g.GetYaxis().SetTitle(' % error ')
    raw_g.GetYaxis().SetTitleSize(0.07)
    raw_g.GetYaxis().SetTitleOffset(0.4)

    #draw error graphs
    raw_g.Draw('al')
    rms_g.Draw('lsame')
    boot_g.Draw('lsame')

    #print canvas depending on job index
    if inc:
	c1.Print(var+'_inc.pdf')

    else:
        c1.Print(var+str(index)+'.pdf')




def A_plot(normgraph, rmsgraph, normratio, rmsratio, A_index,ptedges,yedges):
    '''
    Makes unrolled hists of A_i with ratio plot of errors
    '''

    #make canvas and first pad
    c1=R.TCanvas()
    p1=R.TPad('p1','p1',0.05,0.35,1,1)
    p1.SetGrid()
    p1.SetLogx()
    p1.Draw()
    p1.cd()

    #format colors
    rmsgraph.SetFillColorAlpha(R.kAzure,0.6)	
    rmsratio.SetFillColorAlpha(R.kAzure,0.6)
    rmsgraph.SetLineColor(R.kAzure)
    rmsratio.SetLineColor(R.kAzure)
    normgraph.SetLineColor(R.kBlack)
    normratio.SetLineColor(R.kBlack)

    #draw hists on upper pad
    normgraph.Draw('ape')
    rmsgraph.Draw('3same')

    #format axes
    normgraph.SetTitle(' A_{'+str(A_index)+'} vs. p_{T}^{W}, |y^{W}|; unrolled along p_{T}^{W}, p_{T}^{W} \in '+str([ptedges[0],ptedges[-1]])+' ; A_{'+str(A_index)+'} ')
    normgraph.GetXaxis().SetNdivisions(4)
    normgraph.GetXaxis().SetTitleSize(0.06)
    normgraph.GetXaxis().SetTitleOffset(0.7)
    normgraph.GetYaxis().SetTitleSize(0.06)
    normgraph.GetYaxis().SetTitleOffset(0.6)
    normgraph.GetYaxis().CenterTitle()
    
	    
    #make the legend
    leg1=R.TLegend(0.7,0.1,0.9,0.3)
    leg1.AddEntry(normratio,'Monte Carlo' ,'l' )
    leg1.AddEntry(rmsratio,'Bootstrapped','f' )
    leg1.Draw()

    #make the ratio plot
    c1.cd()

    #lower tpad
    p2=R.TPad('p2','p2',0.05,0,1,0.35)
    p2.SetGrid()
    p2.SetLogx()
    p2.Draw()
    p2.cd()

    # draw graphs
    normratio.Draw('al')
    rmsratio.Draw('lsame')

    #format graphs
    normratio.SetTitle()
    normratio.GetYaxis().SetTitle('Error')
    normratio.GetYaxis().SetTitleSize(0.06)
    normratio.GetYaxis().SetTitleOffset(0.6)
    normratio.GetYaxis().CenterTitle()
    normratio.GetXaxis().SetNdivisions(4)

    #print canvas
    c1.Print('A'+str(A_index)+'unrolled.pdf')





################################ Command line controls ################################


if __name__=='__main__':

    #get variables from the command line
    variables=['eta','pt','phi','cost','gen2d','reco2d','bos2d']
    functions=['smooth','inclusive','all','bootstrapped','ind']
    var=sys.argv[1]
    function=sys.argv[2]

    #screen out bad command line arguments
    if (not (var in variables) ) or (not (function in functions) ):

        #we want a purposeful error to hold the job
        print 'variable or function entry unsupported'
        d1= 'hello'+5

    #choose edges based on number of bins
    n_pt=8
    n_y=4

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


################################ Create seperate unrolled hists ################################

if __name__=='__main__':

    if function=='ind':

        #iter over the number of files there will be
        for i in range(n_pt*n_y):

	    #read the appropriate hists from the file
	    graphs=get4(var,i,ptedges,yedges)
	
	    #reassign
	    raw = graphs[0]
	    smooth = graphs[1]
	    rms = graphs[2]
	    boot = graphs[3]

	    print raw

	    raw.ProjectionX('name',1,5)

	    project(raw, out=True)









 

