import os, sys, operator
from pylab import *

#
# Get the psi values from a file of "old psi" computations.
#
def getSummary(summary_f):
    
    eventToInfo = {}
    for line in open(summary_f):
        if not line.startswith("#"):
            vals = line.strip().split("\t")
            event = vals[0]
            psi = vals[1]
            inc = vals[2]
            exc = vals[3]
            inclen = vals[4]
            exclen = vals[5]
            incdens = vals[6]
            excdens = vals[7]
            
            eventToInfo[event] = {}
            eventToInfo[event]['psi'] = float(psi)
            eventToInfo[event]['inclusion reads'] = int(inc) 
            eventToInfo[event]['exclusion reads'] = int(exc) 
            eventToInfo[event]['inclusion length'] = int(inclen)
            eventToInfo[event]['exclusion length'] = int(exclen)
            eventToInfo[event]['inclusion density'] = float(incdens)
            eventToInfo[event]['exclusion density'] = float(excdens)
        else:
            header = line[1:].strip().split("\t")

    return eventToInfo, header



# Scatter plot of psi values.
def psi1_vs_psi2(psi1_f, psi2_f, label1, label2, out_f):

    eventToInfo1, header1 = getSummary(psi1_f)   
    eventToInfo2, header1 = getSummary(psi2_f)   

    data = [] 
    events = [e for e in eventToInfo1 if e in eventToInfo2]
    for e in events: 
        data.append([eventToInfo1[e]['psi'], eventToInfo2[e]['psi'],\
            eventToInfo1[e]['inclusion reads'] + eventToInfo1[e]['exclusion reads'],\
            eventToInfo2[e]['inclusion reads'] + eventToInfo2[e]['exclusion reads']])
    
    data = array(data, dtype='f')

    figure(figsize=(6,6))
    a = pi / 4.
    R = array([[cos(a) , -sin(a)], [sin(a), cos(a)]])
    rotated = dot(vstack((data[:, 0], data[:, 1])).T, R)
    xvals = linspace(0, sqrt(2), 30)
    plusspread = []
    minusspread = []
    for x in range(len(xvals) - 1):
        idx = where((rotated[:, 0] >= xvals[x]) & \
            (rotated[:, 0] <= xvals[x + 1]))[0]
        plusspread.append([xvals[x:x+2].mean(), rotated[idx, 1].std()])
        minusspread.insert(0, [xvals[x:x+2].mean(), -rotated[idx, 1].std()])

    plusspread = array(plusspread)
    minusspread = array(minusspread)
    
    a = -pi / 4.
    R = array([[cos(a) , -sin(a)], [sin(a), cos(a)]])
    plusrotated = dot(plusspread, R)
    minusrotated = dot(minusspread, R)
    allrotated = vstack((plusrotated, minusrotated))  
 
    fill(allrotated[:, 0], allrotated[:, 1], color='#EEEEEE', zorder=0)
    
    scatter(data[:,0], data[:,1], s=1, color='k', alpha=.5)
 
    xlim(-.1,1.1)
    ylim(-.1,1.1) 
    xlabel(os.path.basename(psi1_f).split("_")[0] + ' $\psi$ ' + '(' + label1 + ')',\
        fontsize=12, horizontalalignment='center')
    ylabel(os.path.basename(psi2_f).split("_")[0] + ' $\psi$ ' + '(' + label2 + ')',\
        fontsize=12, horizontalalignment='center')
  
    savefig(out_f)



def psiTable(psi1_f, psi2_f, ensGeneMap_f, out_f, minct=10):

    import shelve
    eventToInfo = shelve.open(ensGeneMap_f,'r') 
 
    eventToInfo1, header1 = getSummary(psi1_f)   
    eventToInfo2, header1 = getSummary(psi2_f)   

    data = [] 
    events = [e for e in eventToInfo1 if e in eventToInfo2]
    for e in events: 
        data.append([e, eventToInfo1[e]['psi'], eventToInfo2[e]['psi'],\
            eventToInfo1[e]['inclusion reads'] + eventToInfo1[e]['exclusion reads'],\
            eventToInfo2[e]['inclusion reads'] + eventToInfo2[e]['exclusion reads'],\
            eventToInfo1[e]['psi'] - eventToInfo2[e]['psi']])
        
    out = open(out_f, 'w')
    out.write("\t".join(["#event_name",\
        psi1_f + "_psi", psi2_f + "_psi",\
        psi1_f + "_nreads", psi2_f + "_nreads",\
        "delta", "gene", "symb", "desc"]) + "\n")

    data.sort(key=operator.itemgetter(-1)) 
    for item in data:
        if min([item[-2], item[-3]]) >= minct:
            event = item[0]
            try:
                item.extend(eventToInfo[item[0]])
            except: 
                item.extend(["n/a","n/a","n/a"])
            out.write("\t".join(map(str,item))+"\n")
    out.close() 


