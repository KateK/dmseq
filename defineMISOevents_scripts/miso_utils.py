import os, sys, operator, re, glob, subprocess, shelve
from pylab import *

def plotMultiplePosteriors(*args):
    
    files = args[:-2]
    resolution = args[-2]    
    out_f = args[-1]

    n = 1
    figure(figsize=(8.5, 11))
    for f in files:
        subplot(int(ceil(len(files) / 3.)), 3, n) 
        plotPosterior(f, resolution)
        title(f.split("/")[0], fontsize=8)
        n += 1
    subplots_adjust(wspace=.5, hspace=.8)
    savefig(out_f)


def plotPosterior(miso_f, resolution):
  
    resolution = int(resolution) 
    psis = [] 
    for line in open(miso_f):
        if not line.startswith("#") and not line.startswith("sampled"):
            psi, logodds = line.strip().split("\t")
            psis.append(float(psi.split(",")[0]))
  
    ci = .95 
    alpha = 1 - ci
    lidx = int(round((alpha / 2) * len(psis)) - 1)
    # the upper bound is the (1-alpha/2)*n nth smallest sample, where n is
    # the number of samples
    hidx = int(round((1 - alpha / 2) * len(psis)) - 1)
    psis.sort()
    clow, chigh = [psis[lidx], psis[hidx]]
   
    hist(psis, linspace(0, 1, resolution), normed=True, facecolor='k',\
        edgecolor='w') 
    axvline(clow, linestyle='--', color='#CCCCCC')
    axvline(chigh, linestyle='--', color='#CCCCCC')
    axvline(median(psis), color='r')
    text(.95, 14, "$\Psi$ = %.2f\n$\Psi_{low}$ = %.2f\n$\Psi_{high}$ = %.2f"%(median(psis),\
        clow, chigh), fontsize=6, va='top', ha='right')

    ylim(0, 15)
    xlabel("Posterior $\Psi$", fontsize=8)
    ylabel("Probability", fontsize=8)


# Get information from a single summary (non-comparison)
def getSummarySingle(summary_f):
    
    eventToInfo = {}
    for line in open(summary_f):
        if not line.startswith("event_name"):
            try:
                vals = line.strip().split("\t")
                event = vals[0]
                s1m = vals[1]
                s1l = vals[2]
                s1h = vals[3]
                isoforms = vals[4]
                s1c = vals[5]
                s1ac = vals[6]
                s1ac = sum([int(x.split(":")[1]) for x in s1ac.split(",")])
                chrom = vals[7]
                strand = vals[8]
                starts = vals[9]
                ends = vals[9]
                
                eventToInfo[event] = {} 
                eventToInfo[event]['mean'] = map(float, s1m.split(","))
                eventToInfo[event]['low'] = map(float, s1l.split(","))
                eventToInfo[event]['high'] = map(float, s1h.split(","))
                eventToInfo[event]['isoforms'] = isoforms
                eventToInfo[event]['counts'] = s1c
                eventToInfo[event]['assigned counts'] = s1ac
                eventToInfo[event]['chrom'] = chrom 
                eventToInfo[event]['strand'] = strand 
                eventToInfo[event]['starts'] = starts 
                eventToInfo[event]['ends'] = ends 
            except:
                pass
        else:
            header = line[1:].strip().split("\t")

    return eventToInfo, header


# Function to retrieve information from a summarized comparison.
def getSummary(summary_f):
    
    eventToInfo = {}
    for line in open(summary_f):
        if not line.startswith("#"):
            try:
                vals = line.strip().split("\t")
                event = vals[0]
                s1m = vals[1]
                s1l = vals[2]
                s1h = vals[3]
                s2m = vals[4]
                s2l = vals[5]
                s2h = vals[6]
                dpsi = vals[7]
                bf = vals[8]
                isoforms = vals[9]
                s1c = vals[10]
                s1ac = vals[11]
                s1ac = sum([int(x.split(":")[1]) for x in s1ac.split(",")])
                s2c = vals[12]
                s2ac = vals[13]
                s2ac = sum([int(x.split(":")[1]) for x in s2ac.split(",")])
                maxbf = vals[18]
                gene = vals[19]
                symb = vals[20]
                desc = vals[21]
                
                eventToInfo[event] = {} 
                eventToInfo[event]['sample 1 mean'] = map(float, s1m.split(","))
                eventToInfo[event]['sample 1 low'] = map(float, s1l.split(","))
                eventToInfo[event]['sample 1 high'] = map(float, s1h.split(","))
                eventToInfo[event]['sample 2 mean'] = map(float, s2m.split(","))
                eventToInfo[event]['sample 2 low'] = map(float, s2l.split(","))
                eventToInfo[event]['sample 2 high'] = map(float, s2h.split(","))
                eventToInfo[event]['delta psi'] = map(float, dpsi.split(","))
                eventToInfo[event]['bayes factor'] = map(float, bf.split(","))
                eventToInfo[event]['isoforms'] = isoforms
                eventToInfo[event]['sample 1 counts'] = s1c
                eventToInfo[event]['sample 1 assigned counts'] = s1ac
                eventToInfo[event]['sample 2 counts'] = s2c
                eventToInfo[event]['sample 2 assigned counts'] = s2ac
                eventToInfo[event]['max bf'] = float(maxbf)
                eventToInfo[event]['gene'] = gene
                eventToInfo[event]['symb'] = symb
                eventToInfo[event]['desc'] = desc
            except:
                print event
        else:
            header = line[1:].strip().split("\t")

    return eventToInfo, header


# Helper function to parse counts field.
def parseCounts(countsfield):

    cttot = 0
    iterator = re.finditer("\(.+?\)\:[0-9]+", countsfield)
    for item in iterator:
        isoform, ct = countsfield[item.start() : item.end()].split(":")
        cttot += int(ct)
    return cttot


def plotDistributions(bf_f, out_f):
    
    eventToInfo, header = getSummary(bf_f)
    
    fields = ['sample 1 mean', 'sample 2 mean', 'delta psi',\
        'max bf', 'sample 1 assigned counts', 'sample 2 assigned counts']
    logged = [False, False, False, False, True, True]
    islist = [True, True, True, False, False, False]
    mins = [0, 0, -1, -1, 0, 0]
    maxs = [1, 1, 1, 2, 5, 5] 

    figure(figsize=(11, 14))
    n = 1
    for i in range(len(fields)):
        for j in range(len(fields)):
            subplot(len(fields), len(fields), n) 

            if islist[i]:
                xvals = array([eventToInfo[e][fields[i]][0] for e in eventToInfo])
            else:
                xvals = array([eventToInfo[e][fields[i]] for e in eventToInfo])
            if islist[j]:
                yvals = array([eventToInfo[e][fields[j]][0] for e in eventToInfo])
            else:
                yvals = array([eventToInfo[e][fields[j]] for e in eventToInfo])
            if logged[i]:
                xvals = log10(xvals + 1)
            if logged[j]:
                yvals = log10(yvals + 1)
            xidx = where((xvals >= mins[i])&(xvals <= maxs[i]))[0]
            yidx = where((yvals >= mins[j])&(yvals <= maxs[j]))[0]
            idx = array([k for k in xidx if k in yidx])
            xvals = xvals[idx]
            yvals = yvals[idx]
            #scatter(xvals, yvals, s=1, alpha=.2, rasterized=True)
            hexbin(xvals, yvals, bins='log')#, s=1, alpha=.2, rasterized=True)

            xlabel(fields[i], fontsize=6)
            ylabel(fields[j], fontsize=6)
            xticks(fontsize=6, rotation=60)
            yticks(fontsize=6)
            xlim(mins[i], maxs[i])
            ylim(mins[j], maxs[j])

            n += 1
    
    subplots_adjust(hspace=1.0, wspace=1.0)
    savefig(out_f, dpi=150)



def countSigEvents(settings_f):
    
    files = [] 
    for line in open(settings_f):
        files.append(line.strip())

    allEvents = {} 
    for f in files:
        eventToInfo, header = getSummary(f)
        for e in eventToInfo:
            bf = eventToInfo[e]['max bf'] 
            dpsi = eventToInfo[e]['delta psi'][0]
            if bf >= 5 and abs(dpsi) >= 0.05:    
                allEvents[e] = 0
   
    print len(allEvents)



# Scatter plot of psi values.
def oldpsi_vs_MISOpsi(oldpsi_f, bf_f, bfidx, label1, label2, out_f):

    import psi_utils
    eventToInfo1, header1 = psi_utils.getSummary(oldpsi_f)   
    eventToInfo2, header2 = getSummary(bf_f)   

    data = [] 
    events = [e for e in eventToInfo1 if e in eventToInfo2]
    for e in events: 
        data.append([eventToInfo1[e]['psi'],\
            eventToInfo2[e]['sample ' + bfidx + ' mean']])
    
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
    xlabel(os.path.basename(oldpsi_f).split("_")[0] + ' $\psi$ ' + '(' + label1 + ')',\
        fontsize=12, horizontalalignment='center')
    ylabel(header2[(int(bfidx) - 1) * 3 + 1].split("_")[0] + ' $\psi$ ' + '(' + label2 + ')',\
        fontsize=12, horizontalalignment='center')
  
    savefig(out_f)


# Scatter plot of psi values.
# idx1 and idx2 are 1 or 2.
def psi1_vs_psi2(bf1_f, bf2_f, idx1, idx2, label1, label2, out_f=False):

    eventToInfo1, header1 = getSummary(bf1_f)   
    eventToInfo2, header2 = getSummary(bf2_f)   

    fields = ['sample ' + str(idx1) + ' mean',\
        'sample ' + str(idx2) + ' mean']
      
    data = [] 
    events = [e for e in eventToInfo1 if e in eventToInfo2]
    for e in events: 
        data.append([eventToInfo1[e]['sample ' + idx1 + ' mean'][0],\
            eventToInfo2[e]['sample ' + idx2 + ' mean'][0]])
    
    data = array(data, dtype='f')
    print data.shape

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
    xlabel('$\psi$, ' + label1,\
        fontsize=12, horizontalalignment='center')
    ylabel('$\psi$, ' + label2,\
        fontsize=12, horizontalalignment='center')
 
    if out_f is not False: 
        savefig(out_f)
    else:
        return corrcoef(data[:, 0], data[:, 1])[0][1]


# Plot number of significant events vs. read coverage.
# Estimate an asymptote.
def plotAsymptote(bf_f, outdir):

    print "Getting BF info."
    eventToInfo, header = getSummary(bf_f)
  
    data = [] 
    for event in eventToInfo:
        bf = eventToInfo[event]['bayes factor']
        s1c = eventToInfo[event]['sample 1 counts']
        s2c = eventToInfo[event]['sample 2 counts']

        ct1 = 0
        iterator = re.finditer("\(.+?\)\:[0-9]+", s1c)
        for item in iterator:
            isoform, ct = s1c[item.start() : item.end()].split(":")
            ct1 += int(ct)

        ct2 = 0
        iterator = re.finditer("\(.+?\)\:[0-9]+", s2c)
        for item in iterator:
            isoform, ct = s2c[item.start() : item.end()].split(":")
            ct2 += int(ct)

        data.append([log10(ct1 + ct2), max(bf)])

    print "Processing."
    data.sort(key=operator.itemgetter(0))
    data = array(data)
    binres = 100 
    min_bf = 5
   
    plotdata = [] 
    for i in range(data.shape[0]):
        if i%binres==0:
            window = arange(max([0, i - binres]),\
                min([data.shape[0] - 1, i + binres]))
            xval = mean(data[window, 0])
            yval = len(where(data[window, 1] >= min_bf)[0]) /\
                float(len(window))
            plotdata.append([xval, yval])
    
    plotdata = array(plotdata, dtype='f')
    plot(plotdata[:, 0], plotdata[:, 1], 'k.')
    savefig(os.path.join(outdir, os.path.basename(bf_f) + ".asymptote.pdf"))


# Output psi means and stdevs from posteriors (sampled values).
def psiFromPosteriors(posteriorDir, outDir):
 
    samples = os.listdir(posteriorDir)
    for sample in samples:
        print sample
        sampleDir = os.path.join(posteriorDir, sample)
        out_f = os.path.join(outDir, sample)
        out = open(out_f, 'w')
        out.write("#Event\tMean\tSD\tN\n")

        chroms = os.listdir(sampleDir)
        for chrom in chroms:
            if chrom.startswith("chr"):
                print chrom 
                chromDir = os.path.join(sampleDir, chrom)
                files = os.listdir(chromDir)
                for f in files:
                    ename = f.split(".")[0]
                    fname = os.path.join(chromDir, f)
                    i = 0
                    psivals = []
                    for line in open(fname):
                        if i > 1:
                            psis, score = line.strip().split()
                            psi = float(psis.split(",")[0])
                            psivals.append(psi)
                        i += 1

                    if len(psivals) > 100:
                        out.write("\t".join(map(str, [ename, mean(psivals),\
                            std(psivals), len(psivals)])) + "\n")

    out.close() 


# Output a psi table from psi files.
def psiTableFromPsiFiles(psiDir, order_f, lookup_f, out_f):
  
    eventToInfo = shelve.open(lookup_f, 'r')
 
    eventToPsi = {} 
    samples = []
    for line in open(order_f):
        samples.append(line.strip())
    for sample in samples:
        psi_f = os.path.join(psiDir, sample)
        for line in open(psi_f):
            if not line.startswith("#"):
                event, m, sd, n = line.strip().split("\t")
                if event not in eventToPsi:
                    eventToPsi[event] = {}
                eventToPsi[event][sample] = float(m)

    out = open(out_f, 'w')
    out.write("#Event\t")
    out.write("\t".join(samples) + "\n")
    for event in eventToPsi:
        out.write(event + "\t")
        psimeans = []
        for sample in samples:
            if sample in eventToPsi[event]:
                psimeans.append(str(round(eventToPsi[event][sample], 2)))
            else:
                psimeans.append('n/a')
        if event in eventToInfo:
            gene, symb, desc = eventToInfo[event] 
        else:
            gene, symb, desc = ['n/a', 'n/a', 'n/a']
        out.write("\t".join(psimeans) + "\t")
        out.write(gene + "\t" + symb + "\t" + desc + "\n")

    out.close() 

# Plot psi values for various time points.  Use a settings file to 
# define x-axis values for each time point.
def psiTimeCourse(posteriorDir, settings_f, event, out_f):
  
    samples = []
    sampleToX = {} 
    sampleToColor = {}
    for line in open(settings_f):
        sample, xval, color = line.strip().split()
        samples.append(sample) 
        sampleToX[sample] = float(xval)
        sampleToColor[sample] = color
   
    sampleToY = {} 
    for sample in samples:
        print sample
        if sample in os.listdir(posteriorDir):
            chrom = event.split(":")[0]
            fname = event + ".miso" 
            if fname in os.listdir(os.path.join(posteriorDir, sample, chrom)):
                posterior_f = os.path.join(posteriorDir, sample, chrom, fname)
                i = 0
                psivals = []
                for line in open(posterior_f):
                    if i > 1:
                        psis, score = line.strip().split()
                        psi = float(psis.split(",")[0])
                        psivals.append(psi)
                    i += 1
                if len(psivals) > 100:
                    sampleToY[sample] = psivals
      
    figure(figsize=(3, 3)) 
    xToMeans = {}
    for sample in samples:
        if sample in sampleToX and sample in sampleToY:
            try:
                xToMeans[sampleToX[sample]].append(mean(sampleToY[sample]))
            except:
                xToMeans[sampleToX[sample]] = [mean(sampleToY[sample])]
            errorbar(sampleToX[sample], mean(sampleToY[sample]),\
                yerr=std(sampleToY[sample]), marker='o',\
                mfc='w', ecolor='k')

    for x in xToMeans:
        y = mean(xToMeans[x])
        plot(x, y, 'ro')

    ylim(0, 1)
    yticks(fontsize=8)
    xticks(fontsize=8)
    ylabel("$\Psi$", fontsize=8)
    xlabel("Timepoint", fontsize=8)
    savefig(out_f)



# Create a table of psi values, where rows are events and columns
# are samples.
# You can specify the samples to include in the includeList file.
def psiTable(summarydir, out_f, includelist_f=False, minct=10):
 
    minct = int(minct)
    comps = os.listdir(summarydir)
    eventToSample = {}
    eventToGene = {}
    samples = {}
    for comp in comps:
        print comp
        eventToInfo, header = getSummary(os.path.join(summarydir, comp))
        for event in eventToInfo: 
            if event not in eventToSample:
                eventToSample[event] = {}

            sample1 = header[1]
            sample2 = header[4]           
            samples[sample1] = 0
            samples[sample2] = 0
    
            bf = eventToInfo[event]['bayes factor']
            s1c = eventToInfo[event]['sample 1 counts']
            s2c = eventToInfo[event]['sample 2 counts']

            ct1 = parseCounts(s1c)
            ct2 = parseCounts(s2c)

            if min([ct1, ct2]) >= minct:            
                if event not in eventToGene:
                    eventToGene[event] = [eventToInfo[event]['gene'],
                        eventToInfo[event]['symb'],\
                        eventToInfo[event]['desc']]
                if sample1 not in eventToSample[event]:
                    eventToSample[event][sample1] = \
                        eventToInfo[event]['sample 1 mean'][0]
                if sample2 not in eventToSample[event]:
                    eventToSample[event][sample2] = \
                        eventToInfo[event]['sample 2 mean'][0]
  
    if includelist_f is not False:
        samples = []
        for line in open(includelist_f):
            samples.append(line.strip())
    else:
        samples = samples.keys()
        samples.sort()
 
    out = open(out_f, 'w') 
    print len(samples), 'samples'
    out.write("#Event\t" + "\t".join(samples) + "\n")
    for event in eventToSample:
        useme = True
        for s in samples:
            if s not in eventToSample[event]:
                useme = False
                break
        if useme:
            out.write(event + "\t")
            out.write("\t".join(map(str, [eventToSample[event][s] \
                for s in samples])) + "\t")
            out.write("\t".join(eventToGene[event]) + "\n")
    out.close()


# Make a correlation matrix for visualizing correlations between
# samples using psi values.
def correlMatrix(table_f, out_f, includelist_f=False):
  
    if includelist_f is not False:
        labels = []
        for line in open(includelist_f):
            labels.append(line.strip())
        labelToIdx = {}
    
    data = []
    for line in open(table_f):
        if line.startswith("#"):
            header = line.strip().split("\t")[1:]
            if includelist_f is not False:
                for i in range(len(header)):
                    labelToIdx[header[i]] = i
            else:
                labels = header
        else:
            psivals = map(float, line.strip().split("\t")[1:-3])
            if includelist_f is not False:
                data.append([psivals[labelToIdx[labels[i]]] \
                    for i in range(len(labels))])
            else:
                data.append(psivals)

    data = array(data)
    print data.shape, 'dimensions'
 
    cc = corrcoef(data.T)
    my_cmap = jet()
    
    fig = figure(figsize=(8, 6))
    pcolor(cc, cmap=my_cmap)
    xticks(arange(data.shape[1]) + .5, labels, rotation=90, fontsize=6)
    yticks(arange(data.shape[1]) + .5, labels, fontsize=6)
    xlim(0, data.shape[1]) 
    ylim(0, data.shape[1]) 

    subplots_adjust(left=.4, bottom=.4)
    colorbar()
    savefig(out_f)


# Count the number of events regulated.
def countRegulatedFromTable(table_f, groups_f, out_f):
    
    groupToSamples = {} 
    for line in open(groups_f):
        group, sample = line.strip().split()
        try:
            groupToSamples[group].append(sample)
        except:
            groupToSamples[group] = [sample]

    data = [] 
    sampleToIdx = {}
    for line in open(table_f):
        if line.startswith("#"):
            header = line.strip().split("\t")[1:]
            for i in range(len(header)):
                sampleToIdx[header[i]] = i
        else:
            vals = line.strip().split("\t")
            psivals = map(float, vals[1:-3])
            data.append(psivals)
    data = array(data)
    print data.shape, 'events x samples'

    mindeltas = [.1, .15, .2, .3]
    groups = groupToSamples.keys()
    out = open(out_f, 'w')
    out.write("#Group1\tGroup2\tdeltapsi\tUp\tDn\tTot\n")
    for mindelta in mindeltas:
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                samples1 = groupToSamples[groups[i]]
                samples2 = groupToSamples[groups[j]]
                idx1 = [sampleToIdx[s] for s in samples1] 
                idx2 = [sampleToIdx[s] for s in samples2] 
                psi1 = data[:, array(idx1)].mean(axis=1)
                psi2 = data[:, array(idx2)].mean(axis=1)
                deltapsi = psi2 - psi1
                numup = len(where(deltapsi >= mindelta)[0])
                numdn = len(where(deltapsi <= -mindelta)[0])
                
                out.write("\t".join([groups[i], groups[j], '%.2f'%(mindelta),\
                    '%d'%(numup), '%d'%(numdn), '%d'%(data.shape[0])]) + "\n")
    out.close()


# Some parameters for plotting
params = {'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8}
rcParams.update(params)

def plotRegulated(summary_f, control, out_f):
  
    groups = ['control', '12h', '24h', '72h', '7d']
    deltaToInfo = {} 
    for line in open(summary_f):
        g1, g2, mindelta, up, dn, tot = \
            line.strip().split("\t")
        if g1 == control: 
            mindelta = float(mindelta)
            if mindelta not in deltaToInfo:
                deltaToInfo[mindelta] = {}
            deltaToInfo[mindelta][g2] = map(int, [up, dn, tot])

    figure(figsize=(8, 4))
    mindeltas = deltaToInfo.keys()
    mindeltas.sort()
    for mindelta in mindeltas:
        upvals = [0]
        dnvals = [0]
        for g in groups[1:]:
            upvals.append(deltaToInfo[mindelta][g][0]) 
            dnvals.append(deltaToInfo[mindelta][g][1]) 
            
        subplot(1, 2, 1)
        plot(arange(len(upvals)), upvals, 'ro-', alpha = .1 + 3 * mindelta,\
            label='$\Delta\Psi$ = %.2f'%(mindelta))
        
        subplot(1, 2, 2)
        plot(arange(len(dnvals)), dnvals, 'bo-', alpha = .1 + 3 * mindelta,\
            label='$\Delta\Psi$ = %.2f'%(mindelta))
                  
    subplot(1, 2, 1)
    legend(loc='upper left')
    xlim(-1, 5)
    ylim(0, 800)
    ylabel('No. of up-regulated exons')
    xticks(arange(len(groups)), groups, fontsize=10, rotation=30) 

    subplot(1, 2, 2)
    legend(loc='upper left')
    xlim(-1, 5)
    ylim(0, 800)
    ylabel('No. of down-regulated exons')
    xticks(arange(len(groups)), groups, fontsize=10, rotation=30) 

    subplots_adjust(hspace=.3, wspace=.4, bottom=.2)
    savefig(out_f)

# Count the number of events regulated.
def countRegulatedFromTables(table1_f, table2_f, groups1_f, groups2_f, out_f):
    
    group1ToSamples = {} 
    for line in open(groups1_f):
        group, sample = line.strip().split()
        try:
            group1ToSamples[group].append(sample)
        except:
            group1ToSamples[group] = [sample]

    group2ToSamples = {} 
    for line in open(groups2_f):
        group, sample = line.strip().split()
        try:
            group2ToSamples[group].append(sample)
        except:
            group2ToSamples[group] = [sample]

    data1 = [] 
    events1 = []
    sample1ToIdx = {}
    for line in open(table1_f):
        if line.startswith("#"):
            header = line.strip().split("\t")[1:]
            for i in range(len(header)):
                sample1ToIdx[header[i]] = i
        else:
            vals = line.strip().split("\t")
            psivals = map(float, vals[1:-3])
            data1.append(psivals)
            events1.append(vals[0])
    data1 = array(data1)
    print data1.shape, 'events x samples'

    data2 = [] 
    events2 = []
    sample2ToIdx = {}
    for line in open(table2_f):
        if line.startswith("#"):
            header = line.strip().split("\t")[1:]
            for i in range(len(header)):
                sample2ToIdx[header[i]] = i
        else:
            vals = line.strip().split("\t")
            psivals = map(float, vals[1:-3])
            data2.append(psivals)
            events2.append(vals[0])
    data2 = array(data2)
    print data2.shape, 'events x samples'

    mindeltas = [.1, .15, .2, .3]
    groups1 = group1ToSamples.keys()
    #groups1.sort()
    groups2 = group2ToSamples.keys()
    #groups2.sort()
    out = open(out_f, 'w')
    out.write("#Table1Group1\tTable1Group2\tTable2Group1\tTable2Group2\t"+\
        "deltapsi\tUp1Up2\tUp1Dn2\tDn1Up2\tDn1Dn2\tUp1\tDn1\tUp2\tDn2\tTable1Tot\tTable2Tot\n")
    for mindelta in mindeltas:
        for i in range(len(groups1)):
            for j in range(i + 1, len(groups1)):
                for k in range(len(groups2)):
                    for l in range(k + 1, len(groups2)):

                        t1samples1 = group1ToSamples[groups1[i]]
                        t1samples2 = group1ToSamples[groups1[j]]
                        t2samples1 = group2ToSamples[groups2[k]]
                        t2samples2 = group2ToSamples[groups2[l]]

                        t1idx1 = [sample1ToIdx[s] for s in t1samples1] 
                        t1idx2 = [sample1ToIdx[s] for s in t1samples2] 
                        t2idx1 = [sample2ToIdx[s] for s in t2samples1] 
                        t2idx2 = [sample2ToIdx[s] for s in t2samples2] 

                        t1psi1 = data1[:, array(t1idx1)].mean(axis=1)
                        t1psi2 = data1[:, array(t1idx2)].mean(axis=1)
                        t2psi1 = data2[:, array(t2idx1)].mean(axis=1)
                        t2psi2 = data2[:, array(t2idx2)].mean(axis=1)

                        t1deltapsi = t1psi2 - t1psi1
                        t2deltapsi = t2psi2 - t2psi1

                        t1up = where(t1deltapsi >= mindelta)[0]
                        t1dn = where(t1deltapsi <= -mindelta)[0]
                        t2up = where(t2deltapsi >= mindelta)[0]
                        t2dn = where(t2deltapsi <= -mindelta)[0]

                        t1upE = [events1[x] for x in t1up]
                        t1dnE = [events1[x] for x in t1dn]
                        t2upE = [events2[x] for x in t2up]
                        t2dnE = [events2[x] for x in t2dn]

                        up1up2 = len([e for e in t1upE if e in t2upE])
                        up1dn2 = len([e for e in t1upE if e in t2dnE])
                        up2dn1 = len([e for e in t1dnE if e in t2upE])
                        dn1dn2 = len([e for e in t1dnE if e in t2dnE])

                        out.write("\t".join([groups1[i], groups1[j], groups2[k], groups2[l],\
                            '%.2f'%(mindelta), '%d'%(up1up2), '%d'%(up1dn2), '%d'%(up2dn1), '%d'%(dn1dn2),\
                            '%d'%(len(t1up)), '%d'%(len(t1dn)), '%d'%(len(t2up)), '%d'%(len(t2dn)),\
                            '%d'%(data1.shape[0]), '%d'%(data2.shape[0])]) + "\n")

    out.close()



def getMeanUTRlength(eventToInfo):

    eventToLength = {} 
    for event in eventToInfo:
        utrs = event.split("@")
        starts = []
        ends = []
        for utr in utrs:
            chrom, s, e, strand = utr.split(":")
            starts.append(int(s))
            ends.append(int(e))
        if strand == '+':
            utrstart = min(starts)
            lengths = [e - utrstart for e in ends]
        else:
            utrstart = max(ends)
            lengths = [utrstart - s for s in starts]
        lengths = sorted(lengths)[::-1]
              
        m1 = list(eventToInfo[event]['sample 1 mean'])
        m2 = list(eventToInfo[event]['sample 2 mean'])
        if len(m1) == 1:
            m1.append(1 - m1[0])
            m2.append(1 - m2[0])
        m1length = sum([lengths[x] * m1[x] for x in range(len(m1))])
        m2length = sum([lengths[x] * m2[x] for x in range(len(m2))])
        eventToLength[event] = [m1length, m2length, lengths]
    return eventToLength


def compareMeanUTRlengths(summary_f, out_f):
    
    eventToInfo, header = getSummary(summary_f)
    eventToLength = getMeanUTRlength(eventToInfo)

    minbf = 5
    sigevents = [e for e in eventToInfo if eventToInfo[e]['max bf'] >= minbf]
    dpsis = [0, .2, .4, .5]
    cm = get_cmap('hot')
    for dpsi in dpsis:
        data = []
        for e in sigevents:
            m1length, m2length, lengths = eventToLength[e]
            if len(lengths) == 2 and lengths[0] - lengths[1] >= 100 and \
                abs(eventToInfo[e]['delta psi'][0]) >= dpsi:
                data.append([m1length, m2length, eventToInfo[e]['delta psi'][0]]) 

        data = array(data)
        print median(data[:, 0]), median(data[:, 1]), mean(data[:, 2])
        print data.shape
        print median((data[:, 0] - data[:, 1]))
        y, x, p = hist(data[:, 0] - data[:, 1], linspace(-500, 500, 40),\
            visible=False, normed=True, cumulative=True)
        plot(x[1:], y, color=cm(dpsi))
    savefig(out_f)
 

def fancyScatterWrapper(summarydir, outdir):
    
    for f in os.listdir(summarydir):
        fancyScatter(os.path.join(summarydir, f),\
            os.path.join(outdir, os.path.basename(f) + '.pdf'))

def fancyScatter(summary_f, out_f, mindpsi=0.05, minbf=5,\
    xtext=None, ytext=None, figheight=3, figwidth=3):

    eventToInfo, header = getSummary(summary_f) 

    events = [e for e in eventToInfo if len(eventToInfo[e]['delta psi']) == 1] 
    data = [[eventToInfo[e]['sample 1 mean'][0],\
        eventToInfo[e]['sample 2 mean'][0],\
        eventToInfo[e]['delta psi'][0],\
        eventToInfo[e]['max bf']] for e in events]
    data.sort(key=operator.itemgetter(-1), reverse=True)
    data = array(data)
    
    print data.shape[0], 'events'

    figheight = float(figheight)
    figwidth = float(figwidth)
    figure(figsize=(figwidth, figheight))
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

    mindpsi = float(mindpsi)
    minbf = float(minbf)
    fill(allrotated[:, 0], allrotated[:, 1], color='#EEEEEE', zorder=0)
    idx = where((abs(data[:, 2]) < mindpsi)|(data[:, 3] < minbf))[0]
    scatter(data[idx, 0], data[idx, 1], s=.5, color='#CCCCCC',\
        rasterized=True)

    xmin, xmax = [-.1, 1.1]
    ymin, ymax = [-.1, 1.1]

    idx = where((abs(data[:, 2]) >= mindpsi)&(data[:, 3] >= minbf))[0]
    scatter(data[idx, 0], data[idx, 1], s=3 + log(data[idx, 3]) / log(1000),\
        facecolor='w', alpha=.8, lw=.5)

    nup = len(where(data[idx, 2] < 0)[0])
    ndn = len(where(data[idx, 2] > 0)[0])
    text(-.05, 1.05, str(nup), fontsize=8, ha='left', va='top') 
    text(1.05, 0, str(ndn), fontsize=8, ha='right', va='top') 
    xlim(xmin, xmax)
    ylim(ymin, ymax) 
    if xtext is None:
        xtext = header[1].split("_")[0]
    if ytext is None:
        ytext = header[4].split("_")[0]
    xlabel('$\psi$, ' + xtext,\
        fontsize=8, horizontalalignment='center')
    ylabel('$\psi$, ' + ytext,\
        fontsize=8, horizontalalignment='center')
    title("N = %s"%(data.shape[0]), fontsize=8)

    subplots_adjust(bottom=.12)
    savefig(out_f, dpi=300)

def fancyScatterMini(summary_f, out_f, mindpsi=0.05, minbf=5,\
    xtext=None, ytext=None, figheight=3, figwidth=3):

    eventToInfo, header = getSummary(summary_f) 

    events = [e for e in eventToInfo if len(eventToInfo[e]['delta psi']) == 1] 
    data = [[eventToInfo[e]['sample 1 mean'][0],\
        eventToInfo[e]['sample 2 mean'][0],\
        eventToInfo[e]['delta psi'][0],\
        eventToInfo[e]['max bf']] for e in events]
    data.sort(key=operator.itemgetter(-1), reverse=True)
    data = array(data)
    
    print data.shape[0], 'events'

    figheight = float(figheight)
    figwidth = float(figwidth)
    figure(figsize=(figwidth, figheight))
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

    mindpsi = float(mindpsi)
    minbf = float(minbf)
    fill(allrotated[:, 0], allrotated[:, 1], color='#EEEEEE', zorder=0)
    idx = where((abs(data[:, 2]) < mindpsi)|(data[:, 3] < minbf))[0]
    scatter(data[idx, 0], data[idx, 1], s=.1, color='#CCCCCC',\
        rasterized=True)

    xmin, xmax = [-.1, 1.1]
    ymin, ymax = [-.1, 1.1]

    idx = where((abs(data[:, 2]) >= mindpsi)&(data[:, 3] >= minbf))[0]
    scatter(data[idx, 0], data[idx, 1], s=1 + log(data[idx, 3]) / log(100000),\
        facecolor='w', alpha=.8, lw=.25)

    nup = len(where(data[idx, 2] < 0)[0])
    ndn = len(where(data[idx, 2] > 0)[0])
    text(0, 1, str(nup), fontsize=6, ha='left', va='top') 
    text(1, 0, str(ndn), fontsize=6, ha='right', va='bottom') 
    xlim(xmin, xmax)
    ylim(ymin, ymax) 
    if xtext is None:
        xtext = header[1].split("_")[0]
    if ytext is None:
        ytext = header[4].split("_")[0]
    xlabel('$\psi$, ' + xtext,\
        fontsize=8, horizontalalignment='center')
    ylabel('$\psi$, ' + ytext,\
        fontsize=8, horizontalalignment='center')
    title("N = %s"%(data.shape[0]), fontsize=8)
    xticks([0, .25, .5, .75, 1], [0, .25, .5, .75, 1], fontsize=6)
    yticks([0, .25, .5, .75, 1], [0, .25, .5, .75, 1], fontsize=6)

    subplots_adjust(bottom=.2, left=.2)
    savefig(out_f, dpi=300)

def fancyScatterOld(summary_f, out_f):

    eventToInfo, header = getSummary(summary_f) 

    events = [e for e in eventToInfo if len(eventToInfo[e]['delta psi']) == 1] 
    data = array([[eventToInfo[e]['sample 1 mean'][0], eventToInfo[e]['sample 2 mean'][0],\
        eventToInfo[e]['delta psi'][0], eventToInfo[e]['max bf']] for e in events])
    
    print data.shape[0], 'events'

    figure(figsize=(4, 4))
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
 
    dpsis = [.25, .5, 1]
    minbf = 100 
    fill(allrotated[:, 0], allrotated[:, 1], color='#EEEEEE', zorder=0)
    idx = where((data[:, 2] <= dpsis[0])&(data[:, 2] >= -dpsis[0]))[0]
    scatter(data[idx, 0], data[idx, 1], s=.5, color='#CCCCCC')

    xmin, xmax = [-.25, 1.25]
    ymin, ymax = [-.25, 1.25]

    for i in range(len(dpsis) - 1):
        dpsi = dpsis[i]
        plot([xmin, xmax], [xmin - dpsi, xmax - dpsi], 'b--', alpha=dpsi)
        plot([xmin, xmax], [xmin + dpsi, xmax + dpsi], 'r--', alpha=dpsi)
        idx = where((data[:, 2] >= dpsi)&(data[:, 2] < dpsis[i + 1])&\
            (data[:, 3] >= minbf))[0]
        print len(idx)
        scatter(data[idx, 0], data[idx, 1], s=20, color='b', alpha=dpsi, lw=0)
        idx = where((data[:, 2] <= -dpsi)&(data[:, 2] > -dpsis[i + 1])&\
            (data[:, 3] >= minbf))[0]
        print len(idx)
        scatter(data[idx, 0], data[idx, 1], s=20, color='r', alpha=dpsi, lw=0)

    xlim(xmin, xmax)
    ylim(ymin, ymax) 
    xlabel(header[1].split("_")[0] + ' $\psi$ ',\
        fontsize=12, horizontalalignment='center')
    ylabel(header[4].split("_")[0] + ' $\psi$ ',\
        fontsize=12, horizontalalignment='center')
 
    #dpsis = [.25, .5]
    #colors = ['#CCCCCC', 'r']
    #minbf = 5

    #for i in range(len(dpsis)):
    #    dpsi = dpsis[i]
    #    idx = where((abs(data[:, 2]) >= dpsi)&(data[:, 3] >= minbf))[0]
    #    scatter(data[idx, 0], data[idx, 1], color=colors[i],\
    #        s=(abs(data[idx, 2]) * 10) ** 2, alpha=1)

    savefig(out_f)


def plotHeatmap(event, summarydir, names, out_f):
 
    files = os.listdir(summarydir)
    names = names.split(",")
    
    sampleToPsi = {}
    pairToBF = {} 
    for f in files:
        cmd = 'grep ' + event + ' ' + os.path.join(summarydir, f)
        pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        sample1, sample2 = f.split(".")[0].split("_vs_")
        for line in pid.stdout:
            vals = line.strip().split("\t")
            psi1 = float(vals[1])
            psi2 = float(vals[4])
            bf = float(vals[8])
            if bf > 1000:
                bf = 1000
            sampleToPsi[sample1] = psi1
            sampleToPsi[sample2] = psi2
            if sample1 not in pairToBF:
                pairToBF[sample1] = {}
            if sample2 not in pairToBF:
                pairToBF[sample2] = {}
            pairToBF[sample1][sample2] = bf
            pairToBF[sample2][sample1] = bf

    deltas = zeros((len(names), len(names)), dtype='f') 
    bfs = zeros((len(names), len(names)), dtype='f') 
    for i in range(len(names)):
        for j in range(len(names)):
            if i != j:
                deltapsi = sampleToPsi[names[j]] - sampleToPsi[names[i]]
                bf = pairToBF[names[i]][names[j]]
                deltas[j, i] = deltapsi
                bfs[j, i] = bf
                print names[i], names[j], deltapsi, bf

    figure(figsize=(1, .8))
    sigdata = ma.masked_where(bfs < 5, deltas)
    nsigdata = ma.masked_where(bfs > 5, deltas) 
    pcolor(sigdata[::-1, :], cmap=get_cmap('RdBu_r'), vmin=-.5, vmax=.5)
    for i in range(len(names)):
        for j in range(len(names)):
            if i != j:
                if bfs[j, i] > 100:
                    t = '***'
                elif bfs[j, i] > 20:
                    t = '**'
                elif bfs[j, i] > 5:
                    t = '*'
                else:
                    t = ''
                text(i + .5,\
                    bfs.shape[1] - j - .5, t,\
                    fontsize=6, ha='center', va='center')
                    
    xticks(arange(len(names)) + .5, names, fontsize=8, rotation=90)
    yticks(arange(len(names)) + .5, names[::-1], fontsize=8)
    subplots_adjust(left=.2, bottom=.2)
    colorbar()
    savefig(out_f) 

def plotHeatmap2(event, summarydir, names1, names2, out_f):
 
    files = os.listdir(summarydir)
    names1 = names1.split(",")
    names2 = names2.split(",")
    
    sampleToPsi = {}
    pairToBF = {} 
    for f in files:
        cmd = 'grep ' + event + ' ' + os.path.join(summarydir, f)
        pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        sample1, sample2 = f.split(".")[0].split("_vs_")
        for line in pid.stdout:
            vals = line.strip().split("\t")
            psi1 = float(vals[1])
            psi2 = float(vals[4])
            bf = float(vals[8])
            if bf > 1000:
                bf = 1000
            sampleToPsi[sample1] = psi1
            sampleToPsi[sample2] = psi2
            if sample1 not in pairToBF:
                pairToBF[sample1] = {}
            if sample2 not in pairToBF:
                pairToBF[sample2] = {}
            pairToBF[sample1][sample2] = bf
            pairToBF[sample2][sample1] = bf

    deltas = zeros((len(names2), len(names1)), dtype='f') 
    bfs = zeros((len(names2), len(names1)), dtype='f') 
    for i in range(len(names2)):
        for j in range(len(names1)):
            deltapsi = sampleToPsi[names1[j]] - sampleToPsi[names2[i]]
            bf = pairToBF[names2[i]][names1[j]]
            deltas[j, i] = deltapsi
            bfs[j, i] = bf
            print names2[i], names1[j], deltapsi, bf

    figure(figsize=(1, .8))
    sigdata = ma.masked_where(bfs < 5, deltas)
    nsigdata = ma.masked_where(bfs > 5, deltas) 
    pcolor(sigdata[::-1, :], cmap=get_cmap('RdBu_r'), vmin=-.5, vmax=.5)
    for i in range(len(names1)):
        for j in range(len(names2)):
            if bfs[j, i] > 100:
                t = '***'
            elif bfs[j, i] > 20:
                t = '**'
            elif bfs[j, i] > 5:
                t = '*'
            else:
                t = ''
            text(i + .5,\
                bfs.shape[1] - j - .5, t,\
                fontsize=6, ha='center', va='center')
    xticks(arange(len(names2)) + .5, names2, fontsize=8, rotation=90)
    yticks(arange(len(names1)) + .5, names1[::-1], fontsize=8)
    subplots_adjust(left=.2, bottom=.2)
    colorbar()
    savefig(out_f) 

def plotBars(event, summarydir, names1, names2, out_f):
 
    files = os.listdir(summarydir)
    names1 = names1.split(",")
    names2 = names2.split(",")
    
    sampleToPsi = {}
    pairToBF = {} 
    for f in files:
        cmd = 'grep ' + event + ' ' + os.path.join(summarydir, f)
        pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        sample1, sample2 = f.split(".")[0].split("_vs_")
        for line in pid.stdout:
            vals = line.strip().split("\t")
            psi1 = float(vals[1])
            lpsi1 = float(vals[2])
            hpsi1 = float(vals[3])
            psi2 = float(vals[4])
            lpsi2 = float(vals[5])
            hpsi2 = float(vals[6])
            bf = float(vals[8])
            if bf > 1000:
                bf = 1000
            sampleToPsi[sample1] = [lpsi1, psi1, hpsi1]
            sampleToPsi[sample2] = [lpsi2, psi2, hpsi2]
            if sample1 not in pairToBF:
                pairToBF[sample1] = {}
            if sample2 not in pairToBF:
                pairToBF[sample2] = {}
            pairToBF[sample1][sample2] = bf
            pairToBF[sample2][sample1] = bf
            print sample1, sample2, bf


    figure(figsize=(3, 2))
    data1 = []
    data2 = []
    for i in range(len(names1)):
        data1.append(sampleToPsi[names1[i]]) 
    for i in range(len(names2)):
        data2.append(sampleToPsi[names2[i]]) 
    data1 = array(data1)
    data2 = array(data2)

    colors = ['#2D983A', '#9ED29A']

    bar(arange(data1.shape[0]) + .15, data1[:, 1], width=.3, color=colors[0],\
        yerr=[data1[:, 1] - data1[:, 0], data1[:, 2] - data1[:, 1]],\
        ecolor='k')
    bar(arange(data2.shape[0]) + .5, data2[:, 1], width=.3, color=colors[1],\
        yerr=[data2[:, 1] - data2[:, 0], data2[:, 2] - data2[:, 1]],\
        ecolor='k')
    yticks(fontsize=8)
    ylabel("Distal pA usage ($\Psi$)", fontsize=8)
    xticks(arange(data1.shape[0]) + .5, names1, fontsize=8, rotation=60,\
        ha='right') 

    subplots_adjust(left=.2, bottom=.2)
    savefig(out_f) 


def filterForSS(summarydir, eventtype, ss_f, outdir):
 
    eventToInfo = {} 
    for line in open(ss_f):
        if line.startswith("#"):
            header = line.strip().split("\t")
        else:
            vals = line.strip().split("\t")
            event = vals[0]
            seqs = [x.split(":")[0] for x in vals[1:]]
            strengths = map(float, [x.split(":")[1] for x in vals[1:]])

            remove = False
            if eventtype == 'SE' or eventtype == 'SEexon' or eventtype=='Trio':
                if seqs[1][3:5] != 'GT' or seqs[2][-5:-3] != 'AG' or \
                    seqs[3][3:5] != 'GT' or seqs[4][-5:-3] != 'AG':
                    remove = True
            elif eventtype == 'A3SS':
                if seqs[1][3:5] != 'GT' or seqs[2][-5:-3] != 'AG' or \
                    seqs[3][-5:-3] != 'AG':
                    remove = True
            elif eventtype == 'A5SS':
                if seqs[1][3:5] != 'GT' or seqs[2][3:5] != 'GT' or \
                    seqs[3][-5:-3] != 'AG':
                    remove = True
            elif eventtype == 'MXE':
                if seqs[1][3:5] != 'GT' or seqs[2][-5:-3] != 'AG' or \
                    seqs[3][3:5] != 'GT' or seqs[4][-5:-3] != 'AG' or \
                    seqs[5][3:5] != 'GT' or seqs[6][-5:-3] != 'AG':
                    remove = True
            if remove == False:
                eventToInfo[event] = [seqs, strengths]

    for f in os.listdir(summarydir):
        print f
        out_f = os.path.join(outdir, f)
        out = open(out_f, 'w') 
        for line in open(os.path.join(summarydir, f)):
            if not line.startswith("#"):
                event = line.strip().split("\t")[0]
                if event in eventToInfo:
                    out.write(line)
            else:
                out.write(line)
        out.close()

# Make a table.
def makeTable(settings_f, out_f, requireAll=False):

    if isinstance(requireAll, str):
        requireAll = eval(requireAll)

    files = []
    etypes = []
    nameorder = []
    allcomparisons = []
    comparisons = []
    uniquecomparisons = []
    for line in open(settings_f):
        f, etype, name1, name2 = line.strip().split()
        if name1 not in nameorder:
            nameorder.append(name1)
        if name2 not in nameorder:
            nameorder.append(name2)
        files.append(f)
        etypes.append(etype)
        comparison = "_vs_".join([name1, name2])
        comparisons.append(comparison)
        if comparison not in uniquecomparisons:
            uniquecomparisons.append(comparison)
   
    allinfo = {} 
    headers = []
    eventToType = {}
    for i in range(len(files)):
        eventToInfo, header = getSummary(files[i])
        headers.append(header)
        comparison = comparisons[i]
        for e in eventToInfo:
            if e not in allinfo:
                allinfo[e] = {}
                eventToType[e] = etypes[i]
            allinfo[e][comparison] = eventToInfo[e]
 
    out = open(out_f, 'w') 
    line = ["#Event", "Event Type"]
    for name in nameorder:
        line.extend([\
            name + "_mean",\
            name + "_low",\
            name + "_high"])
    for comparison in uniquecomparisons:
        name1, name2 = comparison.split("_vs_")
        line.extend([\
            name2 + "_minus_" + name1,\
            name2 + "_vs_" + name1 + "_BF"])
    line.extend(["MaxAbsDelta", "MaxBF", "Gene", "Symb", "Desc"])
    out.write("\t".join(line) + "\n")
    data = []
    for e in allinfo: 
        etype = eventToType[e]
        line = [e, etype]
        bfs = []
        dpsis = []
        skipme = False 
        nameToInfo = {}
        for comparison in uniquecomparisons:
            name1, name2 = comparison.split("_vs_")
            if comparison not in allinfo[e]:
                if requireAll:
                    skipme = True
                    break
                nameToInfo[name1] = ["NA" for i in range(3)]
                nameToInfo[name2] = ["NA" for i in range(3)]
            else:
                if name1 not in nameToInfo:
                    nameToInfo[name1] = [\
                        allinfo[e][comparison]['sample 1 mean'][0],\
                        allinfo[e][comparison]['sample 1 low'][0],\
                        allinfo[e][comparison]['sample 1 high'][0]]
                if name2 not in nameToInfo:
                    nameToInfo[name2] = [\
                        allinfo[e][comparison]['sample 2 mean'][0],\
                        allinfo[e][comparison]['sample 2 low'][0],\
                        allinfo[e][comparison]['sample 2 high'][0]]
                nameToInfo[comparison] = [\
                    allinfo[e][comparison]['delta psi'][0],\
                    allinfo[e][comparison]['bayes factor'][0]]
                dpsis.append(abs(allinfo[e][comparison]['delta psi'][0]))
                bfs.append(allinfo[e][comparison]['bayes factor'][0])
                gene = allinfo[e][comparison]['gene']
                symb = allinfo[e][comparison]['symb']
                desc = allinfo[e][comparison]['desc']

        if skipme == False:
            maxbf = max(bfs)
            maxdpsi = max(dpsis)
            for name in nameorder:
                line.extend(nameToInfo[name])
            for comparison in uniquecomparisons:
                line.extend(nameToInfo[comparison])
            line.extend([maxdpsi, maxbf, gene, symb, desc])
            data.append(line)

    data.sort(key=operator.itemgetter(-4), reverse=True)

    for line in data:
        out.write("\t".join(map(str, line)) + "\n")
    out.close()


# Read psi table and return eventToInfo dictionary of 
# eventToInfo[event][header name] -> value
def readTable(table_f):
   
    eventToInfo = {} 
    for line in open(table_f):
        if line.startswith("#"):
            header = line.strip().split("\t")
        else:
            vals = line.strip().split("\t")    
            e = vals[0]
            eventToInfo[e] = {}
            for i in range(1, len(header)):
                if header[i] in ['Event Type', 'Gene', 'Symb', 'Desc']:
                    eventToInfo[e][header[i]] = vals[i]
                else: 
                    eventToInfo[e][header[i]] = float(vals[i])
    return eventToInfo


# Generate plots from consolidated data.
# Consolidated data comprises a column for low, mean, and high values,
# BF values for each pairwise comparison, and gene/symb/desc info.
def plotFromConsolidated(consolidated_f, groups_f, event, minbf, out_f=False):
  
    from matplotlib import patches 
    minbf = float(minbf)

    # Get group information 
    groupToSamples = {} 
    groups = []
    samples = []
    for line in open(groups_f):
        sample, group = line.strip().split()
        group = int(group)
        groups.append(group)
        samples.append(sample)
        try:
            groupToSamples[group].append(sample)
        except:
            groupToSamples[group] = [sample]
    groups = list(set(groups))
    
    # Get information for event
    for line in open(consolidated_f):
        if line.startswith("#"):
            header = line.strip().split("\t")
        else:
            vals = line.strip().split("\t")
            if vals[0] == event:
                break
    
    sampleToPsi = {} 
    sampleToBF = {}
    for i in range(len(header)):
        if header[i].endswith("_low") or header[i].endswith("_mean") or \
            header[i].endswith("_high"):
            if header[i].endswith("_low"):
                sample = header[i][:-4]
                idx = 0
            elif header[i].endswith("_mean"):
                sample = header[i][:-5]
                idx = 1
            elif header[i].endswith("_high"):
                sample = header[i][:-5]
                idx = 2
            if sample not in sampleToPsi:
                sampleToPsi[sample] = [0, 0, 0]
            try:
                sampleToPsi[sample][idx] = float(vals[i])
            except:
                sampleToPsi[sample] = [0, -1, 0]
        elif "_vs_" in header[i]:
            s1, s2 = header[i].split("_vs_")    
            if s1 not in sampleToBF:
                sampleToBF[s1] = {} 
            sampleToBF[s1][s2] = float(vals[i])
            if s2 not in sampleToBF:
                sampleToBF[s2] = {} 
            sampleToBF[s2][s1] = float(vals[i])
        elif header[i] == 'gene':
            gene = vals[i]
        elif header[i] == 'symb':
            symb = vals[i]
        elif header[i] == 'desc':
            desc = vals[i]

    xlabels = []
    groups.sort()
    n = 0
    m = 0
    fig = figure(figsize=(6, 6))
    subplot(212)
    sampleToX = {}
    for group in groups:
        subsamples = groupToSamples[group]
        xvals = []
        yvals = []
        for sample in subsamples:
            xlabels.append(sample)
            yvals.append(sampleToPsi[sample]) 
            xvals.append(m)
            sampleToX[sample] = m
            m += 1
        yvals = array(yvals)
        errorbar(array(xvals) + .5, yvals[:, 1], yerr=[yvals[:, 1] - yvals[:, 0],\
            yvals[:, 2] - yvals[:, 1]], fmt='o', color=cm.jet(n / float(len(groups)))) 
        n += 1
    xticks(arange(m) + .5, xlabels, rotation=90, fontsize=8)
    ylabel("$\Psi$") 
    xlim(0, m)
    ylim(0, 1)

    ax = subplot(211, frameon=False)
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            bf = sampleToBF[samples[i]][samples[j]] 
            if bf >= minbf:
                x1 = sampleToX[samples[i]] + .5
                x2 = sampleToX[samples[j]] + .5
                midpt = (x1 + x2) / 2.
                arc = patches.Arc(xy=(midpt, 0.0), width=x2 - x1,\
                    height=(x2 - x1) / len(samples),\
                    angle=0.0, theta1=0.0, theta2=180.0,\
                    color=str(1 - (min(log10(bf) / 10, 1.0))))
                ax.add_patch(arc)

    xlim(0, m)
    ylim(0, 1)
    xticks([])
    yticks([])

    suptitle("Event: %s\nGene: %s\nSymbol: %s\n"%(event, gene, symb),\
        fontsize=8, multialignment='center')
    subplots_adjust(bottom=.5) 
    if not out_f:
        savefig(out_f) 
     

# Collapse CIs (get a weighted average of psi values with sample groups)
def collapseCI(consolidated_f, groups_f, out_f):

    groupToSamples = {}
    sampleToGroup = {}
    groups = []
    for line in open(groups_f):
        sample, g = line.strip().split()
        groups.append(g)
        try:
            groupToSamples[g].append(sample) 
        except:
            groupToSamples[g] = [sample]
        sampleToGroup[sample] = g
    groups = groupToSamples.keys()
    print len(sampleToGroup), 'samples'

    out = open(out_f, 'w')
    out.write("#Event")
    for g in groups:
        out.write("\t" + "\t".join([g + "_low", g + "_mean", g + "_high"])) 
    out.write("\tGene\tSymb\tDesc\n")
    for line in open(consolidated_f):
        if line.startswith("#"):
            header = line.strip().split("\t")
            for i in range(1, len(header)):
                if header[i].endswith("_low"):
                    idx = 0
                    sample = header[i][:-4]
                elif header[i].endswith("_mean"):
                    idx = 1
                    sample = header[i][:-5]
                elif header[i].endswith("_high"):
                    idx = 2 
                    sample = header[i][:-5]
                if "_vs_" in header[i]:
                    break
        else:
            vals = line.strip().split("\t")
            event = vals[0]
            out.write(event)
            groupToVals = {} 
            for i in range(len(sampleToGroup)):
                low = vals[3 * i + 1] 
                mean = vals[3 * i + 2] 
                high = vals[3 * i + 3]
                sample = header[3 * i + 1][:-4]
                group = sampleToGroup[sample]
                if 'n/a' not in [low, mean, high]:
                    try:
                        groupToVals[group].append(map(float, [low, mean, high])) 
                    except:
                        groupToVals[group] = [map(float, [low, mean, high])]
            for g in groups:
                if g in groupToVals:
                    mat = array(groupToVals[g])
                    CIs = mat[:, 2] - mat[:, 0]
                    newmean = (mat[:, 1] * 1 / CIs).sum()  / (1 / CIs).sum() 
                    newCI = (CIs * 1 / CIs).sum()  / (1 / CIs).sum() 
                    #newCI = 1 / (1 / (CIs * CIs)).sum()
                    print mat[:, 1], CIs, newmean, newCI
                    newlow = max([0, round(newmean - newCI / 2, 2)])
                    newhigh = min([1, round(newmean + newCI / 2, 2)])
                    newmean = round(newmean, 2)
                    newmean2 = (mat[:, 1] * 1 / (CIs * CIs)).sum()  / (1 / (CIs * CIs)).sum() 
                    newmean2 = round(newmean2, 2)
                    out.write("\t" + "\t".join(map(str, [newlow, newmean, newhigh])))
                else:
                    out.write("\tn/a\tn/a\tn/a")
            out.write("\t" + "\t".join(vals[-3:]) + "\n")
    out.close()
           





 
