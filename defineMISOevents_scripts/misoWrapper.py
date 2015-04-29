import os, sys, shelve, operator, subprocess, time
import run_events_analysis
import samples_utils
from misopy.settings import Settings, load_settings
import Mypbm, Mysge
from pylab import *
import random


# Function to run miso on all samples in a bam/ directory.
# PEdist_f is 'False' if single end.
# For each bam file, copy to a node and run in non-cluster mode.
def runMISO(bamDir, pickledDir, overhanglen, PEdist_f, outDir,\
    settings_f, queue, ppn, clustertype):

    import pysam
    bamDir = os.path.abspath(os.path.expanduser(bamDir))
    pickledDir = os.path.abspath(os.path.expanduser(pickledDir))
    outDir = os.path.abspath(os.path.expanduser(outDir))

    overhanglen = int(overhanglen)
    ppn = int(ppn)
    fileToPE = {}
    if PEdist_f != 'False':
        for line in open(PEdist_f):
            if not line.startswith("#"):
                file, mean, sd = line.strip().split()[:3]
                fileToPE[file.split(".")[0]] = [mean,sd]   

    files = [f for f in os.listdir(bamDir) if f.endswith(".bam")]
    for f in files:
        fname = f.split(".")[0]
        if fname in fileToPE:
            PEinfo = map(float, fileToPE[fname])
        else:
            PEinfo = False 

        # Get read length
        bam = pysam.Samfile(os.path.join(bamDir, f), 'rb')
        readlen = False 
        for read in bam.fetch():
            readlen = int(read.qlen)
            break 
        if readlen == False:
            exit('Read length could not be obtained!')
 
        outdir = os.path.abspath(os.path.join(outDir, fname))
        if not os.path.exists(outdir):
            if queue != 'False':
                scriptOptions = {'ppn':ppn, 'jobname':f}
                if clustertype == 'torque':
                    cmd = 'python %s/script.py misoWrapper.runMISOsingle %s %s %s %s %s %s %s %s'\
                        %(sys.path[0], os.path.abspath(pickledDir),\
                        os.path.abspath(os.path.join(bamDir,f)), readlen,\
                        overhanglen,\
                        outdir,\
                        PEinfo,\
                        settings_f,\
                        '/scratch/et_wang/')
                    print 'Running', cmd
                    Mypbm.launchJob(cmd, scriptOptions, queue_type=queue, verbose=True)
                elif clustertype == 'sge':
                    cmd = 'python %s/script.py misoWrapper.runMISOsingle %s %s %s %s %s %s %s %s'\
                        %(sys.path[0], os.path.abspath(pickledDir),\
                        os.path.abspath(os.path.join(bamDir,f)), readlen,\
                        overhanglen,\
                        outdir,\
                        PEinfo,\
                        settings_f,\
                        '/data/')
                    print 'Running', cmd
                    Mysge.launchJob(cmd, scriptOptions, queue_type='all', verbose=True)
            else:
                runMISOlocal(os.path.abspath(pickledDir),\
                    os.path.abspath(os.path.join(bamDir, f)), readlen,\
                    overhanglen,\
                    os.path.abspath(os.path.join(outDir, fname)),\
                    PEinfo,\
                    settings_f)


# Run a single bam file with miso.
def runMISOsingle(pickledDir, bamFile, readlen, overhanglen, outdir,\
    paired_end, settings_f, scratchDir):

    if paired_end == 'False':
        paired_end = None

    t = str(time.time()) + str(random.random())

    print os.path.basename(pickledDir)
    if not os.path.exists(scratchDir):
        cmd = 'mkdir ' + scratchDir
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
 
    # Copy pickled dir.
    pickled = os.path.join(scratchDir, os.path.basename(pickledDir) + \
        "." + t)
    cmd = 'mkdir ' + pickled
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    cmd = 'cp -r ' + pickledDir + '/* ' + pickled
    process = subprocess.Popen(cmd, shell=True)
    process.wait()

    # Copy bam file. 
    cmd = 'cp -fL ' + bamFile + ' ' + scratchDir
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    cmd = 'cp -fL ' + bamFile + '.bai ' + scratchDir
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    bam = os.path.join(scratchDir, os.path.basename(bamFile))

    # Give output directory in scratch a timestamp
    out = os.path.join(scratchDir, os.path.basename(outdir + "." + t))

    # LOAD SETTINGS FOR MISO
    Settings.load(settings_f)
 
    run_events_analysis.compute_all_genes_psi(\
        pickled, bam, int(readlen), out, overhang_len=int(overhanglen),\
        paired_end=paired_end, settings_fname=settings_f, prefilter=True)

    # Summarize sample
    #summary_fname = os.path.join(out, os.path.basename(outdir) + '.miso_summary') 
    #samples_utils.summarize_sampler_results(out, summary_fname)

    if not os.path.exists(outdir):
        cmd = 'mkdir -p ' + outdir
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        
    # Copy output back.
    cmd = 'cp -r ' + out + '/* ' + outdir
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
 
    # Remove bam, output, and pickled dir. 
    cmd = 'rm ' + bam
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    cmd = 'rm ' + bam + '.bai'
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    cmd = 'rm -fr ' + out
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    cmd = 'rm -fr ' + pickled
    process = subprocess.Popen(cmd, shell=True)
    process.wait()

# Run a single bam file with miso. Run locally, i.e. do not copy to a node.
def runMISOlocal(pickledDir, bamFile, readlen, overhanglen, outdir,\
    paired_end, settings_f):

    if paired_end == False or paired_end == 'False':
        paired_end = None
    Settings.load(settings_f)

    #if not os.path.exists(outdir):
    #    print 'running', outdir
    run_events_analysis.compute_all_genes_psi(\
        pickledDir, bamFile, int(readlen), outdir,\
        overhang_len=int(overhanglen),\
        paired_end=paired_end, settings_fname=settings_f)


# Function to compare miso distributions between samples.
def compareMISO(eventDir, outDir, samples, usecluster=False):

    eventDir = os.path.abspath(os.path.expanduser(eventDir))
    outDir = os.path.abspath(os.path.expanduser(outDir))

    if isinstance(usecluster, str):
        usecluster = eval(usecluster)   

    if not os.path.exists(outDir):
        os.popen("mkdir "+outDir)

    if samples == 'all':
        samples = sorted(os.listdir(eventDir))
    else:
        samples = samples.split(",") 
    for i in range(len(samples)):
        for j in range(i+1,len(samples)):
            cmd = 'python /home/et_wang/Tools/pythonmodules/lib/python2.7/site-packages/misopy/run_miso.py '+\
                '--compare-samples '+os.path.join(eventDir,samples[i])+' '+\
                os.path.join(eventDir,samples[j])+\
                ' '+outDir

            if usecluster:
                cmd = 'echo "'+cmd+'"' + "| qsub -q short " 
            print cmd
            os.popen(cmd)

# Compute the meta-psi for each group of samples.
# Only compute for events where there is data in every sample
# in the group.
def compute_meta_psi(eventDir, groups_f, outDir):
  
    import samples_utils 
    import miso_meta

    groupToSamples = {} 
    for line in open(groups_f):
        sample, group = line.strip().split()
        try:
            groupToSamples[group].append(sample)
        except: 
            groupToSamples[group] = [sample]

    if not os.path.isdir(outDir):
        os.mkdir(outDir)

    for group in groupToSamples:
        print "Group:", group
        samples = groupToSamples[group]
        sampleToFiles = {}
        for s in samples:
            print "Getting event files:", s
            sampleToFiles[s] = samples_utils.get_samples_dir_filenames(\
                os.path.abspath(os.path.join(eventDir, s)))
        files = sampleToFiles[samples[0]] 
        for f in files:     # iterate through each miso file
            eventname = samples_utils.get_event_name(f)
            psivals = []
            for s in samples:
                fname = filter(lambda filename: samples_utils.get_event_name(filename) ==\
                    eventname, sampleToFiles[s])
                if len(fname) == 0:
                    break
                else: 
                    fname = fname[0] 
                    results = samples_utils.load_samples(fname)
                    posterior = results[0][:, 0]
                    header = results[1]
                    isoforms = samples_utils.get_isoforms_from_header(header[0])
                    counts_info = results[5]
                    psivals.append(posterior) 
                    #print header, isoforms, counts_info 

            if len(psivals) == len(samples):
                means = [x.mean() for x in psivals]
                CIs = [x.var() for x in psivals]
                weights = [1 / x for x in CIs]
                newpsi = miso_meta.addDistributions(psivals, weights)
                print means, CIs, newpsi.mean(), newpsi.var()




# Summarize miso output using summarize-samples
# indir is a directory that contains all the miso directories
def summarize(indir, outdir):
   
    if not os.path.exists(outdir):
        os.mkdir(os.path.abspath(outdir)) 
    misodirs = os.listdir(indir)
    for d in misodirs:
        samples_dir = os.path.abspath(os.path.join(indir, d))
        summary_fname = os.path.abspath(os.path.join(outdir, d + '.miso_summary')) 
        cmd = 'python /home/et_wang/Tools/pythonmodules/lib/python2.7/site-packages/misopy/' +\
            'script.py samples_utils.summarize_sampler_results %s %s'%(samples_dir, summary_fname)
        print 'Running', cmd
        scriptOptions = {'ppn':1, 'jobname':'summarize_miso' + d}
        Mypbm.launchJob(cmd, scriptOptions, queue_type='short', verbose=True)


# Function to summarize a comparison between two samples.
# Outputs a text file with the following fields:
# 1. event name
# 2. sample 1 mean
# 3. sample 1 low
# 4. sample 1 high
# 5. sample 2 mean
# 6. sample 2 low
# 7. sample 2 high
# 8. delta psi
# 9. bayes factor
# 10. isoforms
# 11. sample 1 counts
# 12. sample 1 assigned counts
# 13. sample 2 counts
# 14. sample 2 assigned counts
# 15. sample 2 assigned counts
# 16. max bayes factor (multiple isoform max)
# 17. gene
# 18. symbol
# 19. description
def summarizeVs(indir, ensGeneMap_f, outdir):

    eventToInfo = shelve.open(ensGeneMap_f, 'r') 
 
    dirs = [f for f in os.listdir(indir) if "_vs_" in f]
    for dir in dirs:
        print dir
        fname = os.path.basename(dir)
        sample1, sample2 = fname.split("_vs_")
        bf_f = os.path.join(indir,dir,'bayes-factors',fname+".miso_bf")

        data = []
        for line in open(bf_f):
            if not line.startswith("event_name"):
                vals = line.strip().split("\t")
                bf = vals[8]
                if ',' in bf:
                    maxbf = max(map(float,bf.split(",")))
                else:
                    maxbf = float(bf)
                vals.append(maxbf)
                data.append(vals)
        data.sort(key=operator.itemgetter(-1),reverse=True)
        
        out = open(os.path.join(outdir,fname+".miso_bf"),'w')
        out.write("\t".join(["#event_name",\
            sample1+"_posterior_mean",sample1+"_ci_low",sample1+"_ci_high",\
            sample2+"_posterior_mean",sample2+"_ci_low",sample2+"_ci_high",\
            "diff","bayes_factor","isoforms",\
            sample1+"_counts",sample1+"_assigned_counts",\
            sample2+"_counts",sample2+"_assigned_counts", "chrom", "strand",\
            "mRNA_starts", "mRNA_ends", "max_bf"])+"\n")

        for item in data:
            event = item[0]
            try:
                item.extend(eventToInfo[item[0]])
            except: 
                item.extend(["n/a","n/a","n/a"])
            out.write("\t".join(map(str,item))+"\n")
        out.close() 



def summarizeCts(indir,outdir):

    dirs = [f for f in os.listdir(indir) if ".counts" in f]
    for dir in dirs:
        print dir
        chrdirs = [os.path.join(indir,dir,f) for f in os.listdir(dir) if f.startswith("chr")]
        data = []
        for chr in chrdirs:
            print chr
            events = [os.path.join(chr,f) for f in os.listdir(chr)]
            for event in events:
                name = os.path.basename(event).split(".")[0]
                num = -1
                for line in open(event):
                    vals = line.strip().split("\t")
                    for v in vals:
                        if 'num_reads' in v:
                            label,num = v.split("=")
                            break
                if num>-1: 
                    data.append([name,num])
        out = open(os.path.join(outdir,dir+".summary"),'w')
        out.write("#Event\tnum_reads\n")
        out.write("\n".join(["\t".join(x) for x in data]))
        out.close()

# make symbolic links for bam files (from rnaseqlib pipeline)
def symBAM(mappingdir, bamdir):
    
    dirs = os.listdir(mappingdir)
    for d in dirs:
        oldfile1 = os.path.abspath(os.path.join(mappingdir, d, "accepted_hits.sorted.bam"))
        oldfile2 = os.path.abspath(os.path.join(mappingdir, d, "accepted_hits.sorted.bam.bai"))
        newfile1 = os.path.abspath(os.path.join(bamdir, d + ".bam"))
        newfile2 = os.path.abspath(os.path.join(bamdir, d + ".bam.bai"))
        cmd1 = 'cp -s ' + oldfile1 + ' ' + newfile1
        cmd2 = 'cp -s ' + oldfile2 + ' ' + newfile2
        subprocess.Popen(cmd1, shell=True)
        subprocess.Popen(cmd2, shell=True)


def moveDirs(indir, prefix):
    
    dirs = os.listdir(indir)
    for d in dirs:
        f = os.listdir(os.path.abspath(os.path.join(indir, d)))[0]
        cmd = 'mv ' + os.path.abspath(os.path.join(indir, d, f + '/*'))  + \
            ' ' + os.path.abspath(os.path.join(indir, d)) 
        print cmd
        subprocess.Popen(cmd, shell=True)

# Utility function to shelve the ensGene.map or locuslink.map files.
def shelveEventToInfo(txt_f, db_f):
   
    db = shelve.open(db_f, writeback=True)
    for line in open(txt_f):         
        if not line.startswith("#"):
            event, gene, symb, desc = line.strip().split("\t")
            db[event] = [gene, symb, desc]
    db.close()

def shelveEventToInfo2(map_f, lookup_f, db_f):
    
    db = shelve.open(db_f, writeback=True)
    geneToInfo = {}
    for line in open(lookup_f):
        if not line.startswith("#"):
            vals = line.strip().split("\t")
            symb = 'n/a'
            desc = 'n/a'
            if len(vals) == 3:
                gene, symb, desc = vals
            elif len(vals) == 2:
                gene, symb = vals
            geneToInfo[gene] = [gene, symb, desc]

    for line in open(map_f):
        if not line.startswith("event_id"):
            event, gene = line.strip().split("\t")
            if gene in geneToInfo:
                db[event] = geneToInfo[gene]
    db.close() 

def shelveEventToInfoFromV2(gff_f, lookup_f, db_f):
    
    db = shelve.open(db_f, writeback=True)
    geneToInfo = {}
    for line in open(lookup_f):
        if not line.startswith("#"):
            vals = line.strip().split("\t")
            symb = 'n/a'
            desc = 'n/a'
            if len(vals) == 3:
                gene, symb, desc = vals
            elif len(vals) == 2:
                gene, symb = vals
            geneToInfo[gene] = [gene, symb, desc]

    for line in open(gff_f):
        if not line.startswith("#"):
            vals = line.strip().split("\t")
            if vals[2] == 'gene':
                attr = vals[8].split(";")
                for pair in attr:
                    category, value = pair.split("=")
                    if category == 'Name':
                        event = value
                    if category == 'ensg_id':
                        gene = value
                        if "," in gene:
                            gene = gene.split(",")[0]
                        if gene != 'NA' and gene in geneToInfo:
                            db[event] = geneToInfo[gene]

    db.close() 


# Consolidate all summary information into 1 file for subsequent analysis and plotting.
def consolidateSummaries(summarydir, groups_f, out_f):
    
    import miso_utils
   
    groupToSamples = {}
    samples = []
    for line in open(groups_f):
        sample, g = line.strip().split()
        try:
            groupToSamples[g].append(sample) 
        except:
            groupToSamples[g] = [sample]
        samples.append(sample)
    groups = groupToSamples.keys()
    print len(samples), 'samples'

    # Iterate through all samples and save psi values
    eventMaster = {}    # event -> sample -> psi values and bfs
    for s in samples:
        print s
        f = os.path.join(summarydir, s + ".miso_summary")
        eventToInfo, header = miso_utils.getSummarySingle(f)
        for e in eventToInfo:
            if e not in eventMaster:
                eventMaster[e] = {}
            eventMaster[e][s] = [eventToInfo[e]['low'][0], eventToInfo[e]['mean'][0],\
                eventToInfo[e]['high'][0]]

    # Iterate through all comparisons and save bfs
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            f1 = os.path.join(summarydir, samples[i] + "_vs_" + samples[j] + ".miso_bf")
            f2 = os.path.join(summarydir, samples[j] + "_vs_" + samples[i] + ".miso_bf")
            if os.path.exists(f1):
                f = f1
            else:
                f = f2
            print f
            eventToInfo, header = miso_utils.getSummary(f)
            for e in eventToInfo:
                eventMaster[e][samples[i] + "_vs_" + samples[j]] = eventToInfo[e]['max bf']
                if 'gene' not in eventMaster[e]:
                    eventMaster[e]['gene'] = eventToInfo[e]['gene']
                    eventMaster[e]['symb'] = eventToInfo[e]['symb']
                    eventMaster[e]['desc'] = eventToInfo[e]['desc']

    out = open(out_f, 'w')
    out.write("#Event\t")
    for i in range(len(samples)):
        out.write("%s\t%s\t%s\t"%(samples[i] + "_low", samples[i] + "_mean",\
            samples[i] + "_high"))
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)): 
            out.write(samples[i] + "_vs_" + samples[j] + "\t")
    out.write("gene\tsymb\tdesc\n")
    for e in eventMaster:
        out.write(e + "\t")
        for i in range(len(samples)):
            if samples[i] in eventMaster[e]:
                out.write("\t".join(map(str, eventMaster[e][samples[i]])) + "\t")
            else:
                out.write("n/a\tn/a\tn/a\t")
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                compname = samples[i] + "_vs_" + samples[j]
                if compname in eventMaster[e]:
                    out.write(str(eventMaster[e][compname]) + "\t")
                else:
                    out.write("1\t")
        try:
            out.write(eventMaster[e]['gene'] + "\t" + eventMaster[e]['symb'] + "\t" + \
                eventMaster[e]['desc'] + "\n")
        except:
            out.write("n/a\tn/a\tn/a\n")
    out.close()

# Consolidate all summary information into 1 file for subsequent analysis and plotting.
# Consider isoforms of multi-isoform events separately
def consolidateSummariesMultiIso(summarydir, groups_f, out_f):
    
    import miso_utils
   
    groupToSamples = {}
    samples = []
    for line in open(groups_f):
        sample, g = line.strip().split()
        try:
            groupToSamples[g].append(sample) 
        except:
            groupToSamples[g] = [sample]
        samples.append(sample)
    groups = groupToSamples.keys()
    print len(samples), 'samples'

    # Iterate through all samples and save psi values
    eventMaster = {}    # event -> sample -> psi values and bfs
    for s in samples:
        print s
        f = os.path.join(summarydir, s + ".miso_summary")
        eventToInfo, header = miso_utils.getSummarySingle(f)
        for e in eventToInfo:
            if e not in eventMaster:
                eventMaster[e] = {}
            eventMaster[e][s] = [\
                ",".join(map(str, eventToInfo[e]['low'])),\
                ",".join(map(str, eventToInfo[e]['mean'])),\
                ",".join(map(str, eventToInfo[e]['high']))]

    # Iterate through all comparisons and save bfs
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            f1 = os.path.join(summarydir, samples[i] + "_vs_" + samples[j] + ".miso_bf")
            f2 = os.path.join(summarydir, samples[j] + "_vs_" + samples[i] + ".miso_bf")
            if os.path.exists(f1):
                f = f1
            else:
                f = f2
            print f
            eventToInfo, header = miso_utils.getSummary(f)
            for e in eventToInfo:
                eventMaster[e][samples[i] + "_vs_" + samples[j]] = \
                    ",".join(map(str, eventToInfo[e]['bayes factor']))
                if 'gene' not in eventMaster[e]:
                    eventMaster[e]['gene'] = eventToInfo[e]['gene']
                    eventMaster[e]['symb'] = eventToInfo[e]['symb']
                    eventMaster[e]['desc'] = eventToInfo[e]['desc']

    out = open(out_f, 'w')
    out.write("#Event\t")
    for i in range(len(samples)):
        out.write("%s\t%s\t%s\t"%(samples[i] + "_low", samples[i] + "_mean",\
            samples[i] + "_high"))
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)): 
            out.write(samples[i] + "_vs_" + samples[j] + "\t")
    out.write("gene\tsymb\tdesc\n")
    for e in eventMaster:
        out.write(e + "\t")
        for i in range(len(samples)):
            if samples[i] in eventMaster[e]:
                out.write("\t".join(map(str, eventMaster[e][samples[i]])) + "\t")
            else:
                out.write("n/a\tn/a\tn/a\t")
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                compname = samples[i] + "_vs_" + samples[j]
                if compname in eventMaster[e]:
                    out.write(str(eventMaster[e][compname]) + "\t")
                else:
                    out.write("1\t")
        try:
            out.write(eventMaster[e]['gene'] + "\t" + eventMaster[e]['symb'] + "\t" + \
                eventMaster[e]['desc'] + "\n")
        except:
            out.write("n/a\tn/a\tn/a\n")
    out.close()


# Find events that change monotonically (significantly).
# E.g. only look at significant comparisons
def monotonic(consolidated_f, groups_f, minbf, nshuffles, out_f):
    
    import miso_utils
    minbf = float(minbf)
    nshuffles = int(nshuffles)
   
    groupToSamples = {}
    grouporder = []
    samples = []
    for line in open(groups_f):
        sample, g = line.strip().split()
        samples.append(sample)
        if g not in grouporder:
            grouporder.append(g)
        try:
            groupToSamples[g].append(sample) 
        except:
            groupToSamples[g] = [sample]
    groups = groupToSamples.keys()
    ctllabel = grouporder[0]
    explabel = grouporder[-1]

    # Identify between-group blocks and get their indices for reference 
    groupidx = [] 
    lowidx = 0
    for g in grouporder:
        n = len(groupToSamples[g]) 
        groupidx.append([lowidx, lowidx + n])
        lowidx += n
    compidx = []
    for i in range(len(groups)):
        idx1 = groupidx[i]
        for j in range(i + 1, len(groups)):
            # Analyze each block of groups
            idx2 = groupidx[j]
            for k in range(idx1[0], idx1[1]):
                for l in range(idx2[0], idx2[1]):
                    compidx.append([k, l])
    compidx = zip(*compidx)

    data = []
    print "Reading data..."
    counter = 0
    for line in open(consolidated_f):
        vals = line.strip().split("\t")
        if line.startswith("#"):
            header = vals
        else:
            sampleToBF = {}
            event = vals[0]
            symb = vals[-2]
            sampleToPsi = {}
            skipme = False
            for i in range(len(header)):
                if "_mean" in header[i]:
                    sample = header[i][:-5] 
                    if sample in samples:
                        if vals[i] == 'n/a':
                            skipme = True
                            break
                        else:
                            sampleToPsi[sample] = float(vals[i])
                elif "_vs_" in header[i]:
                    s1, s2 = header[i].split("_vs_")
                    bf = float(vals[i])
                    if s1 not in sampleToBF:
                        sampleToBF[s1] = {}
                    if s2 not in sampleToBF:
                        sampleToBF[s2] = {}
                    sampleToBF[s1][s2] = bf
                    sampleToBF[s2][s1] = bf

            if not skipme:
                mat = zeros((len(samples), len(samples)), dtype='int')
                for i in range(len(samples)):
                    for j in range(i + 1, len(samples)):
                        try:
                            psi1 = sampleToPsi[samples[i]]
                            psi2 = sampleToPsi[samples[j]]
                            bf = sampleToBF[samples[i]][samples[j]]
                            dpsi = psi2 - psi1
                            if bf >= minbf:
                                mat[i, j] = sign(dpsi) 
                        except:
                            pass
               
                shuffledsamples = list(samples)
                mat = zeros((len(samples), len(samples), nshuffles), dtype='int') 
                for n in range(nshuffles + 1):
                    for i in range(len(samples)):
                        for j in range(i + 1, len(samples)):
                            try:
                                psi1 = sampleToPsi[shuffledsamples[i]]
                                psi2 = sampleToPsi[shuffledsamples[j]]
                                bf = sampleToBF[shuffledsamples[i]][shuffledsamples[j]]
                                dpsi = psi2 - psi1
                                if bf >= minbf:
                                    mat[i, j, n] = sign(dpsi) 
                            except:
                                pass
                    shuffle(shuffledsamples)
                
                # Get the mean value in exp vs. ctl
                ctlpsi = array([sampleToPsi[x] for x in sampleToPsi if \
                    x in groupToSamples[ctllabel]])
                exppsi = array([sampleToPsi[x] for x in sampleToPsi if \
                    x in groupToSamples[explabel]])

                signs = mat[compidx[0], compidx[1], :]
                signs = signs.sum(axis=0)
                m = signs[1:].mean()
                stdev = signs[1:].std()
                z = (signs[0] - m) / stdev
                if signs[0] == 0:
                    z = 0
                data.append([event, symb, signs[0], m, stdev,\
                    round(ctlpsi.mean(), 2), round(exppsi.mean(), 2), \
                    round(exppsi.mean() - ctlpsi.mean(), 2), z])
                print counter 
                counter += 1

    data.sort(key=operator.itemgetter(-1), reverse=True)
    out = open(out_f, 'w')
    out.write("#Event\tSymb\tTrueval\tMean\tStd\t%s_psi\t%s_psi\tdelta_psi\tZ-score\n"%(\
        ctllabel, explabel))
    for item in data:
        out.write("\t".join(map(str, item)) + "\n")
    out.close()
  
# Find events that change monotonically (significantly).
# E.g. only look at significant comparisons
def monotonicMultiIso(consolidated_f, groups_f, minbf, nshuffles, out_f):
    
    import miso_utils
    minbf = float(minbf)
    nshuffles = int(nshuffles)
   
    groupToSamples = {}
    grouporder = []
    samples = []
    for line in open(groups_f):
        sample, g = line.strip().split()
        samples.append(sample)
        if g not in grouporder:
            grouporder.append(g)
        try:
            groupToSamples[g].append(sample) 
        except:
            groupToSamples[g] = [sample]
    groups = groupToSamples.keys()
    ctllabel = grouporder[0]
    explabel = grouporder[-1]

    # Identify between-group blocks and get their indices for reference 
    groupidx = [] 
    lowidx = 0
    for g in grouporder:
        n = len(groupToSamples[g]) 
        groupidx.append([lowidx, lowidx + n])
        lowidx += n
    compidx = []
    for i in range(len(groups)):
        idx1 = groupidx[i]
        for j in range(i + 1, len(groups)):
            # Analyze each block of groups
            idx2 = groupidx[j]
            for k in range(idx1[0], idx1[1]):
                for l in range(idx2[0], idx2[1]):
                    compidx.append([k, l])
    compidx = zip(*compidx)

    data = []
    print "Reading data..."
    counter = 0
    for line in open(consolidated_f):
        vals = line.strip().split("\t")
        if line.startswith("#"):
            header = vals
        else:
            sampleToBF = {}
            event = vals[0]
            symb = vals[-2]
            sampleToPsi = {}
            skipme = False
            niso = 0

            # Record each sample
            for i in range(len(header)):
                if "_mean" in header[i]:
                    sample = header[i][:-5] 
                    if sample in samples:
                        if vals[i] == 'n/a':
                            skipme = True
                            break
                        else:
                            if "," not in vals[i]:
                                sampleToPsi[sample] = [float(vals[i])]
                                niso = 1
                            else:
                                sampleToPsi[sample] = map(float, vals[i].split(","))
                                niso = len(sampleToPsi[sample])
                elif "_vs_" in header[i]:
                    s1, s2 = header[i].split("_vs_")
                    if "," not in vals[i]:
                        bf = [float(vals[i])]
                    else:
                        bf = map(float, vals[i].split(","))
                    if s1 not in sampleToBF:
                        sampleToBF[s1] = {}
                    if s2 not in sampleToBF:
                        sampleToBF[s2] = {}
                    sampleToBF[s1][s2] = bf
                    sampleToBF[s2][s1] = bf

            if not skipme:
                signvals = []
                meanvals = []
                stdvals = []
                ctlvals = []
                expvals = []
                deltavals = []
                zvals = []
                for ni in range(niso):
                    mat = zeros((len(samples), len(samples)), dtype='int')
                    for i in range(len(samples)):
                        for j in range(i + 1, len(samples)):
                            try:
                                psi1 = sampleToPsi[samples[i]][ni]
                                psi2 = sampleToPsi[samples[j]][ni]
                                bf = sampleToBF[samples[i]][samples[j]][ni]
                                dpsi = psi2 - psi1
                                if bf >= minbf:
                                    mat[i, j] = sign(dpsi) 
                            except:
                                pass
                   
                    shuffledsamples = list(samples)
                    mat = zeros((len(samples), len(samples), nshuffles), dtype='int') 
                    for n in range(nshuffles + 1):
                        for i in range(len(samples)):
                            for j in range(i + 1, len(samples)):
                                try:
                                    psi1 = sampleToPsi[shuffledsamples[i]][ni]
                                    psi2 = sampleToPsi[shuffledsamples[j]][ni]
                                    bf = sampleToBF[shuffledsamples[i]][shuffledsamples[j]][ni]
                                    dpsi = psi2 - psi1
                                    if bf >= minbf:
                                        mat[i, j, n] = sign(dpsi) 
                                except:
                                    pass
                        shuffle(shuffledsamples)
                    
                    # Get the mean value in exp vs. ctl
                    ctlpsi = array([sampleToPsi[x][ni] for x in sampleToPsi if \
                        x in groupToSamples[ctllabel]])
                    exppsi = array([sampleToPsi[x][ni] for x in sampleToPsi if \
                        x in groupToSamples[explabel]])

                    signs = mat[compidx[0], compidx[1], :]
                    signs = signs.sum(axis=0)
                    m = signs[1:].mean()
                    stdev = signs[1:].std()
                    z = (signs[0] - m) / stdev

                    signvals.append(signs[0])
                    meanvals.append(round(m, 2))
                    stdvals.append(round(stdev, 2))
                    ctlvals.append(round(ctlpsi.mean(), 2))
                    expvals.append(round(exppsi.mean(), 2))
                    deltavals.append(round(exppsi.mean() - ctlpsi.mean(), 2))
                    zvals.append(z)

                if max([abs(x) for x in signvals]) > 0:
                    data.append([event, symb, 
                        ",".join(map(str, signvals)),\
                        ",".join(map(str, meanvals)),
                        ",".join(map(str, stdvals)),
                        ",".join(map(str, ctlvals)),
                        ",".join(map(str, expvals)),
                        ",".join(map(str, deltavals)),
                        ",".join(map(str, zvals)),
                        max([abs(z) for z in zvals])])
                print counter 
                counter += 1

    data.sort(key=operator.itemgetter(-1), reverse=True)
    out = open(out_f, 'w')
    out.write("#Event\tSymb\tTrueval\tMean\tStd\t%s_psi\t%s_psi\tdelta_psi\tZ-score\tMaxZ\n"%(\
        ctllabel, explabel))
    for item in data:
        out.write("\t".join(map(str, item)) + "\n")
    out.close()


# Consolidate all summary information into 1 file for subsequent analysis and plotting.
def consolidateSummariesNoBF(summarydir, eventToGeneInfo_f, out_f):
    
    import miso_utils
    eventToGeneInfo = shelve.open(eventToGeneInfo_f, 'r') 
     
    samples = [s.split(".")[0] for s in os.listdir(summarydir)]
    samples.sort()

    # Iterate through all samples and save psi values
    eventMaster = {}    # event -> sample -> psi values and bfs
    for s in samples:
        print s
        f = os.path.join(summarydir, s + ".miso_summary")
        eventToInfo, header = miso_utils.getSummarySingle(f)
        for e in eventToInfo:
            if e not in eventMaster:
                eventMaster[e] = {}
            eventMaster[e][s] = [eventToInfo[e]['low'][0], eventToInfo[e]['mean'][0],\
                eventToInfo[e]['high'][0]]

    out = open(out_f, 'w')
    out.write("#Event\t")
    for i in range(len(samples)):
        out.write("%s\t%s\t%s\t"%(samples[i] + "_low", samples[i] + "_mean",\
            samples[i] + "_high"))
    out.write("gene\tsymb\tdesc\n")
    for e in eventMaster:
        out.write(e + "\t")
        for i in range(len(samples)):
            if samples[i] in eventMaster[e]:
                out.write("\t".join(map(str, eventMaster[e][samples[i]])) + "\t")
            else:
                out.write("n/a\tn/a\tn/a\t")
        try:
            out.write("\t".join(eventToGeneInfo[e]) + "\n")
        except:
            out.write("n/a\tna/a\tn/a\n")
    out.close()

def addSeqsToConsolidated(consolidated_f, seqs_fa, out_f):
   
    eventToSeqs = {} 
    for line in open(seqs_fa):
        name, seq = line.strip().split("\t")
        eventfields = name.split(".")
        event = eventfields[0]
        region = eventfields[-1]
        if len(event.split("@")) == 3: 
            if event not in eventToSeqs:
                eventToSeqs[event] = ['', '', ''] 
            if region == 'up': 
                eventToSeqs[event][0] = seq
            elif region == 'se': 
                eventToSeqs[event][1] = seq
            elif region == 'dn': 
                eventToSeqs[event][2] = seq
        elif len(event.split("@")) == 4: 
            if event not in eventToSeqs:
                eventToSeqs[event] = ['', '', '', ''] 
            if region == 'up': 
                eventToSeqs[event][0] = seq
            elif region == 'mxe1': 
                eventToSeqs[event][1] = seq
            elif region == 'mxe2': 
                eventToSeqs[event][2] = seq
            elif region == 'dn': 
                eventToSeqs[event][3] = seq
        elif len(event.split("@")) == 2 and len(event.split("|")) == 2 and "|" in event.split("@")[1]:
            if event not in eventToSeqs:
                eventToSeqs[event] = ['', ''] 
            if region == 'up': 
                eventToSeqs[event][0] = seq
            elif region == 'altss' and eventfields[1] == 'B': 
                eventToSeqs[event][1] = seq
        elif len(event.split("@")) == 2 and len(event.split("|")) == 2 and "|" in event.split("@")[0]:
            if event not in eventToSeqs:
                eventToSeqs[event] = ['', ''] 
            if region == 'altss' and eventfields[1] == 'B': 
                eventToSeqs[event][0] = seq
            elif region == 'dn': 
                eventToSeqs[event][1] = seq
        elif len(event.split("@")) == 2:
            if event not in eventToSeqs:
                eventToSeqs[event] = ['', ''] 
            if region == 'up':
                eventToSeqs[event][0] = seq
            elif region == 'dn': 
                eventToSeqs[event][1] = seq
            
    out = open(out_f, 'w')
    for line in open(consolidated_f):
        if line.startswith("#"):
            out.write(line.strip() + "\tSequences\n")
        else:
            vals = line.strip().split("\t")
            event = vals[0]
            if event in eventToSeqs:
                vals.extend(eventToSeqs[event])
                out.write("\t".join(vals) + "\n")
    out.close()



def sashimi(event, indexdir, settings_f, outdir):
   
    cmd = "python /home/et_wang/Tools/pythonmodules/lib/python2.7/site-packages/misopy/sashimi_plot/plot.py "+\
        " --plot-event '" + event + "' " + indexdir + ' ' + settings_f + ' --output-dir ' +\
        outdir
    print cmd
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    
    
def sashimiWrapper(monotonic_f, indexdir, settings_f, outdir, nevents, mindpsi):

    mindpsi = float(mindpsi)
    nevents = int(nevents)
    allevents = []
    for line in open(monotonic_f):
        if not line.startswith("#"):
            event, symb, trueval, m, std, psi1, psi2, deltapsi, z = \
                line.strip().split("\t")
            deltapsi = float(deltapsi)
            z = abs(float(z))
            if abs(deltapsi) >= mindpsi:
                allevents.append([event, z])
    allevents.sort(key=operator.itemgetter(-1), reverse=True)
    for i in range(nevents):
        print allevents[i][0], allevents[i][1]
        sashimi(allevents[i][0], indexdir, settings_f, outdir) 


# Convert IDs from 1 build to another build
# Just requires event ID to be in first column of in_f
def convertIDsBtwBuild(in_f, gff3_f, out_f):

    eventToCoords = {}
    for line in open(gff3_f):
        if not line.startswith("#"):
            vals = line.strip().split("\t")
            chrom, etype, ecat, start, end, score, strand, frame, attr = vals
            if ecat == 'gene':
                event = attr.split(";")[0].split("=")[1]
                eventToCoords[event] = [chrom, start, end, strand]

    out = open(out_f, 'w')
    for line in open(in_f):
        if line.startswith("#"):
            vals = line.strip().split("\t")
            vals.insert(1, 'mm10_strand')
            vals.insert(1, 'mm10_end')
            vals.insert(1, 'mm10_start')
            vals.insert(1, 'mm10_chrom')
            out.write("\t".join(vals) + "\n")
        else:
            vals = line.strip().split("\t")
            event = vals[0]
            chrom, start, end, strand = eventToCoords[event]
            vals.insert(1, strand)
            vals.insert(1, end)
            vals.insert(1, start)
            vals.insert(1, chrom)
            out.write("\t".join(vals) + "\n")
    out.close()


# Simple script to get IDs using a new lookup table.    
def getIDsForMonotonic(old_f, lookup_f, new_f):
   
    eventToInfo = shelve.open(lookup_f, 'r')

    out = open(new_f, 'w') 
    for line in open(old_f):
        if line.startswith("#"):
            out.write(line)
        else:
            vals = line.strip().split("\t")
            event = vals[0]
            if event in eventToInfo:
                gene, symb, des = eventToInfo[event]
                vals[1] = symb
                out.write("\t".join(vals) + "\n")
            else:
                out.write(line) 
    out.close()
 
def getIDsForAFEs(gff3_f, xref_f, out_f):
   
    llToInfo = {} 
    for line in open(xref_f):
        tx, ll, symb, desc = line.strip().split("\t")
        llToInfo[ll] = [symb, desc]

    out = open(out_f, 'w')
    out.write("#Event\tGene\tSymb\tDesc\n") 
    for line in open(gff3_f):
        if not line.startswith("#"):
            vals = line.strip().split("\t")
            if vals[2] == 'gene':
                attr = vals[8]
                event = attr.split(";")[0].split("=")[1]
                ll = event.split("@")[0]
                symb = 'n/a' 
                desc = 'n/a' 
                if ll in llToInfo:
                    symb, desc = llToInfo[ll]
                out.write("\t".join([event, ll, symb, desc]) + "\n")
    out.close()
 

