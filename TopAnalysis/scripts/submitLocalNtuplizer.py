import os
import sys
import subprocess
import optparse
import json

def getDatasetComponents(opt):

    """ query components of the dataset """

    fList=[]

    #get files in the dataset
    print 'Querying',opt.dataset
    p = subprocess.Popen(['dasgoclient --query=\"file dataset=%s\"'%opt.dataset], 
                         stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE,
                         shell=True)
    out, err = p.communicate()

    dataset_files=out.split()
    nfiles=len(dataset_files)
    print 'I have %d jobs to submit...adding extra info'%nfiles

    #get parents, if required
    for i in range(nfiles):
        x=dataset_files[i]

        sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(nfiles))))
        sys.stdout.flush()

        if opt.addParent:
            p = subprocess.Popen(['dasgoclient --query=\"parent file=%s\"'%x],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True)
            out, err = p.communicate()
            fList.append( (x,out.split()) )
        elif opt.addAODParent:
            fList = getListAOD(x,opt.dataset)
        else:
            fList.append( (x,[]) )
    
    return fList


def buildCondorFile(opt,fList,FarmDirectory):

    """ builds the condor file to submit the ntuplizer """

    jobTag=opt.jobTag

    cmssw=os.environ['CMSSW_BASE']

    secFileType=None

    #condor submission file
    condorFile='%s/condor_%s.sub'%(FarmDirectory,opt.jobTag)
    print '\nWrites: ',condorFile
    with open (condorFile,'w') as condor:
        condor.write('executable = {0}/worker_{1}.sh\n'.format(FarmDirectory,opt.jobTag))
        condor.write('output     = {0}/output_{1}.out\n'.format(FarmDirectory,opt.jobTag))
        condor.write('error      = {0}/output_{1}.err\n'.format(FarmDirectory,opt.jobTag))
        condor.write('log        = {0}/output_{1}.log\n'.format(FarmDirectory,opt.jobTag))
        condor.write('+JobFlavour = "nextweek"\n')
        OpSysAndVer = str(os.system('cat /etc/redhat-release'))
        if 'SLC' in OpSysAndVer:
            OpSysAndVer = "SLCern6"
        else:
            OpSysAndVer = "CentOS7"
        condor.write('requirements = (OpSysAndVer =?= "{0}")\n'.format(OpSysAndVer)) 

        for i in range(len(fList)):
            secFileList=','.join(fList[i][1])
            if '/RAW' in secFileList: secFileType='RAW'
            if '/AOD' in secFileList: secFileType='AOD'
            condor.write('arguments=%d %s %s\n'%(i,fList[i][0],secFileList))
            condor.write('queue 1\n')
                                
    #local worker script
    print 'Secondary file type',secFileType
    if secFileType == 'RAW':
        opt.extraOpts += ' redoProtonRecoFromRAW=True'
    if secFileType == 'AOD':
        opt.extraOpts += ' runWithAOD=True'
    print 'Additional extra opts will be set to',opt.extraOpts

    workerFile='%s/worker_%s.sh'%(FarmDirectory,opt.jobTag)
    with open(workerFile,'w') as worker:
        worker.write('#!/bin/bash\n')
        worker.write('WORKDIR=`pwd`\n')
        worker.write('echo "Working directory is ${WORKDIR}"\n')
        worker.write('cd %s\n'%cmssw)
        worker.write('eval `scram r -sh`\n')
        worker.write('cd ${WORKDIR}\n')
        worker.write('opts="lumiJson=%s inputFile=${2}"\n'%opt.lumiMask)
        worker.write('if [ ! -z "${3}" ]; then\n')
        worker.write('  opts="${opts} secInputFile=${3}"\n')
        worker.write('fi\n')
        if opt.proxy:
            worker.write('export X509_USER_PROXY=%s/myproxy509\n'%FarmDirectory)
        else:
            worker.write('echo "No proxy has been configured"\n')
        worker.write('opts="${opts} %s"\n'%opt.extraOpts)
        worker.write('echo "Running cmsRun with the following options: ${opts}"\n')
        worker.write('cmsRun %s/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py ${opts}\n'%cmssw)
        worker.write('xrdcp --force MiniEvents.root root://eoscms//%s/%s/MiniEvents_${1}.root\n'%(opt.output,opt.jobTag))
        worker.write('rm MiniEvents.root\n')

    return condorFile

def getListAOD(filename, miniAODds):
  print 'collect parrents for ',filename

  recotag=filename.split('/')[6]
  AODds=miniAODds.replace('MINI','')
  
  #Get list of runs / lumis for the input file
  p = subprocess.Popen(['dasgoclient --query="lumi file=%s" -json'%filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  out, err = p.communicate()
  data = json.loads(out)
  miniaodruns=[]
  miniaodlumis=[]
  for i in range(len(data)):
    miniaodruns.append(data[i]['lumi'][0]['run_number'])
    miniaodlumis.append(data[i]['lumi'][0]['lumi_section_num'])
	
  #Get list of RAW files
  dascmd='dasgoclient --query="parent file=%s"'%filename
  p = subprocess.Popen([dascmd],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  out, err = p.communicate()
  rawfiles = out.split()
  
  #Loop over all RAW files to locate AODs that match miniAOD lumiblocks  
  AODfiles=[]
  for rawfile in rawfiles:
    dascmd='dasgoclient --query="child file=%s" | grep /AOD/%s'%(rawfile,recotag)
    p = subprocess.Popen([dascmd],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    
	#Check AOD files to match with at least one lumiblock
    if len(out.split()) == 1: #if only one AOD child, skip the lumi-check
      if out.split()[0] not in AODfiles: AODfiles.append(out.split()[0])	
    else:  
      for aodfile in out.split():
        if aodfile in AODfiles: continue
        if(MatchedAOD(aodfile,miniaodruns,miniaodlumis)): AODfiles.append(aodfile);

  return AODfiles

def MatchedAOD(aodfile,miniaodruns,miniaodlumis):
  
  aodruns=json.loads(subprocess.check_output(['dasgoclient --query="run file=%s"'%aodfile], shell=True).strip())
  if len(aodruns)==1 and aodruns[0] not in miniaodruns: return False
  aodlumis=[json.loads(x) for x in subprocess.check_output(['dasgoclient --query="lumi file=%s"'%aodfile], shell=True).split()]
  for idx, aodrun in enumerate(aodruns): 
    # check location of the run in the miniaodruns list
    i_run=-1
    for ii, miniaodrun in enumerate(miniaodruns):
      if aodrun == miniaodrun: i_run=ii; break
    if(i_run<0): continue # run not in the list of miniaodruns
	
	#check if have at least a single lumiblock in list of miniaodlumis
    for lumi in aodlumis[idx]:
      if lumi in miniaodlumis[i_run]: return True
	  
  # if faled to return true -> AOD is not matched to the miniAODds
  return False
  
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--jobTag',      
                      dest='jobTag',      
                      help='job tag [%default]',  
                      default='MiniEvents',   
                      type='string')
    parser.add_option('--proxy',
                      dest='proxy',      
                      help='start proxy [%default]',  
                      default=False,
                      action='store_true')
    parser.add_option('--dryRun',
                      dest='dryRun',      
                      help='dry run (do not submit jobs to condor) [%default]',  
                      default=False,
                      action='store_true')
    parser.add_option('--output',    
                      dest='output',   
                      help='output to store [%default]',      
                      default='/store/cmst3/group/top/RunIIReReco/',
                      type='string')
    parser.add_option('--dataset',  
                      dest='dataset', 
                      help='dataset to process [%default]',   
                      default=None,    
                      type='string')
    parser.add_option('--addParent',   
                      dest='addParent',  
                      help='add parent [%default]',   
                      default=False,   
                      action='store_true')
    parser.add_option('--addAODParent',   
                      dest='addAODParent',  
                      help='add AOD parent [%default]',   
                      default=False,   
                      action='store_true')
    parser.add_option('--extraOpts',    
                      dest='extraOpts',
                      help='extra options to use in the ntuplizer cfg [%default]',
                      default='runOnData=True,era=era2017,runWithAOD=True',
                      type='string')
    parser.add_option('--lumiMask',        
                      dest='lumiMask',    
                      help='json with list of good lumis', 
                      default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt')
    (opt, args) = parser.parse_args()

    opt.extraOpts=' '.join(opt.extraOpts.split(','))

    #prepare directory with scripts
    cmssw=os.environ['CMSSW_BASE']
    FarmDirectory='%s/FarmLocalNtuple'%cmssw
    os.system('mkdir -p '+FarmDirectory)

    #start proxy
    if opt.proxy:
        os.system('voms-proxy-init --voms cms --out %s/myproxy509'%FarmDirectory)
        os.environ['X509_USER_PROXY']='%s/myproxy509'%FarmDirectory

    #get file list
    fList=getDatasetComponents(opt)

    #build condor submission script and launch jobs
    condor=buildCondorFile(opt,fList,FarmDirectory)
    if not opt.dryRun:
        os.system('condor_submit %s'%condor)

if __name__ == "__main__":
    sys.exit(main())
