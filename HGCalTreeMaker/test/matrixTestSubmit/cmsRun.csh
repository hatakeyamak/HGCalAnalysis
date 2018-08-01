#!/bin/tcsh
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
setenv PATH /usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
printenv
source /cvmfs/cms.cern.ch/cmsset_default.csh  ## if a bash script, use .sh instead of .csh
echo "source done"
#
# TO BE CHECKED before submission
#
setenv MYCMSSW CMSSW_10_3_X_2018-07-31-2300 ## <========= TO BE CHECKED
setenv VERSION v01
#
printf "Start time: "; /bin/date
printf "Job is running on node: "; /bin/hostname
printf "Job running as user: "; /usr/bin/id
printf "Job is running in directory: "; /bin/pwd
voms-proxy-info
voms-proxy-info -fqan
#
echo ${MYCMSSW} ${VERSION}
### for case 1. EOS have the following line, otherwise remove this line in case 2.
xrdcp -s root://kodiak-se.baylor.edu//store/user/hatake/condor/tarballs/${MYCMSSW}_HGCAL_condor.tgz .
tar -xf ${MYCMSSW}_HGCAL_condor.tgz
rm ${MYCMSSW}_HGCAL_condor.tgz
#setenv SCRAM_ARCH slc6_amd64_gcc530
ls -R
cd ${MYCMSSW}/src
#
printf "Start time: "; /bin/date
printf "Job is running on node: "; /bin/hostname
printf "Job running as user: "; /usr/bin/id
printf "Job is running in directory: "; /bin/pwd
voms-proxy-info
voms-proxy-info -fqan
#
scramv1 b ProjectRename
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
#bash cmsRun ../../pi50.py maxEvents=2000 skipEvents=`echo ${1}\*2000|bc`
cmsRun ../../step2_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT_condor.py inputFiles='root://kodiak-se.baylor.edu//store/user/hatake/crab_outputs/TTbar_14TeV/CMSSW_10_3_X_2018-07-30-2300_Step1_v1/180731_232243/0000/step1_'${1}'.root' outputFile='file:step2_'${1}'.root' maxEvents=10
cmsRun ../../step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT_VALIDATION_DQM_condor.py inputFiles='file:step2_'${1}'.root' outputFile='file:step3_'${1}'.root' outputMINIAOD='file:step3_inMINIAOD_'${1}'.root' outputDQM='file:step3_inDQM_'${1}'.root'  maxEvents=10
foreach f (`ls *.root`)
   echo $f
   set name=`basename $f .root`
   echo $name
   gfal-copy --just-copy ${f} gsiftp://kodiak-se.baylor.edu/cms/data/store/user/hatake/condor/outputs/${name}_${MYCMSSW}_${1}_${VERSION}.root
end
### gfal-copy --just-copy pi50_trees_MCfull.root gsiftp://kodiak-se.baylor.edu/cms/data/store/user/hatake/condor/pi50_trees_MCfull_${MYCMSSW}_${1}_${VERSION}.root
### remove the output file if you don't want it automatically transferred when the job ends
### rm pi50_trees_MCfull.root
cd ${_CONDOR_SCRATCH_DIR}
rm -rf ${MYCMSSW}
