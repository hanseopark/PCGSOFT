#!/bin/bash 

# Settings
DATE=`date +%y%m%d`
DAUGHTER=0
ANALYSISNAME="$1"
ANALYSISNAME=${ANALYSISNAME:-"NoNamed"}
WORKDIRNAME="${DATE}_${DAUGHTER}_${ANALYSISNAME}"
LATEST="NULL"

# Directory path
INITDIR=`pwd`
LATESTWORKDIR="${INITDIR}/${LATEST}"
WORKDIR="${INITDIR}/${WORKDIRNAME}"
CERNBOX="NULL" #cernbox directory path
CERNBOXLATEST="NULL"

# Flags
OWNGIT="NO"
LOCAL="$2"
PCGINIT="$3"

# Load cernbox local path
if [ -e ${HOME}/.local/share/data/CERNBox/cernbox.cfg ]
then
	TMP=`grep localPath ${HOME}/.local/share/data/CERNBox/cernbox.cfg`
	CERNBOX="${TMP:22}Main_Analysis"
	echo ${CERNBOX}
fi

# Look for the latest working directory to take over analysis codes
LATEST=`ls | grep ${ANALYSISNAME} | tail -1`
if [ "${LATEST}" != "NULL" ]
then
	echo /////////////////////////////////////////////////////////////////
	echo /////////////////////////////////////////////////////////////////
	echo ///////Latest working directory for ${ANALYSISNAME} = ${LATEST}/ 
	echo ///////Files will be taken over to today\'s//////////////////////            
	echo /////////////////////////////////////////////////////////////////
	echo /////////////////////////////////////////////////////////////////
	LATESTWORKDIR="${INITDIR}/${LATEST}"
fi

while : # Search working directory and make daughter
do
	if [ -e $WORKDIRNAME ]
	then
		echo ${WORKDIRNAME} found.
		#DAUGHTER=`expr $DAUGHTER + 1`
		DAUGHTER=$((DAUGHTER + 1))
		WORKDIRNAME="${DATE}_${DAUGHTER}_${ANALYSISNAME}"
	else
	#	echo Generating new working directory ${WORKDIRNAME}
	#	cd ${INITDIR}
	#	mkdir $WORKDIRNAME
		WORKDIR="${INITDIR}/${WORKDIRNAME}" # Renew working directory full path
		break
	fi
done

# Look for the latest working directory to take over analysis codes in CERNBOX directory
if [ "${CERNBOX}" != "NULL" ]
then
	cd ${CERNBOX}
	CERNBOXLATEST=`ls | grep ${ANALYSISNAME} | tail -1`
	CERNBOXLATEST="${CERNBOX}/${CERNBOXLATEST}"
	echo ${CERNBOXLATEST}
else
	echo Couldn\'t access to ${CERNBOX}. Is the path correct?
fi

# Generate new working directory in initial directory
cd ${INITDIR}
mkdir ${WORKDIRNAME}

## Clone analysiscodes from git repository
#git clone git@github.com:mtakamur/MyAliceAnalysis.git
#mv MyAliceAnalysis ${WORKDIRNAME}
#cd ${WORKDIRNAME}
#git checkout ${ANALYSISNAME}

# Take over analysis codes
if [ "${LOCAL}" = "y" ]
then
	cd ${WORKDIR}
	if [ "${CERNBOXLATEST}" = "NULL" ]
	then
		if [ "${LATEST}" != "NULL" ]
		then
			echo Files are copied from ${LATESTWORKDIR}
			# /////////////////////////////////////////
			# // Add file type you want to take over //
			# /////////////////////////////////////////
			#cp -v -b --suffix=.bak ${LATESTWORKDIR}/*.cxx ./
			#cp -v -b --suffix=.bak ${LATESTWORKDIR}/*.h ./
			#cp -v -b --suffix=.bak ${LATESTWORKDIR}/*.C ./
			cp -v -b --suffix=.bak ${LATESTWORKDIR}/*.sh ./
			cp -v -b --suffix=.bak ${LATESTWORKDIR}/*.txt ./
			#cp -v -b --suffix=.bak ${LATESTWORKDIR}/* ./
		fi
	else
		echo Files are copied from ${CERNBOXLATEST}
		# /////////////////////////////////////////
		# // Add file type you want to take over //
		# /////////////////////////////////////////
		#cp -v -b --suffix=.bak ${CERNBOXLATEST}/*.cxx ./
		#cp -v -b --suffix=.bak ${CERNBOXLATEST}/*.h ./
		#cp -v -b --suffix=.bak ${CERNBOXLATEST}/*.C ./
		cp -v -b --suffix=.bak ${CERNBOXLATEST}/*.sh ./
		cp -v -b --suffix=.bak ${CERNBOXLATEST}/*.txt ./
		#cp -v -b --suffix=.bak ${CERNBOXLATEST}/* ./
	fi
else
	echo File are not taken over.
fi

# Check if take over MyPCGAnalysisSoftware
echo "Do you want to set the directory as MyPCGAnalysisSoftware(git)? (y/n)"
read OWNGIT

if [ "${OWNGIT}" = "y" ]
then
	cd ${WORKDIR}
	#git init
	#git remote add origin https://github.com/mtakamur/MyPCGAnalysisSoftware.git
	#git remote add origin git@github.com:mtakamur/MyPCGAnalysisSoftware.git
	#git add 0_Unfold.sh
	#git commit -m "New File1"
	#git status
	#git remote add origin https://github.com/hanseopark/PCGSOFT 
	#git push -u origin master
	#git fetch origin
	git config --global username "hapark"
	git config --global user.email "hapark@cern.ch"
	git pull origin MyAnalysis	

	echo "git test"
fi

if [ "${PCGINIT}" = "y" ]
then
	cd ${WORKDIR}
	bash /home/alidock/alice/pcg/AnalysisSoftware/prepareResultsDirectory.sh park 

	echo "prepareResultDirectory.sh Execute"
else
	echo prepareResultsDirectory.sh won\'t be executed.
fi
