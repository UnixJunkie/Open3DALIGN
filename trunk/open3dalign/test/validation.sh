#!/usr/bin/env bash

#
# Usage:
#
# ./validation.sh [{single | multi | both}]
#
# validation statistics are printed on stdout
#

abrupt_exit()
{
	echo
	echo "Validation aborted."
	exit
}


trap abrupt_exit SIGTSTP SIGINT SIGTERM SIGKILL
awk_exe=`which 2> /dev/null gawk | grep -v 'no gawk'`
if [ -z  $awk_exe ]; then
	awk_exe=`which 2> /dev/null awk | grep -v 'no awk'`
fi
if [ -z $awk_exe ]; then
	echo "Cannot find AWK."
	abrupt_exit
	exit
fi
if [ -z $O3A_EXE ]; then
	O3A_EXE=`which 2> /dev/null open3dalign | grep -v 'no open3dalign'`
fi
if [ -z $O3A_EXE ] || [ ! -e $O3A_EXE ]; then
	echo "Cannot find open3dalign binary. Please set the O3A_EXE environment variable and resubmit."
	abrupt_exit
	exit
fi
O3A_PATH=`dirname 2> /dev/null $O3A_EXE`
if [ -z $BABEL_PATH ]; then
	BABEL_PATH=`which 2> /dev/null babel | grep -v 'no babel'`
	BABEL_PATH=`dirname 2> /dev/null $BABEL_PATH`
fi
if [ -z $BABEL_PATH ]; then
	BABEL_PATH=$O3A_PATH
fi
babel_exe=""
if [ -e $BABEL_PATH/babel ]; then
	babel_exe=$BABEL_PATH/babel
fi
if [ -z $babel_exe ]; then
	echo "Cannot find BABEL binaries. Please set the BABEL_PATH environment variable and resubmit."
	abrupt_exit
	exit
fi
$babel_exe -L sdf >& /dev/null
if [ $? != 0 ]; then
	if [ -z $LD_LIBRARY_PATH ]; then
		export LD_LIBRARY_PATH=$BABEL_PATH
	else
		export LD_LIBRARY_PATH=${BABEL_PATH}:${LD_LIBRARY_PATH}
	fi
	if [ -z $LD_32_LIBRARY_PATH ]; then
		export LD_32_LIBRARY_PATH=$BABEL_PATH
	else
		export LD_32_LIBRARY_PATH=${BABEL_PATH}:${LD_32_LIBRARY_PATH}
	fi
	if [ -z $DYLD_LIBRARY_PATH ]; then
		export DYLD_LIBRARY_PATH=$BABEL_PATH
	else
		export DYLD_LIBRARY_PATH=${BABEL_PATH}:${DYLD_LIBRARY_PATH}
	fi
fi
if [ `$babel_exe -L sdf | grep -c 'not a recognized'` != 0 ]; then
	export BABEL_LIBDIR=$BABEL_PATH/plugins
	export BABEL_DATADIR=$BABEL_PATH/../share/openbabel
fi
$babel_exe -L sdf >& /dev/null
if [ $? != 0 ]; then
	echo "It seems that your OpenBabel installation is not working properly."
	echo "When the following command is issued:"
	echo
	echo "$babel_exe -L sdf >& /dev/null"
	echo
	$babel_exe -L sdf
	abrupt_exit
	exit
fi

if [ -z $1 ]; then
	validation_type=single
else
	if [ $1 = single ]; then
		validation_type=single
	elif [ $1 = multi ]; then
		validation_type=multi
	elif [ $1 = both ]; then
		validation_type="single multi"
	else
		echo "Acceptable validation validation_types are \"single\", \"multi\" and \"both\"."
		abrupt_exit
		exit
	fi
fi
		
for val in $validation_type; do
	echo
	if [ $val = single ]; then
		lines='3,5'
		echo "Single-conformation validation:"
	else
		lines='2,4'
		echo "Multi-conformation validation:"
	fi
	echo "------------------------------"
	for dataset in \
		ace ache bzr cox2 dhfr gpb therm thr; do
		inp=${dataset}/${dataset}_reproduce_orig_alignment_${val}.inp
		out=${dataset}/${dataset}_reproduce_orig_alignment_${val}.out
		if [ -e $out ]; then
			if (grep 'Successful completion' < $out >& /dev/null); then
				continue
			fi
		fi
		$O3A_EXE -i $inp -o $out
	done
	echo
	printf '%-24s%-16s%-16s%-16s\n\n' Dataset Mixed Atom Pharao
	n=1
	while [ $n -le 3 ]; do
		overall[${n}]=0
		let n=n+1
	done
	for dataset in \
		ace ache bzr cox2 dhfr gpb therm thr; do
		out=${dataset}/${dataset}_reproduce_orig_alignment_${val}.out
		error=1
		if [ -e $out ]; then
			if (grep 'Successful completion' < $out >& /dev/null); then
				error=0
			fi
		fi
		if [ $error = 1 ]; then
			echo "Something went wrong during the validation run on the ${i} dataset."
			echo "Please check ${out}, then resubmit."
			abrupt_exit
			exit
		fi
		echo $dataset | tr '[:lower:]' '[:upper:]'
		printf '%-24s%-16s%-16s%-16s\n' 'RMSD (angstrom)' \
			`grep Average < $out | awk '{print $2}' | xargs`
		printf '%-24s' 'Time (HH:MM:SS)'
		n=1
		for elapsed in \
			`grep Elapsed < $out | sed -n ${lines}p \
			| awk '{print $3}' | sed 's/\..*$//'`; do
			rounded_elapsed=`printf %.0f $elapsed`
			let overall[${n}]=overall[${n}]+rounded_elapsed
			hhmmss_elapsed=`printf '%02d:%02d:%02d' \
				$((rounded_elapsed/3600)) \
				$((rounded_elapsed/60%60)) \
				$((rounded_elapsed%60))`
			printf '%-16s' $hhmmss_elapsed
			let n=n+1
		done
		echo
		echo
	done
	echo "Overall time"
	printf '%-24s' '(HH:MM:SS)'
	n=1
	while [ $n -le 3 ]; do
		hhmmss_overall=`printf '%02d:%02d:%02d' \
			$((overall[${n}]/3600)) \
			$((overall[${n}]/60%60)) \
			$((overall[${n}]%60))`
		printf '%-16s' $hhmmss_overall
		let n=n+1
	done
	echo 
	echo
done
echo "Validation succeeded."
echo
