
#S=-0.5
#T=1.5

bc=$1
eps=1.d-2
eps1=$2
eps2=$3
beta=0.1
N=1000
#eps=1.d-8


#dir='S_'$S'_T_'$T'_eps_'$eps
#mkdir $dir

outfin='SD_b_'$bc'_eps1_'$eps1'_eps2_'$eps2'_epsA_'$eps'.dat'

./calc_SD << EOF
$eps
$eps2
$eps1
$bc
$beta
$N
$outfin
EOF


#mv $outf $dir/.



