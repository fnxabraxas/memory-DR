
st1=$1
st2=$2
di=1

outf='out_payM_'$st1'_'$st2'.dat'

./payM << EOF
$st1
$st2
$di
1.d-2
$outf
EOF
