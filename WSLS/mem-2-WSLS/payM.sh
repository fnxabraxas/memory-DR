
st1=0
st2=511
di=0

outf='out_payM_diag.dat'

./payM << EOF
$st1
$st2
$di
1.d-2
$outf
EOF

st1=0
st2=511
di=1

outf='out_payM_'$st1'_'$st2'.dat'

./payM << EOF
$st1
$st2
$di
1.d-2
$outf
EOF
