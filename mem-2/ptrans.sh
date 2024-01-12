
S=0.5
T=1.5

N=100
beta=0.1d0

inpf='out_payM_0_15000.dat'
inpd='out_payM_diag.dat'

outf='ptran_ttt_'$S'_'$T'.dat'

./ptrans << EOF
$S
$T
$N
$beta
$inpd
$inpf
$outf
EOF


