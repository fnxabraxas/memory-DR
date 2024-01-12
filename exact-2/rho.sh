# rho: calculates the fixation probabilities for the same elements that those that are listed in the payoff matrix files

b=2.
beta=1.
N=100

./rho << EOF
out_payM_diag.dat
out_rho_b2_beta1.dat
$b
$beta
$N
3
out_payM_0_15000.dat
out_payM_15000_40000.dat
out_payM_40000_65535.dat
EOF
