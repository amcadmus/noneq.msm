dt=1e-3
nst=1e8
branch_feq=10.0
noneq_check_feq=1.0
noneq_time=30.0
warm_time=20.0
rate_lag=1.0
refe_strength=1.0
pert_strength=2.0
resp_order=1
project_name=splitSingleWell
#project_name=tiltDoubleWell

load_saved_resp=no
saved_resp_dir=~/study/noneq.msm/one.test/saved.corr.result5

gamma=1.0
temperature=300.0
x0=-2.0
x1=2.0
seed=`date +%s`

# tilt double well parameters
double_well_k=8.0
double_well_a=1.0
# split single well parameters
single_well_k=8.0
gaussian_sigma=0.16

if echo $project_name | grep tiltDoubleWell &> /dev/null; then
    command_line_param_print="--double-well-k $double_well_k --double-well-a $double_well_a"
    print_k=`printf "%.1f" $double_well_k`
    print_a=`printf "%.1f" $double_well_a`
    save_corr_param_note="k$print_k.a$print_a"
else if echo $project_name | grep splitSingleWell &> /dev/null; then 
    command_line_param_print="--single-well-k $single_well_k --gaussian-sigma $gaussian_sigma"
    print_k=`printf "%.1f" $single_well_k`
    print_s=`printf "%.3f" $gaussian_sigma`
    save_corr_param_note="k$print_k.s$print_s"
fi
fi

