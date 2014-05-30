#!/bin/bash

source parameters.sh

nbase=`wc -l base.info | awk '{print $1}'`
nbaseUse=8
calAb_command=/home/mi/wanghan/study/noneq.msm/committor/template.alanine//tools/analysis/cal.Ab.noAvg
outfileName=base.new.info

nnotUse=`echo "$nbase - $nbaseUse" | bc`
rangeStart=`echo "$nbase - $nbaseUse + 1" | bc`
rangeCode="$rangeStart:$nbase"

if [ ! -f success.dir.name ]; then
    echo "# no simulation is done, exit"
    exit
fi
if [ ! -f calAb.out ] || [ calAb.out -ot success.dir.name ]; then
    echo "# calculating the matrix A and rsh b"
    $calAb_command -n 100 --target-indicate $fht_coreset_target --input-coreset $fht_coreset_data &> calAb.out
    echo "# done"
fi

echo "# solve w by octave"
grep "a = " calAb.out				> solve.w.m
grep "b = " calAb.out				>> solve.w.m
grep "ae = " calAb.out				>> solve.w.m
grep "be = " calAb.out				>> solve.w.m
echo "relative_a_error=ae./(abs(a))"		>> solve.w.m
echo "relative_b_error=be./(abs(b))"		>> solve.w.m
echo "[v d] = eig(a);"				>> solve.w.m
echo "diag(d)"					>> solve.w.m
echo "rd = d($rangeCode,$rangeCode);"		>> solve.w.m
echo "vb = v'*b';"				>> solve.w.m
echo "rvb = vb($rangeCode);"			>> solve.w.m
echo "rrsh = rd^-1 * rvb;"			>> solve.w.m
echo "rsh = [zeros($nnotUse,1); rrsh];"		>> solve.w.m
echo "w = v * rsh;"				>> solve.w.m
echo "fid=fopen('$outfileName', 'w');"		>> solve.w.m
echo "fprintf(fid,'%.10f\n',w);"		>> solve.w.m
echo "fclose(fid);"				>> solve.w.m
echo "# done"

octave solve.w.m

newNbase=`wc -l $outfileName | awk '{print $1}'`
if [ $newNbase -ne $nbase ]; then
    echo "# the number of bases in the output is not equal to the number of lines in base.info, something must be wrong"
    exit
fi
expectedNbase=`echo "$fht_base_phi_number + $fht_base_psi_number" | bc`
if [ $newNbase -ne $expectedNbase ]; then
    echo "# the number of bases in the output is not equal to the number of phi bases ($fht_base_phi_number) plus the number of psi bases ($fht_base_psi_number), exit"
    exit
fi

rm -f base.k.phi.new base.k.psi.new
countLineNewBase=1
for ii in `seq 1 $fht_base_phi_number`;
do
    make_top_line=`grep -v \# $fht_base_phi_k_file | head -n $ii | tail -n 1`
    make_top_type=`echo $make_top_line | awk '{print $1}'`
    make_top_gamma=`echo $make_top_line | awk '{print $2}'`
    make_top_kk_old=`echo $make_top_line | awk '{print $3}'`
    make_top_kk_new=`grep -v \# $outfileName | head -n $countLineNewBase | tail -n 1`
    make_top_kk=`echo "$make_top_kk_old + $make_top_kk_new" | bc -l`
    echo $make_top_type $make_top_gamma $make_top_kk >> base.k.phi.new
    countLineNewBase=$(($countLineNewBase + 1))
done
for ii in `seq 1 $fht_base_psi_number`;
do
    make_top_line=`grep -v \# $fht_base_psi_k_file | head -n $ii | tail -n 1`
    make_top_type=`echo $make_top_line | awk '{print $1}'`
    make_top_gamma=`echo $make_top_line | awk '{print $2}'`
    make_top_kk_old=`echo $make_top_line | awk '{print $3}'`
    make_top_kk_new=`grep -v \# $outfileName | head -n $countLineNewBase | tail -n 1`
    make_top_kk=`echo "$make_top_kk_old + $make_top_kk_new" | bc -l`
    echo $make_top_type $make_top_gamma $make_top_kk >> base.k.psi.new
    countLineNewBase=$(($countLineNewBase + 1))
done


