#!/bin/bash

source parameters.sh

nbase=`wc -l base.info | awk '{print $1}'`
nbaseUse=6
calAb_command=/home/mi/wanghan/study/noneq.msm/committor/template.alanine//tools/analysis/cal.Ab.noAvg

nnotUse=`echo "$nbase - $nbaseUse" | bc`
rangeStart=`echo "$nbase - $nbaseUse + 1" | bc`
rangeCode="$rangeStart:$nbase"

$calAb_command -n 100 --target-indicate $fht_coreset_target --input-coreset $fht_coreset_data &> calAb.out

grep "a = " calAb.out				> solve.w.m
grep "b = " calAb.out				>> solve.w.m
echo "[v d] = eig(a);"				>> solve.w.m
echo "d"					>> solve.w.m
echo "rd = d($rangeCode,$rangeCode);"		>> solve.w.m
echo "vb = v'*b'"				>> solve.w.m
echo "rvb = vb($rangeCode)"			>> solve.w.m
echo "rrsh = rd^-1 * rvb"			>> solve.w.m
echo "rsh = [zeros($nnotUse,1); rrsh]"		>> solve.w.m
echo "w = v * rsh"				>> solve.w.m




rm -f calAb.out
