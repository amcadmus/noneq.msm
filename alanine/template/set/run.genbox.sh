#!/bin/bash

if test $# -lt 1; then
    echo "tell me the box size by using:"
    echo "./run.genbox.sh boxsize"
    exit
fi
if test ! -d gromacs.conf; then
    echo "no dir gromacs.conf, exit"
    exit
fi
if test ! -d gromacs.top; then
    echo "no dir gromacs.top, exit"
    exit
fi

box=$1
p_box=`printf %.2f $box`
i_box=`echo $p_box | cut -d '.' -f 1`
d_box=`echo $p_box | cut -d '.' -f 2`
p_box=`printf %03d.%02d $i_box $d_box`
box=`printf %d.%02d $i_box $d_box`

target_dir=gromacs.conf.box.$p_box

if test -d $target_dir; then
    echo "find dir $target_dir, do nothing"
    exit
fi

mkdir -p $target_dir

cp gromacs.conf/protein.gro $target_dir
cp gromacs.conf/spc216.gro $target_dir
cp gromacs.top/topol.top $target_dir
cp gromacs.top/grompp.mdp $target_dir

echo "# generate box"
cd $target_dir
sed -e "s/2\.5 2\.5 2\.5/$box $box $box/g" protein.gro > tmp.tmp
mv -f tmp.tmp protein.gro
grompp -c protein.gro &> genbox.log
echo 3 0 | trjconv -f protein.gro -o out.gro -center -pbc whole &> genbox.log
genbox -cp out.gro -cs spc216.gro -o conf.gro &> genbox.log
nsol=`grep "Number of SOL molecules:" genbox.log | tail -n 1 | cut -d ":" -f 2`
sed -e "s/SOL.*/SOL $nsol/g" topol.top > tmp.tmp
mv -f tmp.tmp topol.top

echo "# do NPT run"
sed -e "/^Pcoupl/s/= no/= Parrinello-Rahman/g" grompp.mdp |
sed -e "/^Tcoupl/s/= no/= v-rescale/g" |
sed -e "/^freezegrps/s/= .*/= FIX/g" |
sed -e "/^freezedim/s/= .*/= Y Y Y/g" |
sed -e "/^nsteps/s/= .*/= 100000/g" > tmp.tmp
mv -f tmp.tmp grompp.mdp
echo "a 9" > command.index
echo "name 17 FIX" >> command.index
echo "q" >> command.index
cat command.index | make_ndx -f conf.gro &> genbox.log
rm -f command.index
grompp -n index.ndx &> genbox.log
mdrun -v &> genbox.log

echo "# resize box"
echo 18 | g_energy -b 50 -xvg none &> genbox.log
newbox=`avg_jk -v col=2 energy.xvg | grep -v \# | awk '{print $1}'`
oldbox=`tail -n 1 confout.gro | awk '{print $1}'`
scale=`echo "$newbox / $oldbox" | bc -l`
editconf -scale $scale $scale $scale -f confout.gro -o conf.gro &> genbox.log

echo "# do NVT"
sed -e "/^Pcoupl/s/= Parrinello-Rahman/= no/g" grompp.mdp |
sed -e "/^nsteps/s/= .*/= 100000/g" > tmp.tmp
mv -f tmp.tmp grompp.mdp
grompp -n index.ndx &> genbox.log
mdrun -v &> genbox.log

echo "# cleaning"
mv -f confout.gro conf.gro
grocleanit
rm -f out.gro energy.xvg spc216.gro
