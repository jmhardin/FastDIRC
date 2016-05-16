Energy=5
MaxTheta=12
ThetaStep=.5
phi=$1
NRuns=5000

if [ -a tmp_kine_${phi}.csv ] 
then
	rm tmp_kine_${phi}.csv
fi

echo Time Energy Theta Phi pionxadj pionyadj kaonxadj kaonyadj KDERes MidlineUncalibrated MidlineCalibrated LUT >> tmp_kine_${phi}.csv
tail -n 1 tmp_kine_${phi}.csv

for theta in `seq 4 $ThetaStep $MaxTheta`
do
#	../dircfit -n $NRuns -particle_phi $phi -particle_theta $theta -E $Energy &> /dev/null 
#	TmpResKDE=`root -l -q -b "../graphicHistos.C(\"fitdirc.root\",false,${Energy})" | tail -n 1` > /dev/null

	../dircfit -n 0 -line_recon_n $NRuns -particle_phi $phi -particle_theta $theta -E $Energy
	TmpResMidUncalib=`root -l -q -b "../graphicHistos.C(\"fitdirc.root\",false,${Energy})" | tail -n 1` > /dev/null

	CalibString=`../dircfit -n 0 -fill_d_midline_n $NRuns -particle_phi $phi -particle_theta $theta -E $Energy | tail -n 1`		

	../dircfit -n 0 -line_recon_n $NRuns -particle_phi $phi -particle_theta $theta -line_recon_calib $CalibString -E $Energy
	TmpResMidCalib=`root -l -q -b "../graphicHistos.C(\"fitdirc.root\",false,${Energy})" | tail -n 1` &> /dev/null

	../dircfit -n 0 -lut_sim_n $NRuns -particle_phi $phi -particle_theta $theta -E $Energy &> /dev/null 
	TmpResLut=`root -l -q -b "../graphicHistos.C(\"fitdirc.root\",false,${Energy})" | tail -n 1` > /dev/null

	echo `date +%R:%S` $Energy $theta $phi $CalibString $TmpResKDE $TmpResMidUncalib $TmpResMidCalib $TmpResLut >> tmp_kine_${phi}.csv
	
	tail -n 1 tmp_kine_${phi}.csv
		
done

rm fitdirc.root
mv tmp_kine_${phi}.csv res_scan${phi}.csv

