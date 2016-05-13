Energy=5
MaxTheta=12
ThetaStep=2
MaxPhi=180
PhiStep=10
NRuns=500

if [ -a tmp_kine.csv ] 
then
	rm tmp_kine.csv
fi

echo Time Energy Theta Phi pionxadj pionyadj kaonxadj kaonyadj KDERes MidlineUncalibrated MidlineCalibrated LUT >> tmp_kine.csv
tail -n 1 tmp_kine.csv

for theta in `seq 0 $ThetaStep $MaxTheta`
do
	for phi in `seq 0 $PhiStep $MaxPhi`
	do
		../dircfit -n $NRuns -particle_phi $phi -particle_theta $theta -E $Energy &> /dev/null 
		TmpResKDE=`root -l -q -b "../graphicHistos.C(\"fitdirc.root\",false,${Energy})" | tail -n 1` > /dev/null

		../dircfit -n 0 -line_recon_n $NRuns -particle_phi $phi -particle_theta $theta -E $Energy &> /dev/null 
		TmpResMidUncalib=`root -l -q -b "../graphicHistos.C(\"fitdirc.root\",false,${Energy})" | tail -n 1` > /dev/null

		CalibString=`../dircfit -n 0 -fill_d_midline_n $NRuns -particle_phi $phi -particle_theta $theta -E $Energy | tail -n 1` &> /dev/null 		

		../dircfit -n 0 -line_recon_n $NRuns -particle_phi $phi -particle_theta $theta -line_recon_calib $CalibString -E $Energy &> /dev/null
		TmpResMidCalib=`root -l -q -b "../graphicHistos.C(\"fitdirc.root\",false,${Energy})" | tail -n 1` &> /dev/null

		../dircfit -n 0 -lut_sim_n $NRuns -particle_phi $phi -particle_theta $theta -E $Energy &> /dev/null 
		TmpResLut=`root -l -q -b "../graphicHistos.C(\"fitdirc.root\",false,${Energy})" | tail -n 1` > /dev/null

		echo `date +%R:%S` $Energy $theta $phi $CalibString $TmpResKDE $TmpResMidUncalib $TmpResMidCalib $TmpResLut >> tmp_kine.csv
		
		tail -n 1 tmp_kine.csv
		
	done
done
rm fitdirc.root
mv tmp_kine.csv res_scan.csv

