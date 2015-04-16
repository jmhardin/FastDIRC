arg="-mirror_angle_change_yunc"
typerun=ang440_threeseg_${arg}
numruns=10000

mkdir ${typerun}

amount=0.01
time ./dircfit -n ${numruns} ${arg} ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.02
time ./dircfit -n ${numruns} ${arg} ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.04
time ./dircfit -n ${numruns} ${arg} ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.08
time ./dircfit -n ${numruns} ${arg} ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.16
time ./dircfit -n ${numruns} ${arg} ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.3
time ./dircfit -n ${numruns} ${arg} ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.6
time ./dircfit -n ${numruns} ${arg} ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1.2
time ./dircfit -n ${numruns} ${arg} ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1.6
time ./dircfit -n ${numruns} ${arg} ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=2.5
time ./dircfit -n ${numruns} ${arg} ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

