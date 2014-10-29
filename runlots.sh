typerun=ang440_threeseg_mainmirror_unc

mkdir ${typerun}

amount=0.01
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.02
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.04
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.07
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.1
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.2
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.4
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=0.7
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=2
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

