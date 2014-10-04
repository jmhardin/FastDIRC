typerun=ang440_3seg_focmirror_unc

mkdir ${typerun}

amount=.01
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=.02
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=.04
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=.08
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=.16
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=.32
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

