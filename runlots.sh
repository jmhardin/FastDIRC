typerun=ang440_3seg_water_5000MeV_focy_ang

mkdir ${typerun}

amount=.01
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=.04
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=.1
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=.2
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

