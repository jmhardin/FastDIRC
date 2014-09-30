typerun=ang440_3seg_index_5000MeV

mkdir ${typerun}

amount=1.47
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1.45
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1.43
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1.41
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1.39
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1.37
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1.35
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1.33
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=1.31
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

