typerun=oil_flat_res_perp_water

mkdir ${typerun}

amount=6.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=10.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=15.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=20.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=25.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=30.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=35.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=40.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=50.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=60.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=70.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

amount=80.0
time ./dircfit ${amount}
mv fitdirc.root ${typerun}/fitdirc_${typerun}_${amount}.root

