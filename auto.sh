#!/bin/bash

Z="40"
Zi=$Z
Zf="41"
N="64"
Ni=$N
Nf="64"
#<<comment
while [ $Z -le $Zf ]
do

    if (( $Z % 2 )); then
        sed -i "/proton\_blocking =/ s/proton\_blocking = [0-9]*/proton\_blocking = 1/" hfbtho_NAMELIST.dat
    else
        sed -i "/proton\_number =/ s/proton\_number = [0-9]*/proton\_number = $Z/" hfbtho_NAMELIST.dat
        sed -i "/proton\_blocking =/ s/proton\_blocking = [0-9]*/proton\_blocking = 0/" hfbtho_NAMELIST.dat
    fi

    while [ $N -le $Nf ]
    do
    if (( $N % 2 )); then
        sed -i "/neutron\_blocking =/ s/neutron\_blocking = [0-9]*/neutron\_blocking = 1/" hfbtho_NAMELIST.dat
    else
        sed -i "/neutron\_number =/ s/neutron\_number = [0-9]*/neutron\_number = $N/" hfbtho_NAMELIST.dat
        sed -i "/neutron\_blocking =/ s/neutron\_blocking = [0-9]*/neutron\_blocking = 0/" hfbtho_NAMELIST.dat
    fi
    ./main
    cp thoout.dat Z$Z-N$N.dat
    #sed -n 3p hfbtho_NAMELIST.dat 
    #sed -n 16p hfbtho_NAMELIST.dat
    ##echo "---------------"
    ((N=N+1))
    done
    ((Z=Z+1))
    ((N=$Ni))
done
((Z=$Zi))
#comment
while [ $Z -le $Zf ]
do
    Atom=$(grep "Nucleus:" Z$Z-N$N.dat | awk '{print $2}')
    echo -e "Z:N:BindingEnergy" > $Atom-val_stock
    while [ $N -le $Nf ]
    do
        tEnergy=$(grep -i "tEnergy" Z$Z-N$N.dat | awk '{print $4}')
        sign=-1
        tEnergy=$(echo "$tEnergy * $sign" | bc)
        NrmsRadius=$(grep -i "rms-radius" Z$Z-N$N.dat | awk '{print $3}')
        PrmsRadius=$(grep -i "rms-radius" Z$Z-N$N.dat | awk '{print $4}')
        echo "$Z:$N:$tEnergy" >> $Atom-val_stock
        ((N=N+1))
    done
    column -t -s: $Atom-val_stock > $Atom-Energie_liaison
    ((N=$Ni))
    ((Z=Z+1))
done
((Z=$Zi))
#cat val_stock
cat Energie_liaison
