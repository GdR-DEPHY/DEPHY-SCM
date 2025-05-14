cas=MAIZE
nam=_H0

# avec les flux de surface observés (LE,H,ustar) et tskin
# mais le flux sensible est seuillé à 0 (pas de valeur négative)
# pas de tendances de grande échelle
# pas de rayonnement ni de forçage radiatif

cp ../input_files/CR202308190452.cor.nc input.nc
cp ../forcing_files/flux_MOSAI_$cas.txt flux.txt

python ../REF/driver_DEF.py -c $cas -n $nam -s 0
python ../REF/driver_SCM.py -c $cas -n $nam
