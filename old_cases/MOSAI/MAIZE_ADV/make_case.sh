cas=MAIZE
nam=_ADV

# avec les flux de surface observés (LE,H,ustar) et tskin
# avec des tendances représentant l'advection de grande échelle, prises dans ARPEGE oper
# pas de rayonnement ni de forçage radiatif

# => j'ai essayé de faire tourner avec MNH 1D mais l'outil DEF->MNH a pas
# fonctionné car ne sait pas convertir les tendances de T en theta (nécessaire
# à MNH). Je vais m'y remettre mais du coup je ne sais pas si le cas tourne ou
# pas comme ça => je git quand même en attendant et je vais prévenir Alice &
# Emilie pour qu'elles testent avec LMDZ et ARPEGE.

cp ../input_files/CR202308190452.cor.nc input.nc
cp ../forcing_files/flux_MOSAI_$cas.txt flux.txt
cp ../forcing_files/Arpege-oper-L105_Lannemezan_202308190500-202308201800.nc tendencies.nc

python ../REF/driver_DEF.py -c $cas -n $nam -a
python ../REF/driver_SCM.py -c $cas -n $nam 
