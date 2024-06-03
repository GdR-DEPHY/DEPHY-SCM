ncks -O --mk_rec_dmn time Arpege-oper-L105_Lannemezan_2023081900.nc Arpege-oper-L105_Lannemezan_2023081900.nc 
ncks -O --mk_rec_dmn time Arpege-oper-L105_Lannemezan_2023082000.nc Arpege-oper-L105_Lannemezan_2023082000.nc 
ncrcat Arpege-oper-L105_Lannemezan_2023081900.nc Arpege-oper-L105_Lannemezan_2023082000.nc Arpege-oper-L105_Lannemezan_20230819-20.nc 
ncks -d time,5,42 -d time_f,5,41 Arpege-oper-L105_Lannemezan_20230819-20.nc -O Arpege-oper-L105_Lannemezan_202308190500-202308201800.nc 
