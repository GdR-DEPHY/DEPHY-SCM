### options
cover=MAIZE         # type of surface cover: MAIZE or DECIDUOUS

# initial profile:
initpf=rs_idea      # rs_idea (idealised profile from radiosoundings and surface observations)
                    # or rs_smooth (profile from smoothed radiosoundings)

# large scale advection (taken for ARPEGE-oper or ERA5)
advTq=True
advuv=True
adv_from=ARPEGEoper # advection from ARPEGEoper or ERA5
smooth_adv=0.5      # zero means no smooth
zmax_adv=5000.      # advection is imposed up to this altitude (m)

# geostrophic wind
geo=True            # geostrophic wind taken from VHF radar in P2OA
smooth_geo=0.5      # zero means no smooth
zmax_geo=5000.      # geostrophic wind is imposed up to this altitude (m)

# radiation
rad=True            # to activate radiation

# surface forcing
forc_flux=False     # forcing from observations of sensible and latent heat fluxes
forc_ts=False       # forcing from observations of surface temperature


### get initial profiles and forcings files
# initial profile
if [ $initpf = 'rs_idea' ] || [ $initpf = 'rs_idea_ql' ]
then
  cp ../input_files/RS_202308190500_${initpf}.csv input.csv
elif [ $initpf == 'rs_smooth' ]
then
  cp ../input_files/CR202308190452.cor.nc input.nc
fi

# large scale advection
if [ $advTq == 'True' ] || [ $advuv == 'True' ]
then
  if [ $adv_from == 'ARPEGEoper' ]
  then
    cp ../forcing_files/Arpege-oper-L105_Lannemezan_202308190500-202308201800.nc tendencies.nc
  elif [ $adv_from == 'ERA5' ]
  then
    cp ../forcing_files/advection_ERA5.nc tendencies.nc
  fi
fi

# geostrophic wind
if [ $geo == 'True' ]
then
  cp ../forcing_files/P2OA-CRA-Lannemezan_VHF-RADAR_L2B_2023-08-19-20.nc geostrophic_wind.nc
fi

# surface forcing
if [ $forc_flux == 'True' ] || [ $forc_ts == 'True' ]
then
  cp ../forcing_files/flux_MOSAI_${cover}.txt flux.txt
fi

# run scripts to generate nc files
python ../REF/driver_DEF.py $cover $initpf --advTq $advTq --advuv $advuv $adv_from $smooth_adv $zmax_adv --geo $geo $smooth_geo $zmax_geo --rad $rad --ffx $forc_flux --fts $forc_ts
python ../REF/driver_SCM.py $cover $initpf --advTq $advTq --advuv $advuv $adv_from $smooth_adv $zmax_adv --geo $geo $smooth_geo $zmax_geo --rad $rad --ffx $forc_flux --fts $forc_ts
