set -e
mkdir -p logs

sim=SHORT
sim=NONUDGING
sim=REF
sim=RADFRC
sim=ALLFRC
sim=SFCFRC

cd .. && source setenv && cd - > /dev/null 2>&1

for f in input_files/profiles_member*$1.txt
do
  ncase=$(echo $f|sed -e "s/^.*prof.*member_//g" -e "s/.txt//g")
  cas=$sim$ncase
  echo $cas
  test ! -f input_files/surface_fluxes_${ncase}.txt && continue
  mkdir -p $cas
  cp $f $cas/profiles_init.txt
  cp input_files/profiles_tendencies.txt $cas
  case $sim in 
    ALLFRC|RADFRC) cp input_files/profiles_radtend_${ncase}.txt $cas/profiles_radtend.txt ;;
    ALLFRC|SFCFRC) cp input_files/surface_fluxes_${ncase}.txt $cas/surface_flux_forcings.txt ;;
  esac
  cd $cas
  python ../gen_$sim/driver_DEF.py -n $cas -p -v     >> ../logs/log_driver_DEF_$$ 2>&1
  python ../gen_$sim/driver_SCM.py -n $cas -p -c -v  >> ../logs/log_driver_SCM_$$ 2>&1
  cd - > /dev/null 2>&1
done
