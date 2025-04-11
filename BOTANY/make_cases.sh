set -e

cd .. && source setenv && cd - > /dev/null 2>&1

for f in input_files/profiles_member*txt
do
  ncase=$(echo $f|sed -e "s/^.*prof.*member_//g" -e "s/.txt//g")
  cas=RUN$ncase
  echo $cas
  mkdir -p $cas
  cp $f $cas/profiles_init.txt
  cp input_files/profiles_tendencies.txt $cas
  cd $cas
  python ../REF/driver_DEF.py -n $cas -p
  cd - > /dev/null 2>&1
done
