set -e
sim=SHORT
sim=REF

rm -f RUN*

cd .. && source setenv && cd - > /dev/null 2>&1

for f in input_files/profiles_member*$1*txt
do
  ncase=$(echo $f|sed -e "s/^.*prof.*member_//g" -e "s/.txt//g")
  cas=RUN$ncase
  echo $cas
  mkdir -p $cas
  cp $f $cas/profiles_init.txt
  cp input_files/profiles_tendencies.txt $cas
  cd $cas
  python ../$sim/driver_DEF.py -n $cas -p
  cd - > /dev/null 2>&1
done

case $sim in 
  SHORT) dir=1day;;
  REF) dir=3days;;
esac

test -d $dir && mv $dir /tmp/${dir}_$$

mkdir -p $dir
mv -f RUN* $dir
ln -sf $dir/RUN* .
