set -e

ncas=006
cas=RUN$ncas

cp ../input_files/profiles_member_$ncas.txt profiles_init.txt
cp ../input_files/profiles_tendencies.txt .

python ../REF/driver_DEF.py -n $cas -p
