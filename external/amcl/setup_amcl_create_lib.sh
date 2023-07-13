cd amcl/version3/c
python ../../../config64_build.py
rm -rf /usr/local/include/amcl
mkdir /usr/local/include/amcl
cp *.h /usr/local/include/amcl
rm /usr/local/lib/libamcl.a
cp amcl.a /usr/local/lib/libamcl.a
cd ../../..