rm -rf out
gfortran-mp-11 -c ZalesakIC.f90
gfortran-mp-11 -o out ZalesakIC.o
./out

