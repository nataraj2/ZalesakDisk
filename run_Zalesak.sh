rm -rf out
mpif90 -c ZalesakIC.f90
mpif90 -o out ZalesakIC.o
./out

