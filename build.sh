#gfortran heat_polar.f90 -o heat_polar -lopenblas -lpthread -lm
gfortran heat_polar.f90 -o heat_polar -L/usr/lib -llapack -L/usr/lib -lblas -lm
