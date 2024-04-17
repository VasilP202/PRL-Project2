#!/bin/bash

# Jméno souboru s definicí hracího pole pro hru Life
life_grid_file=$1;

# Počet iterací, které má algoritmus provést
num_iterations=$2;

# Počet procesorů, které má být použito = počet řádků v hracím poli
num_processors=$(sed -n '$=' "$life_grid_file")

# Přeložení zdrojového kódu
mpic++ --prefix /usr/local/share/OpenMPI -o life life.cpp

# Spuštění programu
mpirun --prefix /usr/local/share/OpenMPI --oversubscribe -np $num_processors life $life_grid_file $num_iterations 		

# Úklid
# rm -f life numbers