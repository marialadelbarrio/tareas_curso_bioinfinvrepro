#!/bin/bash

#Este script adapta el tutorial de  http://biocontainers.pro/docs/101/running-example/

#Vamos a adaptar el ejemplo para trabajar con la  beta-caseína

#bajando la secuencia FASTA de la caseína humana desde UniProt

wget https://www.uniprot.org/uniprot/P05814.fasta

#Vamos a compararlo con las secuencias de referencia de las vacas, empezamos bajando la base de datos de NCBI

curl -O ftp://ftp.ncbi.nih.gov/refseq/B_taurus/mRNA_Prot/cow.1.protein.faa.gz 

#y para descomprimir

gunzip cow.1.protein.faa.gz

#preparación de la base de datos, montando un volumen para hacer visibles los resultados blast fuera del contenedor. Se usa la ruta absoluta para esto.

sudo docker run -v /home/maria/Escritorio/bioinfvrepro/ejercicios_tareas:/data/ biocontainers/blast makeblastdb -in cow.1.protein.faa -dbtype prot

#ahora se hacen los alineamientos finales

sudo docker run -v /home/maria/Escritorio/bioinfvrepro/ejercicios_tareas:/data/ biocontainers/blast blastp -query P05814.fasta -db cow.1.protein.faa -out resultadoscaseina.txt



