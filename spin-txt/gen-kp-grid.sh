#!/bin/bash

###Origen da malha
k0x='0.0'                    
k0y='0.0'                    

### Componente z da malha  
kz='0.0000000000'


echo 'generated mesh' > KPOINTS     
echo '    961' >> KPOINTS            #número total de pontos da malha
echo 'Reciprocal lattice' >> KPOINTS 

for nx in {0..30..1}                 #escolha o numero de pontos na direção a1
do
 x=`echo $nx"/30+$k0x" | bc -l`
 for ny in {0..30..1}                #escolha o numero de pontos na direção a2
 do
  y=`echo $ny"/30+$k0y" | bc -l`
  echo $x"     "$y"     "$kz"    1" >> KPOINTS
 done
done 
