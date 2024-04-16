##!/bin/bash
#
base_command="python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/collins.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt ../temp/collins"

for i in {0..19}; do
    command="$base_command$i"
    echo $i
    $command
done

#!/bin/bash

#base_command="python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/krogan_core.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt ../temp/krogan_core"
#
#for i in {0..19}; do
#    command="$base_command$i"
#    echo $i
#    $command
#done

#base_command="python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/krogan_extended.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt ../temp/extended"
#
#for i in {0..19}; do
#    command="$base_command$i"
#    echo $i
#    $command
#done

# biogrid
#base_command="python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/biogrid.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt ../temp/biogrid"
#
#for i in {0..19}; do
#    command="$base_command$i"
#    echo $i
#    $command
#done