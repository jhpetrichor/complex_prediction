# usage
- 依赖于 python 2.7
```
conda activate evaluation

python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/collins.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt /home/jh/code/complex_predict/result/CFinder/cfinder_collins.txt
python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/gavin.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt /home/jh/code/complex_predict/result/CFinder/cfinder_gavin.txt
python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/krogan_core.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt /home/jh/code/complex_predict/result/CFinder/cfinder_krogan_core.txt
python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/krogan_extended.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt /home/jh/code/complex_predict/result/CFinder/cfinder_krogan_extended.txt
``
`
python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/collins.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt ../temp/collins.txt
python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/gavin.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt ../temp/gavin.txt
python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/krogan_core.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt ../temp/krogan_core.txt 
python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/krogan_extended.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt ../temp/krogan_extended.txt
python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/biogrid.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt ../temp/biogrid.txt
````
``
#!/bin/bash
base_command="python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/collins.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt ../temp/collins"
for i in {1..19}; do
command="$base_command$i"
$command
done
``

python match.py -n /home/jh/code/complex_predict/dataset/Yeast/PPI/biogrid.txt /home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt /home/jh/code/complex_predict/result/matched/matched_biogrid.txt
```