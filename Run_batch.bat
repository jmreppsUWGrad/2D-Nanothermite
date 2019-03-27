@echo OFF

set fold1=Tests\Axisymmetric\AlMoO3\2
set fold2=Tests\Axisymmetric\AlMoO3\8
set fold3=Tests\Axisymmetric\AlMoO3\9

python main.py Input_file2.txt %fold1%
python Post-processing.py %fold1% 1

python main.py Input_file8.txt %fold2%
python Post-processing.py %fold2% 1

python main.py Input_file9.txt %fold3%
python Post-processing.py %fold3% 1
