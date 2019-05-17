@echo OFF

set proc=6
set fold1=Tests\Axisymmetric\AlMoO3\13
set fold2=Tests\Axisymmetric\AlMoO3\14
set fold3=Tests\Axisymmetric\AlMoO3\15

REM mpiexec -n %proc% python main.py Input_file13.txt %fold1%
python Post.py %fold1% 0

REM mpiexec -n %proc% python main.py Input_file14.txt %fold2%
python Post.py %fold2% 0

REM mpiexec -n %proc% python main.py Input_file15.txt %fold3%
python Post.py %fold3% 0
