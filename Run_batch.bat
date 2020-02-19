@echo OFF

call conda activate base
set proc=6
set fold1=Tests\5
set fold1=nano\dc\1um\1
set fold2=nano\dc\1um\2
set fold3=nano\dc\1um\3
set fold4=nano\dc\1um\4
set fold5=nano\dc\1um\5
set fold6=nano\base\40nm\6

mpiexec -n %proc% python main.py %fold1%\Start.txt %fold1%
mpiexec -n %proc% python main.py %fold2%\Start.txt %fold2%
mpiexec -n %proc% python main.py %fold3%\Start.txt %fold3%
mpiexec -n %proc% python main.py %fold4%\Start.txt %fold4%
mpiexec -n %proc% python main.py %fold5%\Start.txt %fold5%
REM mpiexec -n %proc% python main.py %fold6%\Start.txt %fold6%

call conda deactivate

SHUTDOWN /l