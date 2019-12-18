@echo OFF

call conda activate base
set proc=6
set fold1=Tests\5
set fold1=nano\axi\1um\1
set fold2=nano\axi\1um\2
set fold3=nano\axi\1um\3
set fold4=nano\axi\1um\4
set fold5=nano\axi\1um\5
set fold6=nano\axi\1um\6

mpiexec -n %proc% python main.py %fold1%\Start.txt %fold1%
REM python Post.py %fold1%

REM mpiexec -n %proc% python main.py %fold2%\Start.txt %fold2%
REM python Post.py %fold2%

mpiexec -n %proc% python main.py %fold3%\Start.txt %fold3%
REM python Post.py %fold3%

mpiexec -n %proc% python main.py %fold4%\Start.txt %fold4%
REM python Post.py %fold4%

mpiexec -n %proc% python main.py %fold5%\Start.txt %fold5%
REM python Post.py %fold5%

REM mpiexec -n %proc% python main.py %fold6%\Start.txt %fold6%
REM python Post.py %fold6%

call conda deactivate