@echo OFF

call conda activate base
set proc=6
set fold1=Tests\5
REM set fold1=Tests\Axisymmetric\Al
set fold2=Tests\Axisymmetric\CuO
set fold3=Tests\Axisymmetric\comp
set fold4=Tests\Axisymmetric\comp-por

mpiexec -n %proc% python main.py %fold1%\Start.txt %fold1%
REM python Post.py %fold1% 0

REM mpiexec -n %proc% python main.py %fold2%\Start.txt %fold2%
REM python Post.py %fold2% 0

REM mpiexec -n %proc% python main.py %fold3%\Start.txt %fold3%
REM python Post.py %fold3% 0

REM mpiexec -n %proc% python main.py %fold4%\Start.txt %fold4%
REM python Post.py %fold4% 0

call conda deactivate