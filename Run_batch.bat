@echo OFF

set proc=6
set fold1=Tests\Axisymmetric\Al
set fold2=Tests\Axisymmetric\CuO
set fold3=Tests\Axisymmetric\comp
set fold4=Tests\Axisymmetric\comp-por

mpiexec -n %proc% python main.py %fold1%\Input_File_axi.txt %fold1%
python Post.py %fold1% 0

mpiexec -n %proc% python main.py %fold2%\Input_File_axi.txt %fold2%
python Post.py %fold2% 0

mpiexec -n %proc% python main.py %fold3%\Input_File_axi.txt %fold3%
python Post.py %fold3% 0

mpiexec -n %proc% python main.py %fold4%\Input_File_axi.txt %fold4%
python Post.py %fold4% 0