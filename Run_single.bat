@echo OFF

call conda activate base
set proc=6
set fold=Tests\Conserv\AlMoO3\16\NoClipping
REM set fold=Tests\Axisymmetric\AlMoO3\2

REM mpiexec -n %proc% python main.py Input_file16.txt %fold%
mpiexec -n %proc% python main.py Input_file_stats.txt %fold%
python Post.py %fold% 0

cd C:\Users\mepps\Documents\Research\2D_Conduction

REM Run_batch.bat

call conda deactivate