@echo OFF

set fold=Tests\Axisymmetric\AlMoO3\18

python main.py Input_File18.txt %fold%
python Post-processing.py %fold% 0

cd C:\Users\mepps\Documents\Research\2D_Conduction

REM Run_batch.bat