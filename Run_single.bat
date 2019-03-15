@echo OFF

set fold=Tests\Conserv\AlCuO\2

python main.py Input_File_nt.txt %fold%
python Post-processing.py %fold%

cd C:\Users\mepps\Documents\Research\2D_Conduction

REM Run_batch.bat