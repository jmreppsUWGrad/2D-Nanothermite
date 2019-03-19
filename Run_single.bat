@echo OFF

set fold=Tests\Conserv\AlMoO3\18

python main.py Input_file18.txt %fold%
python Post-processing.py %fold%

cd C:\Users\mepps\Documents\Research\2D_Conduction

REM Run_batch.bat