@echo OFF

set fold=Tests\Conserv\AlMoO3\16

python main.py Input_file16.txt %fold%
python Post-processing.py %fold%

cd C:\Users\mepps\Documents\Research\2D_Conduction

REM Run_batch.bat