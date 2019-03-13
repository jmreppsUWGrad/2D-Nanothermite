@echo OFF

set fold=Tests\Conserv\AlMoO3\Mesh\8a

python main.py Input_File8a.txt %fold%
python Post-processing.py %fold%

cd C:\Users\mepps\Documents\Research\2D_Conduction

Run_batch.bat