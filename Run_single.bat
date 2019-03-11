@echo OFF

set fold=Tests\Conserv\AlMoO3\17

python main.py Input_File_nt.txt %fold%
python Post-processing.py %fold%