@echo OFF

set fold=Tests\Conserv\AlMoO3\3

python main.py Input_File_nt.txt %fold%
python Post-processing.py %fold%