@echo OFF

set fold=Tests\5\1

python main.py Input_File_nt.txt %fold%
python Post-processing.py %fold%