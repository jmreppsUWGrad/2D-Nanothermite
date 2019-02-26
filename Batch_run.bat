@echo OFF

set fold1=Tests\Conserv\7
python main.py Input_File1.txt %fold1%
python Post-processing.py %fold1%

set fold2=Tests\Conserv\8
python main.py Input_File2.txt %fold2%
python Post-processing.py %fold2%

set fold3=Tests\Conserv\9
python main.py Input_File3.txt %fold3%
python Post-processing.py %fold3%
