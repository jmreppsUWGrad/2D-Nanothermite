@echo OFF

set fold1=Tests\Conserv\4
python main.py Input_File1.txt %fold1%
python Post-processing.py %fold1%

set fold2=Tests\Conserv\5
python main.py Input_File2.txt %fold2%
python Post-processing.py %fold2%

set fold3=Tests\Conserv\6
python main.py Input_File3.txt %fold3%
python Post-processing.py %fold3%
