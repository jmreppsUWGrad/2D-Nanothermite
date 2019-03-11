@echo OFF

set fold1=Tests\Conserv\AlMoO3\Mesh\9a
python main.py Input_File1.txt %fold1%
python Post-processing.py %fold1%

set fold2=Tests\Conserv\AlMoO3\16
python main.py Input_File2.txt %fold2%
python Post-processing.py %fold2%

REM set fold3=Tests\Conserv\9
REM python main.py Input_File3.txt %fold3%
REM python Post-processing.py %fold3%
