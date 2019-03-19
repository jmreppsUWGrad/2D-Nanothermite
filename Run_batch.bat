@echo OFF

set fold2=Tests\Conserv\AlMoO3\20
python main.py Input_file20.txt %fold2%
python Post-processing.py %fold2%

set fold1=Tests\Conserv\AlMoO3\19
python main.py Input_file19.txt %fold1%
python Post-processing.py %fold1%

REM set fold3=Tests\Conserv\AlMoO3\20
REM python main.py Input_file20.txt %fold3%
REM python Post-processing.py %fold3%
