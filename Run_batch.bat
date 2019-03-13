@echo OFF

set fold1=Tests\Conserv\AlMoO3\18
python main.py Input_File18.txt %fold1%
python Post-processing.py %fold1%

set fold2=Tests\Conserv\AlMoO3\19
python main.py Input_File19.txt %fold2%
python Post-processing.py %fold2%

set fold3=Tests\Conserv\AlMoO3\20
python main.py Input_File20.txt %fold3%
python Post-processing.py %fold3%
