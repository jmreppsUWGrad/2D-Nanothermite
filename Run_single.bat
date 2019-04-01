@echo OFF

set fold=Tests\Conserv\AlMoO3\16\NoClipping
REM set fold=Tests\Axisymmetric\AlMoO3\2

REM python main.py Input_file16.txt %fold%
python main.py Input_file_stats.txt %fold%
python Post.py %fold% 0

cd C:\Users\mepps\Documents\Research\2D_Conduction

REM Run_batch.bat