@echo OFF

call conda activate base

cd nano\Post_all

for %%x in (dir *.txt) do (python ..\..\Post.py %%~fx)
cd ..\..

REM python Post.py nano\Post_base.txt
REM python Post.py nano\Post_dc.txt
REM python Post.py nano\Post_Ea.txt
REM python Post.py nano\Post_lambda.txt
REM python Post.py nano\Post_mu.txt
REM python Post.py nano\Post_gas.txt
REM python Post.py nano\Post_A0.txt

call conda deactivate