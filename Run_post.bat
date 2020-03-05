@echo OFF

call conda activate base
set fold1=nano\1_base
set fold2=nano\lambda
set fold3=nano\Ea
set fold4=nano\A0
set fold5=nano\dc
set fold6=nano\gas

python Post.py %fold1%\Post_Input.txt
python Post.py %fold2%\Post_Input.txt
REM python Post.py %fold3%\Post_Input.txt
REM python Post.py %fold4%\Post_Input.txt
python Post.py %fold5%\Post_Input.txt
REM python Post.py %fold6%\Post_Input.txt

call conda deactivate