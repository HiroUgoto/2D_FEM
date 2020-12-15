rem Use in src folder
cd %~dp0input
python mk_mesh.py
timeout /t 5 /nobreak
cd %~dp0
python main.py
timeout /t 30 /nobreak
