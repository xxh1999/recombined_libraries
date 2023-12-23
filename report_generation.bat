@echo off
for %%f in (covint_complex\*.pdb) do (
  echo Processing file: %%f
  python plipcmd.py -f %%f --name %%~nf -o covint_report -t
)
