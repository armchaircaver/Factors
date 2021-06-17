
from factorsdll import factoriseCPP

for n in (6,101, 2**63-1, 2**64-1):
  print( n,"factorises into", factoriseCPP(n) )
