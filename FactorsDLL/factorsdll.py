import ctypes
from os import getcwd


cwd = getcwd()
Factorsdll = ctypes.cdll.LoadLibrary (cwd+"\FactorsDLL.dll")

# array to receive prime factors, large enough for factors of any 64 bit number
resArray = (ctypes.c_ulonglong * 64)() 

res_size = ctypes.c_int()
number_to_factor = ctypes.c_ulonglong()
factorfound = ctypes.c_ulonglong()

# suggested by http://tungwaiyip.info/blog/2009/07/16/ctype_performance_benchmark
Factorsdll_factor = Factorsdll.factor

def factoriseCPP(n):
  number_to_factor.value=n
  Factorsdll_factor(number_to_factor, resArray, ctypes.byref(res_size) )
  return resArray[0:res_size.value]

Factorsdll.isprime.restype = ctypes.c_bool
def isprimeCPP(n):
  number_to_factor.value=n
  return Factorsdll.isprime(number_to_factor, ctypes.byref(factorfound) )

def isprimefCPP(n, factor):
  number_to_factor.value=n
  result =  Factorsdll.isprime(number_to_factor, ctypes.byref(factorfound) )
  factor[0] = factorfound
  return result
  
def nulltestCPP(n):
  number_to_factor.value=n
  Factorsdll.nulltest(number_to_factor, resArray, ctypes.byref(res_size) )
  return resArray[0:res_size.value]

def nothing(n):
  number_to_factor.value=n
  Factorsdll.nothing(number_to_factor, resArray, ctypes.byref(res_size) )
  return []

Factorsdll.GetMethod.restype = ctypes.c_char
def getmethodCPP():
  return Factorsdll.GetMethod()

#---------------- all factors ---------------------
"""
from http://wwwhomes.uni-bielefeld.de/achim/highly.txt
169th HCN is  18401055938125660800 ,<2^64,  with 184320 factors
prime powers 7 4 2 2 1 1 1 1 1 1 1 1 1

170th  HCN is 27601583907188491200 ,> 2^64, with 193536 factors
prime powers 6 5 2 2 1 1 1 1 1 1 1 1 1
"""

allfactorsArray = (ctypes.c_ulonglong * 193536)() # long enough for a 64 bit number

Factorsdll_allfactor = Factorsdll.allfactor

def allfactorCPP(n):
  number_to_factor.value=n
  Factorsdll_allfactor(number_to_factor, allfactorsArray, ctypes.byref(res_size) )
  return allfactorsArray[0:res_size.value]

Factorsdll.SquareFreePart.restype = ctypes.c_ulonglong
def SquareFreePartCPP(n):
  number_to_factor.value=n
  return Factorsdll.SquareFreePart(number_to_factor)

Factorsdll.Totient.restype = ctypes.c_ulonglong
def TotientCPP(n):
  number_to_factor.value=n
  return Factorsdll.Totient(number_to_factor)
 

if __name__ == "__main__":
  for n in (6,101, 2**63-1, 2**64-1):
    print( n,"factorises into", factoriseCPP(n) )
