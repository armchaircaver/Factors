import ctypes

# from http://wolfprojects.altervista.org/articles/dll-in-c-for-python/
mydll = ctypes.cdll.LoadLibrary ("C:\\Users\\Arthur\\source\\repos\\Factors\\x64\\Release\\DLL attempt.dll")
print (mydll)

"""
mydll.ipow.argtypes = [ctypes.c_int, ctypes.c_int]
mydll.ipow.restype = ctypes.c_ulonglong

for i in xrange(2,64):
  p=mydll.ipow(2,i)
  print i, p
  
mydll.ipow2.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_ulonglong) ]
p= ctypes.c_ulonglong()
for i in xrange(2,64):
  mydll.ipow2(2,i, ctypes.byref(p) )
  print i, p.value
"""
for i in range(32):
  print ( mydll.addone(), end="," )

print ("\n")
knuthrand = mydll.knuthrand
knuthrand.restype = ctypes.c_ulonglong
for i in range(32):
  print ( knuthrand(), end="," )
print ("\n")

  
resArray = (ctypes.c_ulonglong * 30)()
print (resArray [0])

res_size = ctypes.c_int()
for i in range(2,64):
  mydll.ipow3(2,i, resArray, ctypes.byref(res_size) )
  print (i, res_size.value,[resArray[i] for i in range(res_size.value)])
