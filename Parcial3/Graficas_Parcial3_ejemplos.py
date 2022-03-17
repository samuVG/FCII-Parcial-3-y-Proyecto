import numpy as np
import matplotlib.pyplot as plt

datos = []

with open("FinDiff_a=1.000000_b=3.000000.txt") as fname:
	for lineas in fname:
		datos.append(lineas.split())

Datos=np.array(datos[2:])

X=Datos[:,0].astype(float)
W=Datos[:,1].astype(float)
Y=Datos[:,2].astype(float)

plt.figure(figsize=(6,6))
plt.plot(X,W,"b",label="Solucion numerica")
plt.plot(X,Y,"r",label="Solucion exacta: y = x^2 + 16/x")
plt.title("y'' = (32 + 2*x^3 - y*y') /8 ")
plt.xlabel("x")
plt.ylabel("Y(x)")
plt.legend()
plt.savefig("Ejemplo1_11.4_Burden.jpg")
plt.show()

datos = []

with open("FinDiff_a=1.000000_b=2.000000.txt") as fname:
	for lineas in fname:
		datos.append(lineas.split())

Datos=np.array(datos[2:])

X=Datos[:,0].astype(float)
W=Datos[:,1].astype(float)
Y=Datos[:,2].astype(float)

plt.figure(figsize=(6,6))
plt.plot(X,W,"b",label="Solucion numerica")
plt.plot(X,Y,"r",label="Solucion exacta: y = (x+1)^-1")
plt.title("y'' = y^3 - y*y' ")
plt.xlabel("x")
plt.ylabel("Y(x)")
plt.legend()
plt.savefig("Exercise4a_11.4_Burden.jpg")
plt.show()


datos = []

with open("FinDiff_a=0.000000_b=5.000000.txt") as fname:
	for lineas in fname:
		datos.append(lineas.split())

Datos=np.array(datos[2:])

X=Datos[:,0].astype(float)
W=Datos[:,1].astype(float)
Y=Datos[:,2].astype(float)

plt.figure(figsize=(6,6))
plt.plot(X,W,"b",label="Solucion numerica")
plt.plot(X,Y,"r",label="Solucion exacta: y = 1 /sqrt( 1 + x^2 /3 )")
plt.title("y'' = -( y^5 + 2*y'/x ) ")
plt.xlabel("x")
plt.ylabel("Y(x)")
plt.legend()
plt.savefig("Aplicacion_Fisica.jpg")
plt.show()

