import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import csv

def spline1(x,y,point):
    f = interpolate.interp1d(x,y,kind="cubic")
    X = np.linspace(x[0],x[-1],num=point,endpoint=True)
    Y = f(X)
    return X,Y

with open('curves.csv') as fin:
    s = csv.reader(fin)
    cnt = 0
    cx = []
    cy1 = []
    cy2 = []
    cy3 = []
    cy4 = []
    for row in s:
        if (cnt>0):
            cx += [cnt]
            cy1 += [1+float(row[1])/100]
            cy2 += [1+float(row[5])/100]
            cy3 += [1+float(row[6])/100]
            cy4 += [1+float(row[7])/2]
        print(cnt,row[1],row[5],row[6],row[7])
        cnt += 1


plt.title('community')
plt.plot(cx,cy1)
plt.show()
plt.title('workplace')
plt.plot(cx,cy2)
plt.show()
plt.title('household')
plt.plot(cx,cy3)
plt.show()
plt.title('obon')
plt.plot(cx,cy4)
plt.show()

x = [1,10,20,30,40,50,127,134,141,148,155,162,169,176,183,190,197,204,211,218,225,232,239,270,280,287]
y = [0,0,0,0,0,0,3,5,7,11,21,33,44,67,79,85,89,90,91,91,92,93,94,95,95,95]
for i in range(len(x)):
    y[i] = float(y[i])/100*0.5 + 1.0
a,b = spline1(x,y,287)

plt.title('delta')
plt.plot(a,b)
plt.show()
