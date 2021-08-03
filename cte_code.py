import numpy as np
import os
import matplotlib.pyplot as plt
import re

filename = str(input('Enter file name: ')) 

with open('D:\\' + filename + '.txt','r') as f:
    x = []
    y = []
    dy = []
    for line in f:
        try:
            line = line.strip()
            x_value, y_value, dy_value = map(float, re.split('\s+', line))
        except ValueError as err:
            continue
        x.append(x_value)
        y.append(y_value)
        dy.append(dy_value)
        
npoint = len(x)
lat = float(input('Lattice constant: '))
v = np.empty((npoint,7))
ee = 2

#SETUPQ
'''Определяет матрицы Q_TRANSP*D и Q_TRANSP*D**2*Q из X и DY,
как и вектор QTY = Q_TRANSP*Y. Далее для выбранного P, вектор U определяется
в CHOL1D как решение линейного уравнения: (6(1-P)Q_TRANSP*D**2*Q + P*Q)U = QTY
Из полученного U, сглаживающий сплайн F (при данном Р) получается так, что
F(X(.)) = Y - 6(1-P)D**2*Q*U
F''(X(.)) = 6*P*U
Q_TRANSP нужна для нахождения параметров кубического сплайна и представляет
собой трехдиагональную матрицу размером (N-2) x N c линией:
1/dx[i-1] , -1/dx[i-1] - 1/dх[i] , 1/дельта_х[i]
D - это диагональная матрица из ошибок dy
'''
npm1 = npoint-1

# Считаем Q_TRANSP*D
while ee != 1:
    v[0][3] = x[1]-x[0]
    for i in range(1,npm1):
        v[i][3] = x[i+1]-x[i] # Здесь считаем dx
        v[i][0] = dy[i-1]/v[i-1][3] # Считаем элементы первого столбца матрицы Q_TRANSP*D
        v[i][1] = -dy[i]/v[i][3]-dy[i]/v[i-1][3] # Считаем элементы второго столбца матрицы Q_TRANSP*D
        v[i][2] = dy[i+1]/v[i][3] # Cчитаем элементы третьего столбца матрицы Q_TRANSP*D  
    v[npm1][0] = 0
    for i in range(1,npm1):
        v[i][4] = v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2] # Считаем элементы первого столбца матрицы Q_TRANSP*D**2*Q
    if npm1<3:
        continue
    else:
        for i in range(2, npm1):
            v[i-1][5] = v[i-1][1]*v[i][0]+v[i-1][2]*v[i][1] # Считаем элементы второго столбца матрицы Q_TRANSP*D**2*Q
    v[npm1-1][5] = 0
    if npm1<4:
        continue
    else:
        for i in range(3,npm1):
            v[i-2][6] = v[i-2][2]*v[i][0] # Считаем элементы третьего столбца матрицы Q_TRANSP*D**2*Q 
    v[npm1-2,6] = 0
    v[npm1-1,6] = 0
    
    # CONSTRUCT Q-TRANSP. * Y in QTY
    prev = (y[1]-y[0])/v[0][3]
    diff1 = np.empty(npoint)
    a = np.empty((npoint,4))
    for i in range(1, npm1):
        diff = (y[i+1]-y[i])/v[i][3]
        diff1[i] = (y[i+1]-y[i])/v[i][3] 
        a[i][3] = diff-prev # QTY(I)
        prev = diff  
    p = float(input('P: '))

# Cholid
    # COSTRUCT 6*(1-P)*Q_TRANSP*(D**2)*Q + P*R
    six1mp = 6*(1-p)
    twop = 2*p
    for i in range(1, npm1):
        v[i][0] = six1mp*v[i][4]+twop*(v[i-1][3]+v[i][3])
        v[i][1] = six1mp*v[i][5]+p*v[i][3]
        v[i][2] = six1mp*v[i][6]
    npm2 = npoint-2
    if npm2>2:
        # FACTORIZATION
        for i in range(1,npm2):
            ratio = v[i][1]/v[i][0]
            v[i+1][0] = v[i+1][0]-ratio*v[i][1]
            v[i+1][1] = v[i+1][1]-ratio*v[i][2]
            v[i][1] = ratio
            ratio = v[i][2]/v[i][0]
            v[i+2][0] = v[i+2][0]-ratio*v[i][2]
            v[i][2] = ratio
        # FORWARD SUBSTITUTION
        a[0][2] = 0 # U(1) = 0
        v[0][2] = 0
        a[1][2] = a[1][3] # U(2) = QTY(2)
        for i in range(1,npm2):
            a[i+1][2] = a[i+1][3]-v[i][1]*a[i][2]-v[i-1][2]*a[i-1][2] # U(I+1) = QTY(I+1) - V(I,2)*U(I) - V(I-1,3)*U(I-1)
        # BACK SUBSTITUTION
        a[npm1][2] = 0 # U(NPOINT) = 0
        a[npm2][2] = a[npm2][2]/v[npm2][0]
        for i in range(1,npm2):
            a[npm2-i][2] = a[npm2-i][2]/v[npm2-i][0]-a[npm2+1-i][2]*v[npm2-i][1]-a[npm2+2-i][2]*v[npm2-i][2]
    else:
        a[0][2] = 0 # U(1) = 0
        a[1][2] = a[1][3]/v[1][0] # U(2) = QTY(2)/V(2,1)
        a[2][2] = 0 # U(3) = 0
    # CONSTRUCT Q*U
    prev = 0
    for i in range(1,npoint):
        a[i][0] = (a[i][2]-a[i-1][2])/v[i-1][3] # QU(I) = (U(I)-U(I-1))/V(I-1,4)
        a[i-1][0] = a[i][0]-prev # QU(I-1) = QU(I) - PREV
        prev = a[i][0] # PREV = QU(I)
    a[npm1][0] = -a[npm1][0] # QU(NPOINT) = -QU(NPOINT)
# Cholid end

    six1mp = 6*(1-p)
    for i in range(0,npoint):
        a[i][0] = y[i]-six1mp*dy[i]*dy[i]*a[i][0]
    sixp = 6*p
    for i in range(0,npoint):
        a[i][2] = a[i][2]*sixp
    for i in range(0,npm1):
        a[i][3] = (a[i+1][2]-a[i][2])/v[i][3]
        a[i][1] = (a[i+1,0]-a[i][0])/v[i][3]-(a[i][2]/2)*v[i][3]-(a[i][3]/6)*v[i][3]*v[i][3]


#Picture 2
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    ax[0].scatter(x,y, color = 'k')
    tha1 = [a[i][1]/lat for i in range(0, npoint)]
    del tha1[-1]
    tha1.append(tha1[-1])
    tha = [a[i][0] for i in range(0, npoint)]
    ax[1].scatter(x,tha1, color = 'k')
    ax[0].plot(x, tha, color = 'r')
    ax[0].set_ylim(ymin=min(y)-abs(y[4]-y[5]), ymax = max(y)+abs(y[4]-y[5]))
    ax[1].set_ylim(ymin=min(tha1)-abs(tha1[4]-tha1[3]), ymax = max(tha1)+abs(tha1[4]-tha1[3]))
    ax[0].errorbar(x,y,yerr = dy, fmt = 'none', elinewidth=0.5, ecolor = 'k')
    ax[0].grid(True, 'major', linewidth=0.3)
    ax[1].grid(True, 'major', linewidth=0.3)
    print('Cмотрите график')
    plt.show()
    ee = float(input('Введите (1 - Хороший фит, 0 - Плохой фит (изменить p): '))
    if ee == 1:
        file = open(f'D:\\{filename}_AFTER_FIT.txt','w')
        #file2 = open('D:\\CTE_BEFORE_FIT.txt','w')
        number_points = int(input('Сколько точек было добавлено: '))
        file.write('0' + '  ' + '0' + '\n')
        print(f'Данные записаны в файл {filename}_AFTER_FIT.txt')
        #file2.write('0' + '  ' + '0' + '\n')
        for i in range(number_points,len(x)):
            file.write(str(x[i]) + '  ' + str(tha1[i]) + '\n')
            #file2.write(str(x[i]) + '  ' + str(y[i]) + '  ' + str(tha[i]) + '\n')
        file.close()
        #file2.close()
