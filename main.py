import numpy as np
import random as rn
if __name__ == "__main__":
    def dist(a,b,h_cntr,):# функция расчёта расстояния от неравенства до центра гиперкуба
        distantion=(abs(np.dot(a,h_cntr)-b))/np.linalg.norm(a,ord=2) 
        return distantion
    def AddSupport(n,alpha): #добавление исходной матрицы
        A=np.eye(n)
        A1=np.eye(n)
        for i in range(0,n): #генерация единчного массива
            A1[i,i]*=-1
        A=np.vstack((A, A1)) 
        A1=np.ones(n)
        A=np.vstack((A, A1))
        B=[] # генерация массива свободных членов
        for i in range(0,n):
            B=np.append(B,alpha)
        for i in range(n,2*n):
            B=np.append(B,0)
        B=np.append(B,alpha*(n-1)+100)
        return A,B
    def projection(h_cntr,a,b): #вычисление проекции
        proj=h_cntr-((np.dot(a,h_cntr)-b)/(np.linalg.norm(a,ord=2)**2))*a
        return proj
    def like(A,B,a,b): #проверка на схожесть
        Backup_A=np.zeros((len(A),n))
        Backup_a=[]
        Backup_B=[]
        Backup_b=0
        for i in range (0,n):
            Backup_a.append(a[i]/np.linalg.norm(a,ord=2))
        Backup_b=b/np.linalg.norm(a,ord=2)
        for i in range (0,len(A)):
            Backup_B.append(B[i]/np.linalg.norm(A[i][:],ord=2))
            for j in range (0,n):
                Backup_A[i][j]=A[i][j]/np.linalg.norm(A[i][:],ord=2)
        if (np.linalg.norm(Backup_A-Backup_a,ord=2) < L_max):
            for i in range(len(A)):
                if abs(Backup_B-Backup_b)< S_min:
                    return True
                else:
                    return False
        else:
            return False
    def point_func(func,c): #расчёт целевой функции для проекции
        state=c*func
        state = np.sum(state)
        return state
    #основаная часть
    n=int(input('Введите n ')) #размерность
    d=int(input('Введите d ')) #количество случайных неравенств
    alpha= int(input('Введите alpha ')) #положительный масштабирующий  коэффицентв
    thetha=int(input('Введите thaeta ')) #полож постоянный множитель должен быть меньше alpha/2
    ro=int(input('Введите ro '))
    S_min= float(input('Введите S_min ')) #парамеры программы для генерации
    L_max=float(input('Введите L_max ') )# L_max<0,7
    a_max=int(input('Введите a_max ')) # полож или отриц будет коэффицент ЛФ
    b_max=int(input('Введите b_max '))# коэффицент ЛФ """
    k=0 #счётчик для сгенеированных переменных
    c=[] # массив для целевой функции
    a=np.zeros(n) #массив для новых переменных функций
    rsign=[-1,1] #рандомизатор для новых переменных
    h_cntr=np.zeros(n) #массив для ценра гиперкуба
    for j in range (n,0, -1): # генерация ЦФ
        c = np.append(c,thetha*j)
    for j in range(0,n): # генерация центровой координаты
        h_cntr[j]=alpha/2
    A,B=AddSupport(n,alpha)

    while k<d : #начало цикла генерации
        for j in range(0,n):
            a[j]=rn.choice(rsign)*rn.randint(0,a_max) # генерация случайного коэффициента
        b=rn.choice(rsign)*rn.randint(0,b_max) # генерация случайного свободного члена
        if np.dot(h_cntr,a) <= b: # валидация сгенерированной переменной
            if ((dist(a,b,h_cntr) > ro) or (dist(a,b,h_cntr)<=thetha)) and (point_func(projection(h_cntr,a,b),c) < point_func(h_cntr,c)) and (like(A,B,a,b)==False) : # т.к Питон жёстко структурированный ЯП нужно проверять всё сразу
                F= np.array([a]) # добавление коэффициентов рпошедших валидацию
                A= np.append(A,F,axis=0)
                B=np.append(B,b)# добавление переменных прошедших валидацию
                k+=1
            else:
                for j in range(0,n):
                    a[j]=rn.choice(rsign)*rn.randint(0,a_max) # генерация случайного коэффициента
                b=rn.choice(rsign)*rn.randint(0,b_max) # генерация случайного свободного
        else:
            for j in range(0,n): # если расстояние от центра до свободного члена больше b то заменяем все переменные a на отрицалтельные 
                a[j]=-a[j]
            b=-b
            if ((dist(a,b,h_cntr) > ro) or (dist(a,b,h_cntr)<=thetha)) and (point_func(projection(h_cntr,a,b),c) < point_func(h_cntr,c)) and (like(A,B,a,b)==False) :  # Такая же проверка что и выше
                F= np.array([a])# добавление коэффициентов рпошедших валидацию
                A= np.append(A,F,axis=0)
                B=np.append(B,b)# добавление переменных прошедших валидацию
                k+=1
            else:
                for j in range(0,n):
                    a[j]=rn.choice(rsign)*rn.randint(0,a_max) # генерация случайного коэффициента
                b=rn.choice(rsign)*rn.randint(0,b_max) # генерация случайного свободного
    print("Вывод значений сгененированных задач ЛП \n")
    print(A) #Вывод полученных значений в частности матрица коэффициентов
    print("\n")
    print("Матрица свободных членов сгененрированых задач ЛП\n")
    B=np.reshape(B,(-1,1))
    print(B) # матрица свободных членов сгенерированных задач
    print("\n")
    summurize = np.append(A,B,axis=1)
    print("\n")
    print("Общая матрица коээфицентов и свободных членов\n")
    print(summurize) 
    c =np.append(c,0)
    print("\n")
    print("Матрица свободных членов ЦФ\n")
    print(c) # матрица свободных членов ЦФ
