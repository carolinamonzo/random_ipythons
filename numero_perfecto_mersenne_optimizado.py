#! /usr/bin/env
# -*- coding:utf-8 -*-

# Cuestionario en grupos: Numeros perfectos, version mejorada con mersenne
# Alumnos:
# - Sancho Blanco Fernandez
# - Celia Borja Almarcha
# - Alejandro Flores Leon
# - Nuria Gomez Cebrian
# - Carolina Monzo Catalunya
# - Lorena E. Rosaleny Peralvo

from time import time
from math import sqrt
from multiprocessing import Process, Queue
from threading import Thread

def es_primo(num):
    '''
    Funcion para comprobar si un numero es primo
    Input: numero de mersenne de interes
    Output: numero primo
    '''
    
    ## CAMBIO RESPECTO A LAS VERSIONES ANTERIORES

    primo = True
    max_num = sqrt(num) + 1
    mit = 2

    # Utilizamos el numero input como dividendo y su raiz cuadrada mas uno
    # como divisor, si el resto es 0, el numero no sera primo

    while mit < max_num:
        
        if num % mit == 0:
            primo = False
            # Si el numero tiene divisor, no es primo
            break
        mit += 1

    return(primo)


def es_primo_mersenne(num):
    '''
    Funcion para comprobar que los numeros primos lo son tambien de mersenne
    Input: numero de interes
    Output: comprobacion y valor
    '''

    # Comprobacion de si el numero de mersenne tambien es primo
    mersenne = False
    num_mersenne = 2**(num) - 1

    if es_primo(num_mersenne):
        mersenne = True

    return(mersenne, num_mersenne)


def es_perfecto(num, final, queue):
    '''
    Funcion para comprobar si un numero primo y de mersenne es perfecto
    Input: contador, inicio y final del intervalo
    Ouptut: numeros perfectos
    '''
    lista_perfectos = []

    # Comprobar un numero es primo
    while num < final:
        primo = es_primo(num)

        # Comprobar si el numero es de mersenne
        if primo:
            mersenne, num_mersenne = es_primo_mersenne(num)

            # Comprobar si el numero de mersenne es perfecto
            if mersenne:
                perfecto = num_mersenne * 2 ** (num - 1)
                lista_perfectos.append(perfecto)

        num += 1
    print(lista_perfectos)
    queue.put(lista_perfectos)


def main():
    '''
    Programa para calcular el numero perfecto de mersenne
    Output: 10 numeros perfectos
    '''

    procesos = int(input('Introduce el numero de procesos:'))
    if procesos > 10:
        print('El numero de procesos es mayor al numero de numeros a calcular')
        exit(1)

    time_start = time()

    inicio = 0
    final_intervalo = 0

    decimo_perfecto = 191561942608236107294793378084303638130997321548169216

    aumento = decimo_perfecto / procesos

    p = []

    # Iniciar la cola
    queue = Queue()
    i = 0

    # Hacer un bucle para crear tantos procesos como se necesiten
    for i in range(procesos):
        
    # Limites del total de valores a calcular por cada proceso
        final_intervalo = inicio + aumento + final_intervalo        
        # Lanzar la funcion en cada proceso
        p.append(Process(target=es_perfecto, args=(inicio, final_intervalo, queue)))
        print(p)

    for i in p:
        print(i)        
        i.start()
        
        # Sacar el valor obtenido por cada proceso por pantalla
    #for i in range(procesos):
      #  p[i].join()

        resultado = queue.get()
        print(resultado)

    # Calcular el tiempo de ejecucion
    time_end = time()
    tiempo_ejecucion = time_end - time_start
    print('El tiempo de ejecucion es: {}'.format(tiempo_ejecucion))



if __name__ == '__main__':

    main()
