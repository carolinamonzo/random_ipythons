#! /usr/bin/env python3.6
# -*- coding:utf-8 -*-

# Cuestionario en grupos: Numeros perfectos, version inicial no mejorada ni optimizada
# Alumnos:
# - Sancho Blanco Fernandez
# - Celia Borja Almarcha
# - Alejandro Flores Leon
# - Nuria Gomez Cebrian
# - Carolina Monzo Catalunya
# - Lorena E. Rosaleny Peralvo

from time import time

def divisores(num):
    '''
    Funcion para calcular los divisores de un numero
    Input: numero de interes
    Output: lista de divisores del numero de interes
    '''
    lista_div = []

    # Para ser divisor, el modulo de la division sera igual a 0
    for i in range(1, num):

        if num % i == 0:
            lista_div.append(i)

    return(lista_div)

def es_perfecto(num, lista_div):
    '''
    Funcion para comprobar si un numero es perfecto
    Input: lista de divisores
    Output: numeros perfectos
    '''

    # Un numero perfecto es aquel para el que la suma de sus divisores es igual a el mismo
    perfecto = False
    
    if sum(lista_div) == num:
        perfecto = True

    return(perfecto)


def main():
    '''
    Funcion para correr el programa, calcula los 10 primeros numeros perfectos
    Input:
    Output: 10 numeros perfectos
    '''
    i = 0
    lista_perf = []
    num = 1

    print('Numeros perfectos:')

    # Vamos a calcular tambien el tiempo de ejecucion
    tiempo_inicial = time()

    # Hasta que no tengamos los primeros 10 numeros perfectos, calculamos
    while i < 10:

        lista_div = divisores(num)

        perfecto = es_perfecto(num, lista_div)

        # Cuando conseguimos un perfecto, lo guardamos y pasamos al siguiente
        if perfecto:
            i += 1
            print(num)

        num += 1

    # Tiempo de ejecucion
    tiempo_final = time()

    tiempo_ejecucion = tiempo_final - tiempo_inicial
    print('El tiempo de ejecucion es: {}'.format(tiempo_ejecucion))

if __name__ == '__main__':

    main()
