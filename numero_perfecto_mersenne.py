#! /usr/bin/env python3.6
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

def es_primo(num):
    '''
    Funcion para comprobar si un numero es primo
    Input: numero de mersenne de interes
    Output: numero primo
    '''

    primo = True

    # Calculamos los divisores hasta la mitad del numero de interes
    # para calcular divisores comprobamos que el modulo de la division es 0

    for mit in range(2, int(num/2)+1):
        
        if num % mit == 0:
            primo = False
            # Si el numero tiene divisor, no es primo
            break

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


def main():
    '''
    Programa para calcular el numero perfecto de mersenne
    Output: 10 numeros perfectos
    '''
    
    i = 0
    num = 2
    max_perfectos = 10

    time_start = time()
    print('Los numeros perfectos calculados a partir de numeros de mersenne son:')
    # Comprobar si un numero es primo

    while i < max_perfectos:
        primo = es_primo(num)
        
        # Compromar si el numero es de mersenne
        if primo:
            mersenne, num_mersenne = es_primo_mersenne(num)

            # Comprobar si el numero de mersenne es perfecto
            if mersenne:
                print(num_mersenne * 2 ** (num - 1))
                i += 1

        num += 1

    # Calcular el tiempo de ejecucion
    time_end = time()
    tiempo_ejecucion = time_end - time_start
    print('El tiempo de ejecucion es: {}'.format(tiempo_ejecucion))



if __name__ == '__main__':

    main()
