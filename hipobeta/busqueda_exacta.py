# !/usr/env/python                                                     #
# -*- coding: utf-8 -*-                                                #
# Título: Busqueda de secuencias empleando el algoritmo de karp-rabin  #
# Autores:                                                             #
#   - Azahara Maria Fuentes Trillo                                     #
#   - Adrián Santiago Ortiz                                            #
#   - Daniel Carrillo Bautista                                         #
#   - Iris Manzano Blasco                                              #
# Fecha: 30/12/2016                                                    #
########################################################################

import math
from multiprocessing import Process, Queue
import string
import time

"""
Función que crea un diccionario con todos los caracteres del abecedario.
"""
def crear_diccionario():

    diccionario = dict.fromkeys(string.ascii_lowercase, 0)
    for i, j in enumerate(diccionario):
        diccionario[j] = i

    return diccionario


"""
Función que calcula el valor asignado al patrón buscado, empleando
los valores asignados a cada una de las letras del diccionario.
"""
def calcular_codigo_cadena(diccionario, cad, p):

    cod = 0
    long = len(cad)
    for i in range(0, long):
        cod = cod + (diccionario[cad[i].lower()] * p ** i)

    return cod


"""
Función que calcula el valor asociado a cada una de las secuencias de
tamaño n, contenidas en la secuencia sujeta a la búsqueda.
"""
def calcular_codigo_encadenado(diccionario, codigo, anterior, cadena, primo):

    cod = ((codigo - diccionario[anterior.lower()])/primo) + \
            (diccionario[cadena[-1].lower()] * primo ** (primo - 1))

    return cod

"""
Función que compara el patrón incluido por el usuario con cada
una de las posiciones encontradas por el algoritmo. Las colisiones
típicas de este procedimiento son así descartadas.
"""
def verificar_coincidencia(i, cod, codigo_sec, secuencia, sec, primo):

    if cod == codigo_sec:
        if secuencia[i : i + primo] == sec:

            return True
        else:

            return False


"""
Función que realiza la búsqueda del patrón en la secuencia mediante el
algoritmo de Karp-Rabin, empleando el número de procesos seleccionado
por el usuario.
"""
class busqueda_procesos(Process):

    def __init__(self, diccionario, secuencia, i1, i2, codigo_sec, \
                 primo, sec, cola):

        Process.__init__(self)
        self.diccionario = diccionario
        self.secuencia = secuencia
        self.i1 = i1
        self.i2 = i2
        self.codigo_sec = codigo_sec
        self.primo = primo
        self.sec = sec
        self.cola = cola


    def run(self):

        l = []
        cod = 0
        pos = 0

        if self.i1 == 0:
            self.i1 = 0
        else:
            self.i1 = self.i1 - len(self.sec)

        intervalo_sup = self.i2 - len(self.sec)
        for i in range(self.i1, intervalo_sup):

            if i == self.i1:

                cod = calcular_codigo_cadena(self.diccionario, \
                                             self.secuencia[i : i + \
                                             self.primo], self.primo)

                if verificar_coincidencia(i, cod, self.codigo_sec, \
                                          self.secuencia, self.sec, \
                                          self.primo):
                    pos = i
                    l.append(pos)

            else:
                cod = calcular_codigo_encadenado(self.diccionario, cod, \
                                                 self.secuencia[i - 1], \
                                                 self.secuencia[i : i + \
                                                 self.primo], self.primo)
                if verificar_coincidencia(i, cod, self.codigo_sec, \
                                          self.secuencia, self.sec, \
                                          self.primo):
                        pos = i
                        l.append(pos)
        self.cola.put(l)


"""
Función que procesa un archivo en formato fasta y devuelve la secuencia
contenida en el mismo para la posterior búsqueda.
"""
def procesar_archivo(fhand):

    cabecera = []
    texto = []
    t = fhand.readlines()
    for j, n in enumerate(t):
        t[j] = n.strip()
        if t[j].startswith(">"):
            cabecera.append(t[j])
            texto.append("I")
        else:
            texto.append(t[j])

    texto = "".join(texto)

    return texto, cabecera


"""
Función que realiza la búsqueda del patrón en la secuencia empleando el
algoritmo de Karp-Rabin de forma secuencial.
"""
def secuencial(diccionario, cadena, codigo_sec, primo, sec):
    l = []
    cod = 0
    pos = 0

    intervalo_sup = len(cadena) - len(sec) + 1
    for i in range(0, intervalo_sup):

        if i == 0:
            cod = calcular_codigo_cadena(diccionario, cadena[i : i + primo], \
                                         primo)
            if cod == codigo_sec:
                if cadena[i : i + primo] == sec:
                    pos = i
                    l.append(pos)
        else:
            cod = calcular_codigo_encadenado(diccionario, cod, cadena[i - 1], \
                                             cadena[i : i + primo], primo)
            if cod == codigo_sec:
                if cadena[i : i + primo] == sec:
                    pos = i
                    l.append(pos)

    return l

"""
Función que imprime y da formato a la solución obtenida por la búsqueda.
"""

def printar_solucion(posicion, cabecera, secuencia):

    index_cabecera = secuencia[:posicion].count("I") - 1
    posicion_local = posicion - secuencia.rfind("I", 0, posicion) - 1

    return "\t\t" + cabecera[index_cabecera] + ": " + str(posicion_local)

"""
Función principal.
"""
def main():


    nombre = input("Introduzca el nombre del archivo con la secuencia: ")

    try:
        fhand = open(nombre, "r")
    except:
        print ("Error. El fichero indicado no existe o no está" + \
               " en el directorio actual.")
    else:
        sec = input("Introduzca el patron a buscar: ")
        n_procesos = int(input("Introduzca el numero de procesos: "))
        secuencia, cabecera = procesar_archivo(fhand)

        primo = len(sec) # longitud del patrón

        # Se calcula el valor asociado a la cadena a buscar
        diccionario = crear_diccionario()
        codigo_sec = calcular_codigo_cadena(diccionario, sec, primo)

        ## 1. Proceso secuencial

        print ("1. Versión secuencial:")

        t1_sec = time.time()
        lista_posiciones = (secuencial(diccionario, secuencia, codigo_sec, \
                            primo, sec))
        t2_sec = time.time()

        if not lista_posiciones:
            solution = False
            print ("\tPatron no encontrado.")
        else:
            print ("\tSe ha encontrado en las posiciones: ")
            solution = True
            for posicion in lista_posiciones:
                print (printar_solucion(posicion, cabecera, secuencia))

        print ("-" * 80)

        ## 2. Proceso concurrente

        print ("2. Versión concurrente (" + str(n_procesos) + " procesos)")
        procesos = ["t"] * n_procesos
        q = Queue()
        for i in range(n_procesos):

            interval1 = i * int(math.ceil(len(secuencia) / float(n_procesos)))
            interval2 = interval1 + int(math.ceil(len(secuencia) / \
                                                      float(n_procesos)))
            procesos[i] = busqueda_procesos(diccionario, secuencia, \
                                            interval1, interval2, \
                                            codigo_sec, primo, sec, q)

        t1_procesos = time.time()
        for i in range(n_procesos):
            procesos[i].start()

        for j in range(n_procesos):
            procesos[j].join()

        t2_procesos = time.time()


        if solution == False:

            print ("\tPatron no encontrado.")
        else:
            print ("\tSe ha encontrado en las siguientes posiciones: ")
            l = []
            while not q.empty():
                l += q.get()

            l = sorted(l)
            for posicion in l:
                print (printar_solucion(posicion, cabecera, secuencia))


        print ("-" * 80)

        ## 3. Tiempos resultantes

        t_sec = t2_sec - t1_sec
        t_conc =  t2_procesos - t1_procesos

        print ("Tiempos resultantes")
        print ("\t1. Versión secuencial: " + str(round(t_sec, 5)) + " s.")
        print ("\t2. Versión concurrente: " + str(round(t_conc, 5)) + " s.")
        print ("\t3. Aceleración:" + str(round(t_sec / t_conc, 2)))


if __name__ == "__main__":
    main()
