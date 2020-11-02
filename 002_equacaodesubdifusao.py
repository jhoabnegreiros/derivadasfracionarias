import numpy as np
from sympy import gamma, sin
import json

def gravar(log):
    with open('log.json', 'w+') as arq:
        json.dump(log, arq, indent=True)

def condicaoinicial(xnos):
    cinicial = np.linspace(0, np.pi, xnos + 1)
    cinicial = np.array(list(map(lambda x: 0 * sin(x*np.pi), cinicial)))
    return cinicial

def omega(tnos, oderi):
    omega = (-1)**(tnos)*gamma(2-oderi)/(gamma(tnos+1)*gamma(2-oderi-tnos))
    return omega

def termofonte(xi, tj, oderi):
    return (2.5*tj**1.5+(gamma(3.5)/gamma(3.5-oderi)*tj**(2.5-oderi)))*sin(xi)


#  Entrada dos Dados

exnos = int(input('Digite o número de intervalos em x = '))
etnos = int(input('Digite o número de intervalos em t = '))
eoderi = float(input('Digite a ordem da derivada = '))

#  Programa Principal
Dt = 1/etnos
mu = 1*Dt**eoderi
U = np.zeros((etnos+1, exnos+1))

n = 0
vetor = np.zeros(etnos+1)
log = {'vetor': {}, 'U': {}}
while n <= etnos-1:
    vetor[n] = omega(n, eoderi)
    #log['vetor'][f'{n}'] = str(vetor[n])  # Add o vetor no dicionário para salvar no arquivo JSON
    for i in range(1, exnos):
        aux = 0
        for k in range(0, n+1):
            aux = aux + vetor[n-k]*U[k][i]
        U[n+1][i] = U[n][i] + mu * aux + Dt * termofonte(np.pi*i/exnos, 1*(n+1)/etnos, eoderi)
        #log['U'][f'{n+1}-{i}'] = U[n+1][i]  # Add U no dicionário para salvar no arquivo JSON
    n = n + 1
    #gravar(log)  # Guarda o vetor com os valores omega e U no arquivo 'log.json'


print(U[etnos])
