# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 10:23:33 2019

@author: vitor
"""

import numpy as np

################### MÉTODO JACOBI #####################
def met_jacobi(A, b, x0, tol, max_iter):
    n = np.alen(A) # n é o número de linhas
    x = np.zeros(n, dtype='double') # Solução atual
    x[0] = b[0]
    xant = x0 # Solução da iteração anterior
    for k in range(max_inter): # Iterações do método
        for i in range(n): # Iterações para cada incógnita
            x[i] = b[i] # Termo constante
            for j in range(i): # Incógnitas anteriores
                x[i] -= A[i,j]*xant[j]
            for j in range(i+1, n): # Incógnitas posteriores
                x[i] -= A[i,j]*xant[j]
            x[i] /= A[i,i] # Divisão pelo coeficiente da incógnita atual
        erro = np.linalg.norm(x-xant, np.inf) # Cálculo do erro
        print("Iteração {k:3d}: ".format(k=k+1) +
              "x={x}, ".format(x=np.round(x,8)) +
              "Erro={e:+5.5f}".format(e=erro))
        if (erro < tol): 
            return x
        xant = np.copy(x)
        
# Pegar o maior valor da linha e da coluna para o pivo
def linha_pivo (M,i,j):
    n = np.alen(M)
    for k in range (n):
            if (M[i][j] <= M[k][j]): # verifica as celulas 
              T= k
    return T    

# troca o pivo com a linha seguinte
def trocar_linha (M,pivo,lin):
    
    if pivo != lin:
        aux = np.copy(M[pivo,:])
        M[pivo,:] = np.copy(M[lin,:])
        M[lin,:] = aux
       
# Função para escalonar a matriz  
def escalonamento (M):
    n = np.alen(M)
    for k in range (n):
        l = linha_pivo(M,k,k)
        trocar_linha(M,l,k)
        pivo = M[k,k]
        if pivo == 0:
            return False
        for l in range (k+1,n):
            M[l,:] = M[l,:] - M[k,:] * M[l,k]/pivo
    return True

# Calculo do determinante
def determinante (M):
    if escalonamento(M) == True:
        n = np.alen(M)
        determinante = 1
        for i in range (n):
            determinante *= M[i,i] # calculo da diagonal principal
        return determinante
#Função necessária se o sistema convergirá
def raio_espectral(M):
    av, _ = np.linalg.eig(M)
    return max(abs(av))  
  
# tenste de convergencia
def converge(M):
    D = np.diag(np.diag(M)) #diagonal principal
    L = np.tril(M)-D #diagonal inferior
    U = np.triu(M)-D #diagonal superior
    T = -np.linalg.inv(D).dot(L+U)
    if raio_espectral(T) <= 1:
        return print(" \n >>>>>>>>>> Sistema possui solução <<<<<<<<<< ")
    else:
        return print(" \n >>>>>>>>>> Sistema não possui solução <<<<<<<<<< ")
#Informações colhidas do usuário
n = int(input('Informe o número de variáveis do sistema: \n ' )) 
max_inter= int(input('Informe o número máximo de interações: \n '))
tol = int(input('Informe a tolerância: \n '))
constantes = input ('O termos das constantes: \n ')
const = constantes.split()
const = [float(x) for x in const]
M = []
#Séra onde o usuário informará cada valor de cada posição da matriz
print('\n\nInforme a matriz estendida do sistema (linha/Coluna): \n ' )
for i in range (n):
    linha = input ("Informe a linha "+str(i+1)+ ": ")
    lin_num = linha.split()
    lin_num = [float(x) for x in lin_num]
    M.append(lin_num)
M = np.array(M, dtype = "double")
const = np.array(const, dtype = "double")
print(' \n ********** Matriz Estendida: ********** \n',M)

#Opção do usuário em continuar ou sair
print (""" \n Deseja Prossequir?: 
        [ 1 ] Sim
        [ 2 ] Não """)
opcao = int(input("Sua opção: "))

if opcao == 1: 
     Mat = np.copy(M)
     x0 = np.array([1, 1, -1], dtype='double') 
    #Chamada da Função Escalonamento
     escalonamento(M)
     print("\n\n ********** Matriz escalonada ********** ")
     print(M)
    #Armazenamento da determinante do sistema em uma variável
     d=determinante (M)
     print ("\n\n ********** Determinante ********** ")
     print(d)
     converge(M)
     print("\n\n ********** Metodo de Jacobi ********** \n ")
     resultado = met_jacobi (Mat, const , x0, tol, max_inter)
     print("\n\n ********** SOLUÇÃO DO SISTEMA ********** ")
     print(resultado)    
elif opcao == 2:
    print('Programa Encerrado!')