#!/usr/bin/env python
# coding: utf-8

# # Primeira Lei da termodinâmica

# ## Resolvendo problemas via integração numérica
# 
# <font size = "4">Como percebido nos exemplos anteriores, dependendo da equação de estado utilizada, pode ser necessário a resolução de integrais complicadas. Por este motivo, muitos problemas em termodinâmica são resolvidos por integração numérica dentro do intervalo de valores desejados. A linguagem Python contem uma série de bibliotecas e métodos computacionais que permitem a resolução de uma integral definida numericamente sem previamente resolvê-la manualmente. Este método será ilustrado repetindo o cálculo do trabalho realizado nos Exemplos 3-4, mas desta vez utilizando um método de integração numérica.
# 
# <font>

# Os seguintes módulos precisam ser carregados para resolução dos exercícios:

# In[1]:


#bibliotecas necessárias para resolução de exercícios
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt


# ## Exemplo 1: Cálculo do trabalho em uma expansão reversível de um gás considerando diferentes equações de estado
# 
# <font size="4"> Considere a compressão reverssível de 1.0 mol de um gás ideal, de 22.4 L para 1.0L mantendo T constante a 0ºC. Calcule o trabalho realizado neste processo considerando que o gás obedece a: a) equação de van der Waals, b) a equação do virial, e c) a equação de Redlich-Kowng.
# 
# ### Solução 
# 
# O trabalho reversível é o trabalho teórico máximo que o sistema pode realizar. Este trabalho é teórico pois considera que o trabalho é feito em n etapas infinitesimais. O valor do trabalho dependerá do sistema e da equação de estado que descreve este sistema. Por este motivo, em primeiro lugar é preciso encontrar um expressão para o trabalho p-V de um gás ideal. No trabalho reversível podemos considerar que a pressão de oposição se iguala a pressão do sistema, sendo esta calculada pela equação de estado. 
#     
# A biblioteca SciPy possibilita o uso da função quad, baseada na biblioteca Fortran QUADPACK, para calcular integrais definidas no intervalo de "a" até "b". Para utilizar esta função é preciso definir a função matemática, segundo a sintaxe do Python, e posteriormente proceder com integração numérica. 
# 
# $$ dW = -pdV $$
#     
# Sintaxe a ser usada: resultado, erro = quad(função,valor inicial, valor final)
# 
#        

# In[2]:


#a) Trabalho considerando a eq. de van der Waals

# Primeiro definimos a função segundo a sintaxe do Python
#Note que será resolvida a integral de -pdV.

R = 8.314          # m^3 Pa / K mol   
a = 0.1380         # m^3 Pa / mol^2
b = 0.0319e-3      # m^3 / mol  

V1 = 22.4e-3       # m^3
V2 = 1e-3          # m^3
T = 273            # K

def p(V):
    
    return -(R*T/(V-b) - a/(V**2))

# Integração numérica da função p(V)

Wquad, err = quad(p, V1, V2)  # Wqad é o objeto que armazenará o resultado, e o objeto err armazenará o erro da integração numérica


print('Trabalho (J/ mol)= ' + format(Wquad , '6.3f'))
print('erro da integração = ' + format(err, '6.3e'))


# <font size = "4">Neste caso ocorre uma pequena diferença entre os valores calculados usando o resultado analítico e o método de integração numérica, mas o erro é aceitável. <font>

# In[3]:


#b) Trabalho considerando a eq. do virial

# Primeiro definimos a função segundo a sintaxe do Python

B = -22.0e-6        # m**3 / mol 
C = 1100e-12        # m**6 / mol**2


def p(V):
    
    return -R*T/V *(1 + B/V + C/V**2) 

# Agora é calculado a integral (Wquad), o erro associado a integração numérica (err), da função definida anteriormente
# nos intervalos de V1 a V2

Wquad, err = quad(p, V1, V2)



# Verifique que o resultado é basicamente o mesmo do obtido anteriormente.

print('Trabalho (J/ mol)= ' + format(Wquad , '6.3f'))
print('erro da integração = ' + format(err, '6.3e'))


# In[4]:


#c) Trabalho considerando a eq. de Redlich-Kwong

# Primeiro definimos a função segundo a sintaxe do Python

a = 17.16*101.3e-3        # m**6 K^0,5 Pa / mol^2 
b = 0.0221e-3             # m**3 / mol

def p(V):
    
    return -(R*T/(V-b) - a/(V*(V+b)*T**0.5))

Wquad, err = quad(p, V1, V2)




print('Trabalho (J/ mol)= ' + format(Wquad , '6.3f'))
print('erro da integração = ' + format(err, '6.3e'))


# <font size = "4">Os métodos de integração numérica simplificam o  resolução de integrais complicadas e podem levar a resultados comparativamente iguais ao do método analítico. Como pode ser observado, não foi preciso integrar analiticamente a função desejada, definir a equação resultante no script e, finalmente, fazer a substituição dos valores nas variáveis da equação. O que levou a uma grande simplificação na resolução do exercício. No entanto, métodos de integração numérica introduzem um erro no resultado da integração e devem ser escolhidos de acordo com o problema matemático a ser resolvido. <font>

# <font size ="4"> Todos os exemplos anteriores que envolvem o polinômio da capacidade calorífica e obtenção do polinômio da capacidade calorífica poderiam ser resolvidos por integração numérica. No entanto, é importante conhecer as duas formas de resolução dos problemas para que seja escolhido o método mais adequado de acordo com o problema abordado.<font>

# ## Sugestão de exercício:
# 
# <font size ="4">Resolva os exemplos 8 e 9 do capítulo anterior usando o método de integração numérica e compare os resultados obtidos com a solução analítica do problema. <font>
