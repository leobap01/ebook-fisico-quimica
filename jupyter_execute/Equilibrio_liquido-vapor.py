#!/usr/bin/env python
# coding: utf-8

# ## Exemplo de como resolver sistema de equações não lineares em python
# 
# Este exemplo servirá como base para resolução de exercícios de equilíbrio químico. Uma das formas de se calcular a composição de equilibrio de reações químicas é pela solução de um sistema de equações não lineares. O exemplo abaixo mostrará duas maneiras de se resolver sistemas de equações não lineares em python. O primeiro faz uso do módulo fsolve e o segundo do módulo least_squares.

# In[1]:


# Módulos necessários para resolver sistema de equações
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import least_squares
import matplotlib.pyplot as plt


# ## Ex.) Resolva o sistema de equações abaixo:
# 
# $$ 4x^2 + y^2  = 1 $$
# $$ x^2 + 9y^2  = 1 $$
# 

# In[2]:


# 1) Definir o sistema de equações. Isto é feito escrevendo as funções no formato do python. Neste exemplo
# será resolvido o sistema:
# 4x**2 + y**2 = 1
# x**2 + 9y**2 = 1
# Cada função é escrita de forma que sejam calculadas suas raízes

def equations(vars):
    x, y = vars
    eq1 = 4*x**2 + y**2 - 1
    eq2 = x**2 + 9*y**2 - 1
    return [eq1, eq2]

#Resolvendo com fsolv. Os números entre parênteses são o chute inicial para x e y.

x, y =  fsolve(equations, (3, 2))
print('Solução com fsolve= ' +  format(x , ' 6.5f') + format(y , ' 6.5f'))

print()

#Resolvendo com lest_squares. Os números entre parênteses são o chute inicial para x e y.
res = least_squares(equations, (2, 0.5))

print('Solução com least_squares= ' , res.x)

  


# ## Exemplo 8.1 (Levine)
# Use a equação de eestado RK para estimar a pressão de vapor e os volumes molares de saturação das fases líquidas e gasosas do $C_3 H_8$ a 25ºC. 

# Para estudar o equilíbrio líquido vapor, as três equações abaixo precisam ser satisfeitas:
# 
# $$ p = \frac{1}{V_m^v - V_m^l} \left( RT ln \frac{V_m^v - b}{V_m^l -b} -\frac{a}{bT^{1/2}} ln \frac{V_m^v(V_m^l + b)}{V_m^l(V_m^v + b)} \right) $$
# 
# $$ p = \frac{RT}{V_m^v - b}  - \frac{a}{V_m^v(V_m^v + b)T^{1/2}} $$
# 
# $$ p = \frac{RT}{V_m^l - b}  - \frac{a}{V_m^l(V_m^l + b)T^{1/2}} $$
# 
# Também usaremos os valores tabelados de a e b para o propano.

# In[3]:


#Variáveis

R = 82.06 # atm cm³/K mol
T = 298.15 # K
a = 1.80e8 # cm^6 atm K^1/2 mol^-2
b = 62.7   # cm^3/mol



# definindo o sistema de equações



def equations(vars):
    pvap, Vl, Vv = vars
    
    eq1 = 1/(Vv - Vl)*(R*T*np.log((Vv - b)/(Vl - b)) - a/(b*T**0.5)*np.log((Vv*(Vl+b)/(Vl*(Vv +b))))) - pvap
    
    eq2 = R*T/(Vv - b) - a/(Vv*(Vv + b)*T**0.5) - pvap 
    
    eq3 = R*T/(Vl - b) - a/(Vl*(Vl + b)*T**0.5) - pvap
    
    return [eq1, eq2, eq3]

#Resolveremos com lest_squares, pois desta maneira podermos incluir algumas restrições.
#Usaremos as seguintes restrições:
# 1 atm < pvap < pc
# b < Vl < Vc
# Vc < Vv


res = least_squares(equations, (30, 80, 1200), bounds=((1, 70, 200 ),(40, 200, 3000 )) )

print('Solução com least_squares= ' , +  res.x)


# ## **Exercício em substituição ao sexto teste. Data da entrega: 07 de abril, até 20 horas.**
# 
# 
# 
# 
# 1.   Repita o exercício anterior para o gás de van der Waals
# 2.   Resolva o exercício 8.15b. É o mesmo processo anterior, mas usando a equação de estado de Soave-Redlich-Kwong
# 
# 

# In[4]:


#Questão 2
#variáveis
R=82.06 # atm cm³/K mol
T=298.15 # K
Tc=369.8 # K
Pc=41.9 # atm


w=0.153
m=0.480+1.574*w-0.176*w**2
a=0.42748*((R**2*Tc**2)/(Pc))*(1+m*(1-(T/Tc)**0.5))**2

print('O valor de a(T) é', a)
print('O resultado bateu com o livro')


b = 0.08664*R*Tc/Pc

print('b=', b)


# In[5]:


#Variáveis


#a= 10082000

# definindo o sistema de equações


def equations(vars):
    pvap, Vl, Vv = vars
    
    eq1 = 1/(Vv - Vl)*(R*T*np.log((Vv - b)/(Vl - b)) - a/b * np.log((Vv*(Vl+b)/(Vl*(Vv +b))))) - pvap
    
    eq2 = R*T/(Vv - b) - a/(Vv*(Vv + b)) - pvap 
    
    eq3 = R*T/(Vl - b) - a/(Vl*(Vl + b)) - pvap
    
    return [eq1, eq2, eq3]

#Resolveremos com lest_squares, pois desta maneira podermos incluir algumas restrições.
#Usaremos as seguintes restrições:
# 1 atm < pvap < pc
# b < Vl < Vc
# Vc < Vv


res = least_squares(equations, (11, 100, 1789), bounds=((1, 65, 200 ),(30, 200, 3000 )) )

print('Solução com least_squares= ' , +  res.x)


# In[ ]:




