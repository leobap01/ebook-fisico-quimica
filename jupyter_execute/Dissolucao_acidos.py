#!/usr/bin/env python
# coding: utf-8

# # Equilíbrio de fases em misturas
# ## Dissolução de gases ácidos em água

# In[1]:


#módulos para resolução dos exercícios:
import numpy as np
from scipy.optimize import fsolve


# <font size="4"> A dissolução de espécies gasosas em solução pode ser descrita a partir do conhecimento da constante da Lei de Henry de cada substância. Uma aplicação típica da Lei de Henry é estudar a dissolução de gases encontrados na atmosfera em góticulas do aerossol atmosférico ou em corpos hídricos. Por exemplo, é possível estimar o pH médio do aerossol atmosférico (e consequentemente  de nuvens) se soubermos a concentração de $CO_2$ na atmosfera. Sabendo que a concentração média de $CO_2$ na atmosfera é 360 ppm (mol/mol) e a constante da Lei de Henry para dissolução em água é $3,4 \times 10^{-2} ~mol~ L^{-1}~ atm^{-1}$, estime o pH médio da chuva considerando dois casos: a) Considerando que apenas o íon bicarbonato, oriundo da dissolução do $CO_2$, afeta o pH; b) Considerando que os íons carbonato e bicarbonato são importantes para determinar o pH da chuva. </font> 

# ## Solução
# 
# <font size="4"> Os três equilíbrios que precisam ser estudados são:
# $$ CO_{2(g)} \rightarrow CO_{2(aq)} ~~~ K_H = \frac{[CO_2]}{pCO_2} = 3.4 \times 10^{-2}$$
# $$ CO_{2(aq)} \rightarrow HCO_{3(aq)} ^{-} + H_(aq) ^{+} ~~~ K_{a1} = \frac{[HCO_3 ^{-}][H^{+}]}{[CO_2]} = 4.5 \times 10^{-7}$$
# $$ HCO_{3(aq)} \rightarrow CO_{3(aq)} ^{2-} + H_(aq) ^{+} ~~~ K_{a2} = \frac{[CO_3 ^{2-}][H^{+}]}{[HCO_3 ^{-}]} = 7.0 \times 10^{-11}$$
# 
# Os valores de Ka foram obtidos na literatura. Sistema a ser resolvido assumindo inicialmente que o segundo equilíbrio pode ser negligenciado
# 
# $$ K_H \times pCO_2 - [CO_2] = 0 $$
# $$ K_{a1} \times [CO_2] - [HCO_3 ^{-}][H^{+}] = 0 $$
# $$ [HCO_3 ^{-}] - [H^{+}] = 0 $$ 
# 
# A condição de eletroneutralidade foi usada para fechar o sistema de equações </font> 

# In[2]:


#Definição das variáveis

Kh = 3.4e-2 # mol/L atm
Ka1 = 4.5e-7 
Ka2 = 7.0e-11


# pressão de CO2 é obtida a partir da razão de mistura informada. Considerado que o sistema está a 1 atm, 
# 360 ppm de CO2 equivalem a 360E-6 atm de CO2 na atmosfera.
pCO2 = 360e-6 # atm (pressão parcial de CO2 no vapor)

#definição do sistema de equações

def equations(vars):
    CO2, HCO3, H = vars
    
    eq1 = Kh*pCO2 -CO2 
    eq2 = Ka1*CO2 - HCO3*H
    eq3 = HCO3 - H
    
    return [eq1, eq2, eq3]

# resolução do sistema com a estimativa inicial das concetrações em mol/L
CO2, HCO3, H =  fsolve(equations, (1e-2, 1e-3,1e-3))

print('Solução com fsolve (mol/L)= ' + format(H , ' 6.5e'))

pH = -np.log10(H)
print('pH = ', pH)


# In[3]:


#Resolução da parte b incluindo o terceiro equilíbrio
# Note que todas variáveis foram definidas na primeira célula

def equations(vars):
    CO2, HCO3, H, CO3 = vars
    
    eq1 = Kh*pCO2 -CO2 
    eq2 = Ka1*CO2 - HCO3*H
    eq3 = Ka2*HCO3 - CO3*H
    eq4 = HCO3 + CO3 - H
    
    return [eq1, eq2, eq3, eq4]

CO2, HCO3, H, CO3 =  fsolve(equations, (1e-2, 1e-6,1e-6, 1e-10 ))




print('[H+] (mol/L)= '  + format(H , ' 6.5e'))

#Podemos imprimir a concentração de carbonato para verificar que o terceiro equilíbrio é desprezível

print('[CO3] (mol/L)= '  + format(CO3 , ' 6.5e'))

pH = -np.log10(H)
print('pH = ', pH)


# Podemos verificar que o pH em ambos os casos é praticamente igual e que a concentração de carbonato é muito baixa, confirmando que o terceiro equilíbrio pode ser desprezado.

# <font size="4"> Doyle e colaboradores (Environmental Science & Technology, Volume 13, Number 11, November 1979, p 1416) determinaram a concentração de $HNO_3$ na atmosfera da cidade de Riverside, Califórnia, usando espectrometria na região do infravermelho. Os valores máximos de $HNO_3$ obtidos em cada dia estão apresentados na tabela abaixo.
# 
# |   |   |   |   |
# |---|---|---|---|
# |Dia | 21/07/1977 | 25/07/1977 | 11/08/1977
# |$HNO_3$ (ppb (mol/mol) | 10 | 20 | 13
# 
# Sabendo que a constante da Lei de Henry para solubilização do $HNO_3$ em água é $2,1 \times 10^5 mol ~L^{-1} ~atm^{-1}$ e que o Ka deste ácido é 20, repita o cálculo anterior e calcule o pH da chuva na presença destas concentração de ácido nítrico. DICA: Não podemos disconsiderar a contribuição do $CO_2$. </font>

# In[ ]:





# In[ ]:





# In[4]:


#teste

Ka = 20
Kh_hno3 = 2.1e5
pHNO3 = 10e-9


def equations(vars):
    CO2, HCO3, H, HNO3, NO3 = vars
    eq1 = Kh*pCO2 -CO2 
    eq2 = Ka1*CO2 - HCO3*H
    eq3 = Ka*HNO3 - NO3*H
    eq4 = Kh_hno3*pHNO3 - HNO3
    eq5 = HCO3 + NO3 - H
    return [eq1, eq2, eq3, eq4, eq5]

CO2, HCO3, H, HNO3,NO3 =  fsolve(equations, (1e-2, 1e-6,1e-3, 1e-20,1e-3 ))
print('Solução com fsolve= ' +  format(CO2 , ' 6.5e') + format(H , ' 6.5e'))
print('concentração de nitrato=', NO3, 'mol/l' )

pH = -np.log10(H)
print('pH = ', pH)


# In[ ]:




