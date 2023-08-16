#!/usr/bin/env python
# coding: utf-8

# # Equilíbrio de fases em misturas

# ## Fundamentos Teóricos
# 
# Quando um soluto é adicionado em um solvente, este leva ao abaixamento do potencial químico do solvente, afetando assim as propriedades físicas do solvente. Como primeira aproximação, podemos assumir a formação de uma solução ideal. Neste caso, $ \Delta_{mis}G$ e $\Delta_{mis}S$ podem ser calculados pelas equações:
# 
# $ \Delta_{mis}G = RT \sum_{i} n_i lnx_i $
# 
# $ \Delta_{mis}S = -R \sum_{i} n_i lnx_i $
# 
# De forma que em uma solução ideal $ \Delta_{mis}H = 0$  e $ \Delta_{mis}V = 0 $
# 
# Em uma solução ideal soluto e solvente obedecerão a Lei de Raoult, e a pressão de vapor de cada espécie que faz parte da solução dependerá apenas da composição e pressão de vapor dos componentes puros:
# 
# $ p_i = x_i p_i^* $     Lei de Raoult
# 
# onde $p_i^* $ é a pressã ode vapor da espécie pura.
# 
# 
# No entanto poucos sistemas formam soluções ideais. No entanto, em soluções muito diluídas verifica-se que solvente continua a obedecer a Lei de Raoult. Em contrapartida o soluto, por estar em um ambiente químico em que interage majoritariamente com moléculas do solvente, terá sua pressão de vapor calculada pela Lei de Henry:
# 
# $ p_A = x_A p_A^* $     Lei de Raoult
# $ p_i = K_i p_i^* $     Lei de Henry
# 
# onde $K_i$ é a constante de Lei de henry que é função do solvente e temperatura. Esta constante substitui a pressã ode vapor do soluto como constante de proposrcionalidade na Lei de Raoult. Uma solução com estas características é chama de solução diluída ideal.
# 
# 

# A resolução dos exercícios a seguir só necessitam do módulo Numpy. 

# In[1]:


import numpy as np


# ## Exemplo 1
# Estime a fração molar de vapor d'água no ar atmosférico em equilíbrio com água líquida a 20ºC e 1 atm. 

# ## Solução 
# 
# A resolução deste exercídio fará uso da Lei de Raoult: 
# 
# $$ p_{H_2 O} = x_{H_2 O}p_{vap} ^*$$
# $$ y_{H_2 O}p_{total} = x_{H_2 O}p_{vap} ^*$$
# 
# Nestas condições: $p_{vap} ^* = 0.023 atm$
# Assume-se que o ar é composto de oxigênio e nitrogênio e que estes são pouco solúveis em água: $x_{H_2 O} = 1$

# In[2]:


P = 1.          #[atm]
p_vap = 0.023     #[atm]
x_water = 1

y_water = x_water*p_vap/P

print( "A fração molar de vapor d'água em equilíbrio com líquido é %f"%(y_water))

# umidade relativa

Um = (y_water * P / p_vap)*100
print( "Umidade relativa nestas condições  %f"%(Um))


# ## Exemplo 2
# Estime a concentração de oxigênio dissolvido em água quando água e ar atmosférico estão em equilíbrio a 20 ºC e 1 atm. Utilize os resultados do exemplo anterior. Dados: $ \rho (H_2O)= 998.2  g/L $ e  $H(O_2) = 40100 atm $ a 1 atm e 20ºC

# ## Solução
# 
# De acordo com o exercício anterior, $y_{(H_2 O)} = 0.023$. O ar seco é composto por 21% de $O_2$ e 79% de $N_2$, logo, $ y_{(N_2 )} +  y_{(O_2 )} = 1 - y_{(H_2 O)}$. Devemos assumir que a solução aquosa é majoritariamente água, devido a baixa solubilidade dos gases.

# In[3]:


#Inicialização ads variáveis do problema

P = 1.          #[atm]
H_O2 = 40100    #[atm] Constante da Lei de Henry do O2
y_agua = 0.023

#Calculando a fração molar de O2 no vapor
y = 1-y_agua
y_O2 = y*0.21

x_O2 = y_O2*P/H_O2 # fração molar de O2 dissolvido em água

c = x_O2*998.2/18   #[(mol O2)/(L solução)]

#Assumindo comportamento de gás ideal para o O2
V = c*0.082*293         #[(L O2)/(L solução)]
V = V*1000              #[(ml O2)/(L solução)]

print( "Concentração de oxigênio dissolvido em água no equilíbrio é %f mL O2)/L solução)"%(V))


# ## Exemplo 3
# 
# Repita o primeiro exercício, mas levando em consideração os três equilíbrios existentes simultaneamente:a) Equilíbrio líquido-vapor da água, b) Dissolução de $O_2$ do ar em água e c) Dissolução de $N_2$ do ar em água. Dado:  $H(N_2) = 80400 $ atm a 1 atm e 20ºC.

# ## Solução
# 
# Para resolver o problema devemos montar um sistema de equações lineares que descreve o problema. Os três equilíbrios existentes são:
# 
# $$ y_{agua}P = x_{agua}p_{vap}(H_2 O) $$
# $$ y_{O_2}P = x_{O_2}H_{O_2} $$
# $$ y_{N_2}P = x_{N_2}H_{N_2} $$
# 
# Temos um conjunto de 3 equações com 6 incógnitas. As três equações que completam o sistema são:
# 
# $$ y_{agua} + y_{O_2} + y_{N_2} = 1 $$
# $$ x_{agua} + x_{O_2} + x_{N_2} = 1 $$
# $$ \frac{y_{O_2}}{y_{N_2}} = 0.21/0.79 = 0.266 $$
# 
# A resolução deste sistema pode ser feita pela regra de Cramer: $ A = A^{-1} B $. Escrevendo as equações acima no formato de um sistema de equações, temos:
# 
# $$ x_{agua}p_{vap}(H_2 O) + 0x_{O_2} + 0x_{N_2} - y_{agua}P + 0y_{O_2} + 0y_{N_2}= 0$$
# $$ 0x_{agua} + x_{O_2}H_{O_2} + 0x_{N_2} + 0y_{agua} -y_{O_2}P + 0y_{N_2} = 0 $$
# $$ 0x_{agua} + 0x_{O_2} + x_{N_2}H_{N_2} + 0y_{agua} +0y_{O_2}P - y_{N_2}P = 0 $$
# $$ 0x_{agua} + 0x_{O_2} + 0x_{N_2} + y_{agua} + y_{O_2} + y_{N_2} = 1 $$
# $$ x_{agua} + x_{O_2} + x_{N_2} + 0y_{agua} + 0y_{O_2} + 0y_{N_2}= 1 $$
# $$ 0x_{agua} + 0x_{O_2} + 0x_{N_2} + 0y_{agua} + 079y_{O_2} - 0.21y_{N_2} = 0 $$
# 

# In[4]:


#Inicializando as variáveis:

P = 1.0             #[atm]
p_agua = 0.023      #[atm] Pressão de vapor da água
H_o = 40100         #[atm] Constante da Lei de Henry
H_n = 80400.        #[atm] Constante da Lei de Henry




A = np.array([[0.023, 0, 0, -1, 0, 0],[0, 40100, 0, 0, -1, 0],[0, 0 ,80400, 0, 0, -1],[0, 0, 0, 1, 1 ,1],[1 ,1, 1, 0, 0 ,0],[0, 0, 0, 0, 0.79, -0.21]])

B = np.array([[0],[0],[0],[1],[1],[0]])

#X = np.linalg.inv(A)

#C = np.dot(X, B) # produto das duas matrizes

C = np.linalg.solve(A, B)



print( " A composição nas fases líquidas e vapor são :")

print( "    y_agua       \t %f"%(C[3]))
print( "    y_O2         \t %f"%(C[4]))
print("    y_N2          \t %f"%(C[5]))
print( "    x_agua       \t %f"%(C[0]))
print( "    x_O2         \t %e"%(C[1]))
print( "    x_N2         \t %e"%(C[2]))


# ## Exemplo 4
# Estime a pressão de vapor e a composição do vapor em equilíbrio com o líquido de composição 80%, em mol, de benzeno e 20% de tolueno a 20ºC, assumindo que a solução é ideal. Use a equação de Antoine para calcular a pressão de vapor das substâncias puras.
# 

# ## Solução
# Uma solução ideal é aquela que solvente e soluto obedecem a Lei de Raoult. Para isso é necessário saber a pressão de vapor das substâncias puras, que será cálculada pela equação de Antoine: $ log(p(torr)) = A - \frac{B}{T(ºC)+C} $.
# Os parâmetros A, B e C são tabelados para cada substância.
# 
# Benzeno:
# A = 6.90565
# B = 1211.033
# C = 220.79
# 
# Tolueno:
# A = 6.95334
# B = 1343.943
# C = 219.377

# In[5]:


#Inicializando as variáveis
T = 20.         #[C]
x_b = 0.80
x_t = 0.20
A_b = 6.90565
B_b = 1211.033
C_b = 220.79

p_b = 10**(A_b-B_b/(T+C_b)) #pressão de vapor do benzeno em torr


A_t = 6.95334
B_t = 1343.943
C_t = 219.337

p_t = 10**(A_t-B_t/(T+C_t)) # pressão de vapor do tolueno em torr

p_1 = x_b*p_b # pressão parcial do benzeno no vapor
p_2 = x_t*p_t # pressão parcial do tolueno no vapor


P = p_1+p_2 #cálculo da pressão total

#composição do vapor
y_b = x_b*p_b/P
y_t = x_t*p_t/P



print( " Pressão da mistura na fase vapor         %f torr"%(P))
print( " Fração molar do benzeno na fase vapor    %f"%(y_b))
print( " Fração molar do tolueno na fase vapor    %f"%(y_t))


# ## Exemplo 5
# 
# Em qual temperatura a mistura anterior de benzeno e tolueno terá pressão de vapor de uma atmosfera? Assuma qua a solução é ideal.
# 

# ## Solução
# 
# Podemos escrever a lei de Raoult para mistura:
# 
# $$ (y_{benzeno} + y_{tolueno}) P = x_{benzeno}p_{vap}(b) + x_{tolueno}p_{vap}(t) $$
# 
# É possível calcular a composição do vapor que que tem pressão de 1 atm, a partir do valor da pressão de vapor de cada substância pura, pois esta é função da temperatura. No entanto, as equações usadas para o cálculo da pressão de vapor são equações logarítmicas e equações contendo dois logaritmos diferentes não tem solução analítica. Logo, este é um problema de tentativa e erro, que é resolvido facilmente com auxílio de computadores. 

# In[6]:


P = 760.            #[mm Hg]
x_b = 0.8           # Fração molar de benzeno no líquido
x_t = 0.2           # Fração molar de tolueno no líquido

#parâmetros da equação de Antoine para benzeno
A_b = 6.90565
B_b = 1211.003
C_b = 220.79

#parâmetros da equação de Antoine para tolueno
A_t = 6.95334
B_t = 1343.943
C_t = 219.337



#vamos criar um loop para ir testando diferentes valores de temperatura
#Inicialmente criaremos um parâmetro de erro que será testado a cada iteração

err = 1.
T = 50.              #[C] temperatura inicial


while err > 10**(-3):
    p_b = 10**(6.90565 - 1211.003/(T + 220.79))
    p_t = 10**(6.95334 - 1343.943/(T + 219.337))
    y_b = x_b*p_b/P
    y_t = x_t*p_t/P
    err = abs((y_b + y_t) - 1) # este valor deveria ser zero, mas o tempo de computação seria muito longo
    T = T + 0.01 # se escolher mal este incremento terá um erro considerável ou não conseguirá encontrar a solução.

print( "A temperatura no qual a mistura benzeno-tolueno tem pressão de vapor de 1 atm é %0.3f deg C"%(T))


# ## Exemplo 6
# 
# O sangue de um mergulhador está saturado com nitrogênio a 5 atm de pressão e fração molar de 0.79. Assim que o mergulhador se direciona a superfície, seu sangue precisa expelir o nitrogênio a medida que seu corpo entra em equlíbrio com a atmosfera. Quanto nitrogênio o corpo deve expelir? Assuma que o mergulhador pese 55 kg, que 75% de seu corpo seja água, com temperatura corporal de 37 ºC, e que a solubilidade do nitrogênio nos fluidos corporais é a mesma que em água pura. Dados: $ H(N_2) = 10.05 \times 10^4$ atm a 37ºC.

# ## Solução
# 
# O número de mols de nitrogênio expelido do corpo do mergulhador é  proporcional diferença entre a fração molar de nitrogênio dissolvido a 5 atm e a fração molar a 1 atm. É necessário calcular a quantidade de matéria dos fluidos corporais e multiplicá-lo pela fração de nitrogênio expelida pelo corpo para encontrar a quantidade total de nitrogênio expelido.
# 
# $$ N_2(expelido) = n_{fluidos} (x_{N_2, 5atm} - x_{N_2, 1atm} )$$
# 
# Consideraremos que a fração molar de $N_2$ (0.79), não é alterada com a variação e pressão.
# 

# In[7]:


P_1 = 5.                #[atm]
y_n = 0.79              # Fração molar de N2 na atmosfera
P_2 = 1.0               #[atm]
M = 55.                 #[kg] Massa do mergulhador
x_w = 0.75              # fração de água no corpo humano
T = 37                  #[C] Temperatura corporal do mergulhador 

H_n = 10.05*10**4       # [atm]

n_rej = (M*1000*x_w/18)*( P_1*y_n/H_n - P_2*y_n/H_n)         #[mol]

M_rej = n_rej * 28

V_n = n_rej*0.082*293/1            #[L]

print(" Massa de nirogênio expelida pelo corpo %0.2f g"%(M_rej))
print(" Volume de nirogênio expelida pelo corpo %0.2f L"%(V_n))


# In[ ]:




