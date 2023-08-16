#!/usr/bin/env python
# coding: utf-8

# # Equilíbrio Químico Entre Gases Ideais

# ## Fundamentos teóricos
# 
# <font size="4">  Considere a reação $aA + bB \rightarrow cC + dD$ que ocorre em fase gasosa. Podemos escrever a energia livre de reação em função do quociente reacional
# 
# $\Delta G = \Delta Gº -RT lnQ$
# 
# e no equilíbrio a equação se reduz a $\Delta Gº = -RT lnK$, esta constante pode ser escrita em função das frações molares e concentrações molares da seguinte forma:
# 
# $K = K_x p^{\Delta \nu}$
# 
# $K = K_c (RT)^{\Delta \nu}$
# 
# Esta constante de equilíbrio depende de T segundo a equação de van’t Hoff: 
# 
# $ \frac{\partial lnK}{\partial 1/T} = \frac{\Delta_r Hº}{R T^2} $
# 
# $ ln \frac{K_2}{K_1} = -\frac{\Delta_r Hº}{R} \left (\frac{1}{T_2} – \frac{1}{T_1} \right )$ <font > 

# Para a resolução dos exercíos faremos uso das bibliotecas Numpy e Matplolib.

# In[1]:


#módulos para resolução dos exercícios:
import math
from scipy.optimize import fsolve


# # Exemplo 1

# <font size="4"> a) Use tabelas de dados termodinâmicos para calcular a constante de equilíbrio da reação $ N_2(g) + O_2(g) \rightarrow 2NO(g) $ em T = 298.15 K e 1 bar (0.987 atm). Calcule a concentração de equilíbrio de NO no ar atmosférico em b)  T = 298.15 K e 1 atm e c) T = 2000 K e 1 atm. <font> 

# ##  Solução:
# <font size="4"> Em primeiro lugar deve ser calculado a energia livre de reação. Em seguida a relação $ K = -RTln(\Delta_r G^0)$ para o cálculo da constante de equilíbrio a 298K. Para este tipo de reação, a expressão da cosntante de equilíbrio em termos das pressões parciais é igual a constante em termos de frações molares. Assim podemos assumir a composição aproximada do ar como $N_2$ = 0.78 e $O_2$ = 0.21 e que devido a magnitude da constante de equilíbrio estes valores não se alteram apreciavelmente. Será usada a relação $ K = K_x \times p^{\Delta \nu}$ para resolver o exercício. Por último usamos a equação $ ln(\frac{K_2}{K_1}) = -\frac{\Delta H^0}{R} \times (\frac{1}{T_2}-\frac{1}{T_1}) $para obter a constante a 2000 K e estimar a concentração de equilíbrio de NO.  <font>

# In[2]:


T = 298.15             #[K] temperatura
P = 0.987               #[atm] pressão
g_0_NO = 86.6           #[kJ/mol] G padrão de formação do NO
R = 8.314               #[J/(mol*K)]

g_0_O2 = 0.00           #[kJ/mol] G padrão de formação do O2    
g_0_N2 = 0.00           #[kJ/mol] G padrão de formação do N2    


# cálculo do DeltaG
delta_g_0 = 2*g_0_NO - g_0_O2 - g_0_N2          #[kJ/mol]
delta_g_01 = delta_g_0*1000                     #[J/mol]

# Cálculo de K
K_298 = math.exp((-delta_g_01)/(R*T))

#composição do ar atmosférico

y_N2 = 0.78
y_O2 = 0.21

y_NO_298 = math.sqrt(K_298*y_N2*y_O2)


print( " Constante de equilíbrio a 298.15 K é \t\t\t %.1e"%(K_298))
print(" A concentração de NO no equilíbrio a 298K é \t\t %.1e"%(y_NO_298))


# In[3]:


#item (c)

T_1 = 2000                  #[K]

h_0_NO = 90.25  # kJ/mol H padrão de formação do NO

delta_h_0 = 2*h_0_NO - 0 - 0          #[kJ/mol]
delta_h_01 = delta_h_0*1000           #[J/mol]


#cálculo K
K_2000 = math.exp(math.log(K_298) -delta_h_01/R *(1/T_1 -1/T))

y_NO_2000 = math.sqrt(K_2000*y_N2*y_O2)       
print( " Constante de equilíbrio a  2000 K é \t\t\t %.1e"%(K_2000))
print( " A fração molar de NO no equilíbrio a 2000 K é \t\t %.1e "%(y_NO_2000))


# ## Exemplo 2
# 
# <font size="4"> Refaça o cálculo da concentração de equilíbrio NO, do exercício anterior, a 2000 K e 1 atm, mas desta vez sem negligenciar a variação na concentração de $O_2$ e $N_2$.<font>

# ## Solução
# 
#  <font size="4">  Monte o quadro do equilíbrio e resolva a equação quadrática resultante. Uma alternativa, que será usada neste exemplo, é a minimização da energia livre de Gibbs de reação a partir das equações que descrevem o modelo. 
# 
# |   |   |   |   |
# |---|---|---|---|
# |  |$O_2$ |$N_2$ |2NO|
# |ini |0.21| 0.78 | 0 |
# |eq | 0.21-x | 0.78-x | 2x |    <font>

# In[4]:


Temp = 2000.                #[K]



K_2000 = 3.7e-4 

#escrevendo uma função no formato do python para cálculo do mínimo da energia livre 
# a função fsolve vai encontrar as raízes de uma função. Pode ser utilizada para resolução de sistemas. 

def f(x): 
	 return  (2*x)**2 - K_2000*(0.78-x)*(0.21-x)

x = fsolve(f,0.01) # encontrar o valor de x que zera a função


c_NO = 2*x


print(" Fração molar do NO calculado neste exemplo %.1e"%(c_NO))


# ## Exemplo 3
# 
# <font size="4">Determine a composição de equilíbrio para isomerisação do isobutano:
# 
# butano <=> isobutano  K = 4.52 <font >

# <font size="4"> Pra resolver o exercício basta escrever a constante de equilíbrio em termos da frção molar de isobutano
# 
# $$ K = \frac{x_{iso}}{1-x_{iso}}$$
# 
# Depois usamos o *fsolve* para minimizar a energia de Gibbs. <font >

# In[5]:


K = 4.52


def f(x_iso): 
	 return  x_iso/(1-x_iso)-K

x_iso = fsolve(f,0)

print(" Fração molar do isobutano %0.2f"%(x_iso))



# ## Exemplo 4
# 
# <font size="4"> Plantas modernas para produção de amônia trabalham em torno de 400ºC e 150 atm de pressão. Considerando estas condições, estime a conversão no equilíbrio. Assuma comportamento de gás ideal. K(400ºC, 1 bar) = 0.0137 <font >

# ## Solução
# 
# <font size="4"> A reação em questão é $ 0.5N_2(g) + 1.5H_2(g) \rightarrow NH_3(g) $.
# Inicialmente deve-se montar o quadro de equlíbrio e substituir as frações na equação da constante de equilíbrio
# 
# |          |          |          |          |          |
# |----------|----------|----------|----------|----------|
# |          |0.5$N_2$     |1.5$H_2$ |$NH_3$| |
# |ini       |0.5       | 1.5 | 0 |
# |eq        | 0.5-0.5x | 1.5-1.5x | x |
# |x         | $ \frac{0.5-0.5x}{2.0-x}$ | $ \frac{1.5-1.5x}{2.0-x}$ | $ \frac{x}{2.0-x } $ |
# 
# $ K_x = \frac{\frac{x}{2.0-x}}{(\frac{0.5-0.5x}{2.0-x})^{0.5} \times (\frac{1.5-1.5x}{2.0-x})^{1.5}}  $
# 
# Como foi fornecido a constante a 1 bar, é preciso recalcular o valor da constante de equilíbrio a 150 atm. Lembre-se que $K_p$ é independente da pressão.
# 
# $K_{400°C,1bar} = K_{x} \times P^{\Delta \nu_i}  = K_{x} \times P[bar]^{(1-0.5-1.5)}$
#     <font >

# In[6]:


Temp = 273.15+400           #[K]
P = 150*1.01325             #[bar]

K_673 = 0.013

#calculando a constante a 150 atm
K = K_673*P**(1.5+0.5-1)

#escrevendo a função que descreve o problema de equilíbrio

def f(x): 
	 return  (x/(2-x))/(((0.5-0.5*x)/(2-x))**(0.5)*((1.5-1.5*x)/(2-x))**(1.5))-K

x=fsolve(f,0.1)

y_NH3 = x/(2-x)

print("A fração molar de amônia no equilíbrio é  %0.2f"%(y_NH3))


# ## Exemplo 5
# <font size="4"> Gillespie e Beattie (Phys. Rev. 36:743–753 (1930)) computaram o parâmetro $K_{\phi}$ para síntese da amônia, representando o volume específico de vários gases pela equação de estado de Bettie-Bridgeman. Os resultados tem a seguinte forma:
# 
# $$ log(\frac{1}{K_{\phi}}) = (\frac{0.119184}{T}+\frac{91.87212}{T^2}+\frac{25122730}{T^4}) \times P $$ 
# 
# 
# Onde T e P são informados em K e atm, respectivamente.  Levando em consideração os desvios da idealidade, repita o cálculo da conversão no equilíbrio. <font >

# ## Solução
# 
# <font size="4"> Em altas pressões o modelo do gás ideal não pode ser usado. Logo, a constante de equilíbrio deve ser escrita em termos das fugacidades dos gases:
# 
# $$ K = \frac{f_{NH_3}}{f_{N_2} ^{0.5} f_{H_2} ^{1.5}} =  \frac{\phi_{NH_3}}{\phi_{N_2} ^{0.5} \phi_{H_2} ^{1.5}} \times \frac{p_{NH_3}}{p_{N_2} ^{0.5} p_{H_2} ^{1.5}} = K_{\phi} K_p$$
# 
# A partir da nova constante de equilíbrio é possível calcular a conversão no equilíbrio usando o mesmo quadro de equilíbrio do exercício anterior: 
# 
# |       |       |       |       |       |
# |:-----:|:-----:|:-----:|:-----:|:-----:|
# |  |0.5$N_2$ |1.5$H_2$ |$NH_3$| |
# |ini |0.5| 1.5 | 0 |
# |eq | 0.5-0.5x | 1.5-1.5x | x |
# |x | $ \frac{0.5-0.5x}{2.0-x}$ | $ \frac{1.5-1.5x}{2.0-x}$ | $ \frac{x}{2.0-x}  $  |
# 
# $ K_x = \frac{\frac{x}{2.0-x}}{(\frac{0.5-0.5x}{2.0-x})^{0.5} \times (\frac{1.5-1.5x}{2.0-x})^{1.5}}  $ <font>

# In[7]:


#Inicialização das variáveis
T = 273.15+400              #[K] 
P = 150                     #[atm]
P_bar = 150*1.01325 


K_673 = 0.013              

K_fi = (10**((0.1191849/T + 91.87212/T**2 + 25122730/T**4)*P))**(-1)

K = (K_673/K_fi)*P**(1.5+0.5-1)


def f(x): 
	 return  (x/(2-x))/(((0.5-0.5*x)/(2-x))**(0.5)*((1.5-1.5*x)/(2-x))**(1.5))- K
        
x = fsolve(f,0.5)

x_NH3 = x/(2-x)

print("A fração molar de amônia no equilíbrio é  %0.2f"%(x_NH3))


# In[ ]:




