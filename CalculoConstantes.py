import pandas as pd
import numpy as np

# Esta función recibe como parámetros la temperatura de una reacción, sus
# coeficientes estequiométricos (positivo si es producto, negativo si es reactivo),
# y los datos de sus especies a partir de una base de datos y devuelve el valor de
# la constante de equilibrio, como logaritmo y como logaritmo natural. El valor de
# la temperatura se introduce en Kelvin.
# Para calcularlo se considera la ecuación de Van't Hoff, las entalpías de formación y
# sus energías libres respectivas. Por el momento, solamente se consideran las especies 
# gaseosas. Si encuentran forma de considerar las otras especies, me avisan. 

def ConstanteQuimica(T, Species, Coefficients):
    R=8.31447 
    if len(Coefficients)!=len(Species):
        raise Exception('La cantidad de especies no coincide con el número de coeficientes estequiométricos. Compruebe de nuevo')
    df = pd.read_csv('DataPure.csv',index_col=0)
    # Mientras encuentre la forma de arreglar esto, las únicas sustancias reconocidas en el 
    # programa serán gases
    df = df[(df['State'] == 'g')]    
    df[['A','B','C','D']].fillna(0)

    # Sonlas constantes del polinomio por especie, proceden a sumarse para formar el polinomio
    ABCD = pd.DataFrame.multiply(df.loc[Species,['A','B','C','D','Enthalpy','Gibbs']],Coefficients,axis='rows')  
    
    # Se multiplican las constantes por los valores reportados en el Smith Van Ness
    PolyH = np.multiply(ABCD.sum(axis='rows').values,[1,1e-3,1e-6,1e5,1,1])
    # La derivación de estas fórmulas se encuentra en el documento titulado "Readme"
    Lambda = -298.15*PolyH[0]-(298.15**2)*PolyH[1]/2-(PolyH[2]/3)*298.15**3+PolyH[3]/298.15
    lnk = ((Lambda+PolyH[4]/R)*(1/298.15-1/T)-PolyH[5]/(R*298.15)+PolyH[0]*np.log(T/298.15)+(PolyH[1]/2)*(T-298.15)+(PolyH[2]/6)*(T**2-298.15**2)+(PolyH[3]/2)*(1/T**2-1/298.15**2))
    
    lnk_indep = -PolyH[5]/(R*298.15) + PolyH[4]/R*(1/298.15-1/T)
    k_indep = np.exp(lnk_indep)
    pk_indep = -np.log10(k_indep)

    k = np.exp(lnk)
    pk = -np.log10(k)
    
    print('Asumiendo que la constante depende de la temperatura\n')
    print('ln K = ', lnk)
    print('K = ', k)
    print('pK = ', pk)
    print('\nAsumiendo que la constante es independiente de la temperatura:\n')
    print('ln K independiente = ', lnk_indep)
    print('K independiente = ', k_indep)
    print('pK independiente = ', pk_indep)

def SpHeat(t1, t2, Species, Coefficients=1):
    R=8.31447 

    if len(Coefficients)!=len(Species):
        raise Exception('La cantidad de especies no coincide con el número de coeficientes estequiométricos. Compruebe de nuevo')

    df = pd.read_csv('DataPure.csv',index_col=0)
    # Mientras encuentre la forma de arreglar esto, las únicas sustancias reconocidas en el 
    # programa serán gases
    df = df[(df['State'] == 'g')]    
    df[['A','B','C','D']].fillna(0)
    ABCD = pd.DataFrame.multiply(df.loc[Species,['A','B','C','D']],Coefficients,axis='rows')  

    # Ajustando los valores de acuerdo al Smith Van Ness:
    ABCD = np.multiply(ABCD.sum(axis='rows').values,[1,1e-3,1e-6,1e5])
    # La siguiente es la parte constante
    Lambda = -t1*ABCD[0]-(t1**2)*ABCD[1]/2-(ABCD[2]/3)*t1**3+ABCD[3]/t1
    deltaH = (Lambda + ABCD[0]*t2 + ABCD[1]*t2**2/2 + ABCD[2]*t2**3/3 - ABCD[3]/t2)*R
    print("El calor es: {} J".format(str(deltaH)))
    return deltaH

ConstanteQuimica(300+273.15,['Acetylene','Nitrogen','Hydrogen cyanide'],[-1,-1,2])

################################################################################################
##--------------------------------------------------------------------------------------------##
##--------------------------------------Código muerto-----------------------------------------##
##--------------------------------------------------------------------------------------------##
################################################################################################


#pd.set_option("display.max_rows",9,"display.max_columns",10)

#df = pd.read_csv('DataPure.csv',index_col=0)
# Mientras encuentre la forma de arreglar esto, las únicas sustancias reconocidas en el programa
# # serán gases
#df = df[(df['State']=='g')]    
#df.update(df[['A','B','C','D']].fillna(0))
#print(df)
#w = pd.DataFrame.multiply(df.loc[['Methane','Water'],['A','B','C','D']],[1,2,],axis='rows')

#print(w.sum(axis='rows'))
#print(df.loc[['Water','Methane'],['A','B','C','D']])

# La siguiente es la solución del problema 2 de la tarea:
# A 145 ºC:
# ConstanteQuimica(145+273.15,['Ethylene','Water','Ethanol'],[-1,-1,1])

# A 320 ºC
# ConstanteQuimica(320+273.15,['Ethylene','Water','Ethanol'],[-1,-1,1])
# SpHeat(298,850+273.15,['Water','Carbon dioxide','Carbon monoxide','Hydrogen'],[0.1725,0.0275,0.1725,0.6275])

