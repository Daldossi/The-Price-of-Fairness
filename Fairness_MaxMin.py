# -*- coding: utf-8 -*-
"""
OBIETTIVO: creare un modello di programmazione lineare che restituisca 
un'allocazione equa delle risorse secondo lo schema della equità maxmin (MMF)
"""

from pyomo.environ import *
import numpy as np

# Parse the training set
def ParseData(filename):
    doc = open(filename, 'r', encoding="utf-8")
    Ns = [] # nomi
    Ps = [] # presenze
    Ss = [] # surplus
    Fs = [] # consumo fisso
    Es = [] # copertura fotovoltaico
    for line in doc:
        row = line.split(',')        
        Ns.append( line[1] )
        Ps.append( int(float(row[0])) )
        Ss.append( float(line[8]) )
        Fs.append( float(line[9]) )
        Es.append( float(line[7]) )
    return Ps, Ns, Es, Ss, Fs 

def New(kWh, Fs, Es, Ps):
    """
    Sarà uguale a quella nel codice di Fairness_Proportional    
    """
    return kWh_new, rc

def MMF(Ks, RC):
    """
    max-min scheme
    """    
    
    n = len(Ks) # numero dei dati (= numero di appartamenti)

    # Build ILP Model
    model = ConcreteModel()
    
    # Indici
    model.N = RangeSet(n)

    # Variabili
    model.u = Var(model.N, domain = NonNegativeReals)

    # Funzione obiettivo
    ### ?Devo massimizzare lessicograficamente?
    model.obj = Objective(expr = sum((model.u[i]) for i in model.N), 
                          sense = maximize)
    
    # Vincoli
    kWh_covered = []       
    for i in model.N:
        kWh_covered.append( Ks[i-1]*model.u[i]/100) # kWh coperti dal fotovoltaico
    # 1. massima disponiiblità totale della risorsa
    model.maxtot = ConstraintList()   
    model.maxtot.add( expr = sum(kWh_covered) <= RC )
    # 2. massima copertura individuale della risorsa
    model.maxind = ConstraintList()            
    for i in model.N:
        model.maxind.add( expr = kWh_covered[i-1] <= Ks[i-1] )
    
    # Risoluzione con gurobi
    solver = SolverFactory('gurobi')
    sol = solver.solve(model, tee=False)

    # check feasibility
    sol_json = sol.json_repn()
    # Check solution status
    if sol_json['Solver'][0]['Status'] != 'ok':
        return None    
    
    perc_covered = model.u
    
    # dict_gironi = {}
    # for k in model.K:
    #     girone_k = []
    #     for i in model.N:
    #         if model.y[i,k]() > 0.5:
    #             girone_k.append( (Ls[i-1][0], i) )
    #     dict_gironi['Girone' + str(k)] = girone_k
        
    return perc_covered


# -----------------------------------------------
#   MAIN function
# -----------------------------------------------
if __name__ == "__main__":
    
    ### DATI
    ### Ps : lista delle presenze delle famiglie Ps[i]=0,1
    ### Ns_old : lista dei nomi delle famiglie
    ### Es_old : lista dell'energia dal fotovoltaico (fv)
    ### kWh_old : lista dei consumi-surplus di ogni famiglia
    ### Fs_old : lista dei consumi fissi (l'appartamento della famiglia i che è
    ###          in vacanza consuma Fs[i])
    ### costo : prezzo in €/kWh dell'energia
    ## Dal documento
    Ps, Ns_old, Es_old, kWh_old, Fs_old  = ParseData('condominio2.txt')
    print('Presenze: {}'.format(Ps))
    print('Nomi: {}'.format(Ns))
    print('Dal fotovoltaico: {}'.format(Fs))
    print('Surplus: {}'.format(kWh))
    print('Fisso: {}'.format(Fs))
    ## A mano
    # Ns = ['Bianchi','Rossi','Verdi','Smith','Otto','Nove']
    # kWh = [1,1,1,2,2,3] # consumo surplus
    # Ps = [1,1,0,1,0,1] # presenze
    # Fs = [0.5,0.5,0.5,0.5,1,1] # consumo fisso
    # Es = [0.8,0.8,0.8,0.9,1,1.5] # energia di base coperta dal fotovoltaico
    ##
    costo = 0.277 
    
    Ns_new, Es_new, kWh_new, Fs_new, RC = New(Ps, Ns_old, Es_old, kWh_old, Fs_old) 
    l = len(Ns_new)
    
    ### SOLUZIONE
    perc_covered = MMF(kWh_new, RC) 
    
    ### PRINT
    kWh_covered = list(map(lambda x,y: x[i]*y[i]/100, perc_covered, kWh_new))
    kWh_uncovered = list(map(lambda x,y: x-y , kWh_new, kWh_covered))
    costi = list(np.multiply( kWh_uncovered, costo))
    print('lung Ns: {}, lung kWh_covered: {}, lung costi: {}'.format(len(Ns), len(kWh_covered), len(costi)))
    for i in range(1,l+1):
        print('Appartamento {}, kWh coperti: {}, Costo: {} \n'.format(Ns[i-1], kWh_covered[i-1], costi[i-1]))
