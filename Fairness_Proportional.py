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
    
    Parameters
    ----------
    Ps : TYPE
        DESCRIPTION.
    Fs : lista
        alla posizione i-esima ha il consumo fisso dell'appartamento i-esimo.

    Returns
    -------
    Ss_new, Rc_new, prezzi_new
    RC : reale non negativo
        totale nuova capacità a disposizione.
    
    """
    n = len(Ps)
    kWh_new = []
    rc = 0
    for i in range(1,n):
        if Ps[i-1] == 1:
            kWh_new.append(kWh[i-1])
        else:
            kWh_new.append(Fs[i-1])
            rc += (Es[i] - Fs[i])
    return kWh_new, rc

def MMF(Ks, RC, prezzi):
    """
    Risolve il problema lineare per trovare la soluzione ottima che ...

    Parameters
    Ss : lista
        lista di surplus (output di parsefile)
    Ps : lista
        lista di presenze (output di parsefile)
    prezzi : lista
        alla posizione i-esima c'è il prezzo in €/kWh per la famiglia i-esima
    RC : intero
        massima capacità a disposizione
    

    Returns
    -------
    Ds : lista
        alla posizione i-esima c'è la percentuale sul surplus totale che viene
        coperta dal fotovoltaico per l'appartamento i-esimo

    """    
    
    n = len(Ks) # numero dei dati (= numero di appartamenti)

    # Build ILP Model
    model = ConcreteModel()
    
    # Indici
    model.N = RangeSet(n)

    # Variabili
    model.u = Var(model.N, domain = NonNegativeReals)

    # Funzione obiettivo
    model.obj = Objective(expr = sum((model.u[i]) for i in model.N), 
                          sense = maximize)
    
    # Vincoli
    kWh_covered = []       
    for i in model.N:
        kWh_covered.append( Ks[i-1]*model.u[i]/100) # kWh coperti dal fotovoltaico
    # 1.Vincolo della massima capacità totale della risorsa
    model.maxtot = ConstraintList()   
    model.maxtot.add( expr = sum(kWh_covered) <= RC )
    # 2.Vincolo della massima capacità individuale della risorsa
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
        
    return perc_covered, kWh_covered


# -----------------------------------------------
#   MAIN function
# -----------------------------------------------
if __name__ == "__main__":
    
    ### DATI
    # # Ns lista dei nomi delle famiglie
    # # Ss lista dei consumi-surplus di ogni famiglia
    # # Ps lista delle presenze delle famiglie Ps[i]=0,1
    # # prezzi è la lista dei prezzi in €/kWh di ciascuna famiglia
    Ps, Ns, Es, kWh, Fs  = ParseData('condominio2.csv')
    print('Presenze: {}'.format(Ps))
    print('Nomi: {}'.format(Ns))
    print('Dal fotovoltaico: {}'.format(Fs))
    print('Surplus: {}'.format(kWh))
    print('Fisso: {}'.format(Fs))
    # Ns = ['Bianchi','Rossi','Verdi','Smith','Otto','Nove']
    # kWh = [1,1,1,2,2,3] # consumo surplus
    # Ps = [1,1,0,1,0,1] # presenze
    # Fs = [0.5,0.5,0.5,0.5,1,1] # consumo fisso
    # Es = [0.8,0.8,0.8,0.9,1,1.5] # energia di base coperta dal fotovoltaico
    # costo = 0.277 # 0,277€/kWh
    kWh_new, RC = New(kWh, Fs, Es, Ps) # lista tc l'i-esima entrata ha il consumo surplus
                            # dell'appartamento i-esimo
    
    # ### SOLUZIONE
    # perc_covered, kWh_covered = MMF(kWh_new, RC, costo) # lista tc l'i-esima entrata ha la percentuale  
    #                        # della i-esima richiesta che rappresenta la copertura 
    #                        # dal fotovoltaico
    
    # ### PRINT
    # # kWh_covered = list(map(lambda x,y: x[i]*y[i]/100, perc_covered, kWh_new))
    # kWh_uncovered = list(map(lambda x,y: x-y , kWh_new, kWh_covered))
    # costi = list(np.multiply( kWh_uncovered, costo))
    # print('lung Ns: {}, lung kWh_covered: {}, lung costi: {}'.format(len(Ns), len(kWh_covered), len(costi)))
    # for i in range(1,6):
    #     print('Appartamento {}, kWh coperti: {}, Costo: {} \n'.format(Ns[i-1], kWh_covered[i-1], costi[i-1]))
