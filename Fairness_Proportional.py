# -*- coding: utf-8 -*-
"""
OBIETTIVO: creare un modello di programmazione lineare che restituisca 
un'allocazione equa delle risorse secondo lo schema della equità proporzionale (PF)
"""

from pyomo.environ import *
import numpy as np

# Parse the training set
# def ParseData(filename):
#     fh = open(filename, 'r', encoding="utf-8")
#     Xs = []
#     Ys = []
#     for line in fh:
#         row = line.split('\t')        
#         Xs.append( list(map(float, row[:-1])) )
#         Ys.append( int(row[-1]) )
#     return Xs, Ys  

def ParseData(filename):
    doc = open(filename, 'r', encoding="utf-8")
    Ps = [] # presenze
    Ns = [] # nomi
    Ss = [] # surplus
    Fs = [] # consumo fisso
    Es = [] # copertura fotovoltaico
    for line in doc:
        row = line.split(',')        
        Ps.append( int(float(row[0])) )
        Ns.append( line[1] )
        Ss.append( int(line[8]) )
        Fs.append( int(line[9]) )
        Es.append( int(line[7]) )
    return Ps, Ns, Es, Ss, Fs 

def New(Ps, Ns, Es, kWh, Fs):
    """
    Parameters
    ----------
    Ps : lista : ogni entrata è binaria (1 se la fam è presente, 0 altrimenti)
    Ns : lista : ogni entrata è il nome della famiglia
    Es : lista : alla posizione i-esima ha la quantità di energia fornita dal FV alla fam i-esima
    Fs : lista : alla posizione i-esima ha il consumo fisso dell'appartamento i-esimo

    Returns
    -------
    Ns_new, Es_new, kWh_new, Fs_new : liste precedenti senza le famiglie in
        vacanza che hanno il consumo fisso completamente coperto dall'energia 
        che gli viene fornita dal FV normalmente
    RC : reale non negativo : totale nuova capacità a disposizione 
    
    """
    n = len(Ps)
    Ns_new = []
    Es_new = []
    kWh_new = []
    Fs_new = []
    rc = 0
    for i in range(1,n):
        if Ps[i-1] == 1: # se la fam è presente considero tutto
            Ns_new.append(Ns[i-1])
            Es_new.append(Es[i-1])
            kWh_new.append(kWh[i-1])
            Fs_new.append(Fs[i-1])
        else: # se la fam è assente
            diff = Es[i] - Fs[i] # energia scoperta dal FV
            if (diff) <= 0: # se il fotovoltaico non copre il consumo fisso, 
                            # allora kWh != 0
                Ns_new.append(Ns[i-1])
                Es_new.append(Es[i-1])
                kWh_new.append(-diff)
                Fs_new.append(Fs[i-1])
            else: # se il fotovoltaico copre il consumo fisso e avanza anche energia,
                  # allora quell'energia in più va in rc (a disposizione degli altri)
                rc += (diff)
    return Ns_new, Es_new, kWh_new, Fs_new,  rc

def PF(Ks, RC):
    """
    Risolve il problema lineare
    
    Parameters
    ----------
    Ks : lista : lista di surplus (output di New)
    prezzo: intero : prezzo dell'energia in €/kWh 
    RC : intero : massima capacità a disposizione
    
    Returns
    -------
    Ds : lista : alla posizione i-esima c'è la percentuale sul surplus totale che viene
        coperta dal fotovoltaico per l'appartamento i-esimo

    """    
    
    n = len(Ks) # numero dei dati (= numero appartamenti che hanno un surplus)

    # Build ILP Model
    model = ConcreteModel()
    
    # Indici
    model.N = RangeSet(n)

    # Variabili
    # model.u[i] è la percentuale di surplus della famiglia i-esima che viene 
    # coperta dal fv 
    model.u = Var(model.N, domain = NonNegativeReals)

    # Funzione obiettivo
    ### ?Linearizzare la somma di logaritmi?
    model.obj = Objective(expr = sum(np.log(model.u[i]) for i in model.N), 
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
    # 3. ogni appartamento deve avere una porzione non nulla di surplus coperto
    model.minind = ConstraintList()
    for i in model.N:
        model.minind.add( expr = model.u[i] > 0 )
    
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
    perc_covered = PF(kWh_new, RC) 
    
    ### PRINT
    kWh_covered = list(map(lambda x,y: x[i]*y[i]/100, perc_covered, kWh_new))
    kWh_uncovered = list(map(lambda x,y: x-y , kWh_new, kWh_covered))
    costi = list(np.multiply( kWh_uncovered, costo))
    print('lung Ns: {}, lung kWh_covered: {}, lung costi: {}'.format(len(Ns), len(kWh_covered), len(costi)))
    for i in range(1,l+1):
        print('Appartamento {}, kWh coperti: {}, Costo: {} \n'.format(Ns[i-1], kWh_covered[i-1], costi[i-1]))
