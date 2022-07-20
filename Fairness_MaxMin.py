# -*- coding: utf-8 -*-
"""
OBIETTIVO: creare un modello di programmazione lineare che restituisca 
un'allocazione equa delle risorse secondo lo schema della equità maxmin (MMF)
"""

from pyomo.environ import *
import numpy as np

from itertools import combinations

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
    rc = 0
    for i in range(0,n):
        if Ps[i] == 1: # se la fam è presente considero tutto
            Ns_new.append(Ns[i])
            Es_new.append(Es[i])
            kWh_new.append(kWh[i])
        else: # se la fam è assente
            diff = Es[i] - Fs[i] # energia scoperta dal FV
            if (diff) <= 0: # se il fotovoltaico non copre il consumo fisso, 
                            # allora kWh != 0
                Ns_new.append(Ns[i])
                Es_new.append(Es[i])
                kWh_new.append(-diff)
            else: # se il fotovoltaico copre il consumo fisso e avanza anche energia,
                  # allora quell'energia in più va in rc (a disposizione degli altri)
                  rc += (Es[i] - Fs[i])
    return Ns_new, Es_new, kWh_new, rc

def MMF(Ks, RC):
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
    
    n = len(Ks) # numero dei dati (= numero di appartamenti che hanno un surplus)

    # Build ILP Model
    model = ConcreteModel()
    
    # Indici
    model.N = RangeSet(n)

    # Variabili
    # model.u[i] è la percentuale di surplus della famiglia i-esima che viene 
    # coperta dal fv
    model.u = Var(model.N, domain = NonNegativeReals)

    # Funzione obiettivo
    ### ?Devo massimizzare lessicograficamente? 
    ### Multiple Objectives: use n=6 objectives with a lexicographic approach:
    ### 1. massimizzo la percentuale che possono avere tutti,
    ### 2. massimizzo la percentuale che hanno tutti tranne uno,
    ### 3. ... così via
    ### 6. massimizzo la percentuale che ha una sola famiglia,
    ### Allo step k-esimo non posso abbassare il valore dei precedenti k-1 step.
    ### Ogni volta cambio la funzione obiettivo, tenendo le precedenti come vincoli
    ### Duale?
    model.obj = Objective(expr = sum(model.u[i] for i in model.N), 
                          sense = maximize)
    
    # Vincoli
    Ws = []       
    for i in model.N:
        Ws.append( Ks[i-1]*model.u[i]/100) # kWh coperti dal fotovoltaico
    # 1. massima disponiiblità totale della risorsa
    model.maxtot = ConstraintList()   
    model.maxtot.add( expr = sum(Ws) <= RC )
    # 2. massima copertura individuale della risorsa
    # [Infatti, se avanza energia del fotovoltaico, questa viene rimessa in circolo e venduta a Enel]
    model.maxind = ConstraintList()            
    for i in model.N:
        model.maxind.add( expr = Ws[i-1] <= Ks[i-1] )
    
    # Risoluzione con gurobi
    solver = SolverFactory('gurobi')
    sol = solver.solve(model, tee=False)

    # check feasibility
    sol_json = sol.json_repn()
    # Check solution status
    if sol_json['Solver'][0]['Status'] != 'ok':
        return None    
    
    perc_covered = [model.u[j]() for j in model.N]
        
    return perc_covered


# -----------------------------------------------
#   MAIN function
# -----------------------------------------------
if __name__ == "__main__":
    
    ### DATI
    ### Ps : lista delle presenze delle famiglie Ps[i]=0,1
    ### Ns : lista dei nomi delle famiglie
    ### Es : lista dell'energia dal fotovoltaico (fv)
    ### kWh : lista dei consumi-surplus di ogni famiglia
    ### Fs : lista dei consumi fissi (l'appartamento della famiglia i che è
    ###          in vacanza consuma Fs[i])
    ### costo : prezzo in €/kWh dell'energia
    ## Dal documento
    # Ps, Ns_old, Es_old, kWh_old, Fs_old  = ParseData('condominio2.txt')
    ## A mano
    Ns_old = ['Bianchi','Rossi','Verdi','Longo','Costa','Gatti']
    kWh_old = [1.452941176, 4.164705882, 1.970588235, 3.117647059, 3.529411765, 5.764705882]
    Ps = [1,1,1,0,0,1]
    Fs = [1.5, 1.5, 1.5, 1.5, 2, 2]
    Es_old = [2.647058824, 3.235294118, 3.529411765, 5.882352941, 6.470588235, 8.235294118]
    print('DATI')
    print('Presenze: {}'.format(Ps))
    print('Nomi: {}'.format(Ns_old))
    print('Dal fotovoltaico: {}'.format(Es_old))
    print('Surplus: {}'.format(kWh_old))
    print('Fisso: {} \n'.format(Fs))
    ##
    costo = 0.277 
    ##
    Ns_new, Es_new, kWh_new, RC = New(Ps, Ns_old, Es_old, kWh_old, Fs) 
    l = len(Ns_new)
    print('TAGLIO')
    print('Nomi: {}'.format(Ns_new))
    print('Dal fotovoltaico: {}'.format(Es_new))
    print('Surplus: {}'.format(kWh_new))
    print('Nuova energia: {} \n'.format(RC))
    
    ### SOLUZIONE
    perc_covered = MMF(kWh_new, RC) 
    
    ### PRINT
    kWh_covered = list(map(lambda x,y: np.multiply(x,y)/100, perc_covered, kWh_new))
    kWh_uncovered = list(map(lambda x,y: x-y , kWh_new, kWh_covered))
    costi = list(np.multiply(kWh_uncovered, costo))
    print('SOLUZIONE')
    # print('lung Ns: {}, lung kWh_covered: {}, lung costi: {}'.format(len(Ns_new), len(kWh_covered), len(costi)))
    for i in range(1,l+1):
        print('Appartamento {}, kWh coperti: {}, Costo: {}'.format(Ns_new[i-1], kWh_covered[i-1], costi[i-1]))
    
    ### POF
    FAIR = sum(perc_covered)
    print('FAIR: {}'.format(FAIR))