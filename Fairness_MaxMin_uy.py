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
    model.y = Var(model.N, model.N, domain = Binary)
    model.u = Var(model.N, domain = NonNegativeReals, bounds = (0,100) )

    # Funzione obiettivo
    ### Devo massimizzare lessicograficamente
    model.obj = Objective(expr = sum((2*n-i)*sum(model.y[i,j]*model.u[j] for j in model.N)  for i in model.N), 
                          sense = maximize)
    ##
    # model.funob = ObjectiveList()
    # for i in model.N:
    #     model.funob.add( (sum(model.y[i,j]*model.u[i] for j in model.N)), i-1, n-i)
    # model.ModelSense = GRB.MAXIMIZE
    ##
    # IDEA: massimizzare la percentuale della prima riga, poi quella della seconda e così via fino all'ultima;
    # ogni step (ogni riga i di y), rispetto allo step precendete, ha almeno una famiglia in più (un valore di y[i,j]) 
    # che ha una percentuale nulla (cioè non viene aumentata la sua porzione che verrà coperta).
    # model.setObjectiveN((sum(model.y[1,j]*model.u[1] for j in model.N)), 0, 4)
    # model.setObjectiveN((sum(model.y[2,j]*model.u[2] for j in model.N)), 1, 3)
    # model.setObjectiveN((sum(model.y[3,j]*model.u[3] for j in model.N)), 2, 2)
    # model.setObjectiveN((sum(model.y[4,j]*model.u[4] for j in model.N)), 3, 1)
    
    # Vincoli
    Ws = []       
    for j in model.N: # kWh coperti dal fotovoltaico per appartamento
        Ws.append( Ks[j-1]*sum(model.y[i,j]*model.u[i] for i in model.N)/100) 
    # 1. massima disponiblità totale della risorsa
    model.maxtot = ConstraintList()   
    model.maxtot.add( expr = sum(Ws) <= RC )
    # 2. massima copertura individuale della risorsa
    # [Infatti, se avanza energia del fotovoltaico, questa viene rimessa in circolo e venduta a Enel]
    model.maxind = ConstraintList()            
    for i in model.N:
        model.maxind.add( expr = Ws[i-1] <= Ks[i-1] )
    # 3. ordine
    model.zeros = ConstraintList()
    for i in model.N: # per ogni riga
        for j in model.N: 
            # se c'è un elemento della riga nullo allora tutte le righe 
            # successive di quella colonna sono nulle
            model.zeros.add( expr = sum(model.y[k,j] for k in range(i,n+1)) <= n*model.y[i,j] )
    # 4. minima percentuale
    model.ones = ConstraintList()
    model.ones.add( expr = sum(model.y[1,j] for j in model.N) == n )
    for i in range(2,n+1):
        model.ones.add( expr = sum(model.y[i,j] for j in model.N) <= sum(model.y[i-1,j] for j in model.N) )
    
    
    
    # Risoluzione con gurobi
    solver = SolverFactory('gurobi')
    sol = solver.solve(model, tee=False)

    # check feasibility
    sol_json = sol.json_repn()
    # Check solution status
    if sol_json['Solver'][0]['Status'] != 'ok':
        return None    
    
    print('DEBUGGING')
    U = [model.u[j]() for j in model.N] # matrice delle percentuali
    Y = [[model.y[i,j]() for j in model.N] for i in model.N] # matrice binaria
    print('Matrice binaria: {}'.format(Y))
    print('Vettore perc: {}'.format(U))
    A = []
    for j in range(0,n):
        A.append( list( np.multiply([Y[j][i] for i in range(0,n)], U[j]) ) )
    print('Matrice perc: {}'.format(A))
    perc_covered = [] # vettore delle quantità di energia corrispondenti alle percentuali 
    for j in range(0,n):
        perc_covered.append( sum(Y[i][j]*U[i] for i in range(0,n)) )
    print('Energia: {}'.format(perc_covered))
        
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