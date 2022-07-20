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
    for i in range(1,n+1):
        if Ps[i-1] == 1: # se la fam è presente considero tutto
            Ns_new.append(Ns[i-1])
            Es_new.append(Es[i-1])
            kWh_new.append(kWh[i-1])
        else: # se la fam è assente
            diff = Es[i] - Fs[i] # energia scoperta dal FV
            if (diff) <= 0: # se il fotovoltaico non copre il consumo fisso, 
                            # allora kWh != 0
                Ns_new.append(Ns[i-1])
                Es_new.append(Es[i-1])
                kWh_new.append(-diff)
            else: # se il fotovoltaico copre il consumo fisso e avanza anche energia,
                  # allora quell'energia in più va in rc (a disposizione degli altri)
                  rc += (Es[i] - Fs[i])
    return Ns_new, Es_new, kWh_new, rc