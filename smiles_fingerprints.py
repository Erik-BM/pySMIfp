


def islower(string):
    return string.lower() == string

def isupper(string):
    return string.upper() == string

def isnumber(string):
    return '1' <= string and string <= '9'

def smiles_fingerprint(smiles):
    BrkRndOp = BrkSqrOp = 0 #Brackets
    SqrBrk = False #Open Brackets
    SB=DB=TB = 0 #Bond types

    ChargePlus = ChargeMinus = 0 # Charges
    AtomB = AtomC = AtomN = AtomO = AtomP = AtomS = AtomF = AtomCl = AtomBr = AtomI = AtomH = 0
    Atomc =  Atomo = Atoms = Atomn = Atomp = 0
    AtomOther = 0
    RingIdx = [0 for _ in range(11)]
    PerCent = 0
    
    for i in range(len(smiles)):
        
        currChar = smiles[i]
        nextChar = ' '
        if i != len(smiles) - 1:
            nextChar = smiles[i+1]
        
        if not SqrBrk:
            if currChar == 'C':
                if nextChar == 'l':
                    AtomCl += 1
                else:
                    AtomC += 1
            elif currChar == 'B':
                if nextChar == 'r':
                    AtomBr += 1
                else:
                    AtomB += 1
            elif currChar == '[':
                BrkSqrOp += 1
                SqrBrk = True
            else:
                for count,v in zip([AtomN,AtomO,AtomP,AtomS,AtomF,AtomI,Atomc,Atomo,Atoms,Atomn,Atomp,BrkRndOp,DB,TB,PerCent,SB],
                                   ['N','O','P','S','F','I','c','o','s','n','p','(','=','#','-']):
                    if currChar == v: count += 1
                for j,v in enumerate(['1','2','3','4','5','6','7','8','9','0']):
                    if currChar == v: RingIdx[j] += 1
        elif SqrBrk:
            if currChar == ']':
                SqrBrk = False
            elif currChar == '-':
                if isnumber(nextChar):
                    ChargeMinus += int(nextChar.replace('0',''))
                else:
                    ChargeMinus += 1
            elif currChar == '+':
                if isnumber(nextChar):
                    ChargePlus += int(nextChar.replace('0',''))
                else:
                    ChargePlus += 1
                    
            elif islower(currChar):
                for count,v in zip([Atomo,Atomc,Atomn,Atomp,Atoms],
                                   ['o','c','n','p','s']):
                    if currChar == v: count += 1
                        
            elif isupper(currChar):
                if currChar == 'C':
                    if islower(nextChar):
                        if nextChar == 'l':
                            AtomCl += 1
                            i += 1
                        else:
                            AtomOther +=1
                            i += 1
                    else:
                        AtomC += 1
                elif currChar == 'B':
                    if islower(nextChar):
                        if nextChar == 'r':
                            AtomBr += 1
                            i +=1
                        else:
                            AtomOther +=1
                            i += 1
                    else:
                        AtomB += 1
                elif currChar == 'H':
                    if islower(nextChar):
                        AtomOther += 1
                        i += 1
                    elif isnumber(nextChar):
                        AtomH += int(nextChar.replace('0',''))
                    else:
                        AtomH += 1
                
                elif currChar in ['N','O','P','S','F','I']:
                    for count,v in zip([AtomN,AtomO,AtomP,AtomS,AtomF,AtomI],
                                    ['N','O','P','S','F','I']):
                        if currChar == v and islower(currChar): 
                            AtomOther += 1
                            i += 1
                        else:
                            count += 1
                else:
                    AtomOther += 1
                    if islower(nextChar):
                        i += 1
    
    SMIfp = [0 for _ in range(34)]
    for ri,smi in zip([0,1,2,3,4,5,6,7,8],[8,7,10,12,18,22,23,31,32]):
        SMIfp[smi] = RingIdx[ri]/2
    
    for i,v in zip([0,13,19,2,9,26,21,20,28,5,3,1,27,11,24,17,25,30,4,14,16,6,33,15,29],
                   [BrkRndOp,BrkSqrOp,SB,DB,TB,PerCent,ChargePlus,ChargeMinus,AtomB,AtomC,AtomN,AtomO,AtomP,AtomS,AtomF,AtomCl,AtomBr,AtomI,Atomc,Atomo,Atoms,Atomn,Atomp,AtomH,AtomOther]):
        SMIfp[i] = v
        
    return SMIfp
                








 
