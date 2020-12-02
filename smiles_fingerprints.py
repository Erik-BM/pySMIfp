


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
            elif currChar == 'N':
                AtomN += 1
            elif currChar == 'O':
                AtomO += 1
            elif currChar == 'P':
                AtomP += 1
            elif currChar == 'S':
                AtomS += 1
            elif currChar == 'F':
                AtomF += 1
            elif currChar == 'I':
                AtomI += 1
            elif currChar == 'c':
                Atomc += 1
            elif currChar == 'o':
                Atomo += 1
            elif currChar == 's':
                Atoms += 1
            elif currChar == 'n':
                Atomn += 1
            elif currChar == 'p':
                Atomp += 1
            elif currChar == '(':
                BrkRndOp += 1
            elif currChar == '[':
                BrkSqrOp += 1
                SqrBrk = True
            elif currChar == '=':
                DB += 1
            elif currChar == '#':
                TB += 1
            elif currChar == '1':
                RingIdx[0] += 1
            elif currChar == '2':
                RingIdx[1] += 1
            elif currChar == '3':
                RingIdx[2] += 1
            elif currChar == '4':
                RingIdx[3] += 1
            elif currChar == '5':
                RingIdx[4] += 1
            elif currChar == '6':
                RingIdx[5] += 1
            elif currChar == '7':
                RingIdx[6] += 1
            elif currChar == '8':
                RingIdx[7] += 1
            elif currChar == '9':
                RingIdx[8] += 1
            elif currChar == '0':
                RingIdx[9] += 1
            elif currChar == '%':
                PerCent += 1
            elif currChar == '-':
                SB += 1
        elif SqrBrk:
            if currChar == ']':
                SqrBrk = False
            elif currChar == '-':
                if '1' <= nextChar and nextChar <= '9':
                    ChargeMinus += nextChar - '0'
                else:
                    ChargeMinus += 1
            elif currChar == '+':
                if '1' <= nextChar and nextChar <= '9':
                    ChargePlus += nextChar - '0'
                else:
                    ChargePlus += 1
            elif currChar.lower() == currChar:
                if currChar == 'o':
                    Atomo += 1
                elif currChar == 'c':
                    Atomc += 1
                elif currChar == 'n':
                    Atomn += 1
                elif currChar == 'p':
                    Atomp += 1
                elif currChar == 's':
                    Atoms += 1
            elif currChar.upper() == currChar:
                if currChar == 'C':
                    if nextChar.lower() == nextChar:
                        if nextChar == 'l':
                            AtomCl += 1
                            i += 1
                        else:
                            AtomOther +=1
                            i += 1
                    else:
                        AtomC += 1
                elif currChar == 'B':
                    if nextChar.lower() == nextChar:
                        if nextChar == 'r':
                            AtomBr += 1
                            i +=1
                        else:
                            AtomOther +=1
                            i += 1
                    else:
                        AtomB += 1
                elif currChar == 'N':
                    if nextChar.lower == nextChar:
                        AtomOther += 1
                        i += 1
                    else:
                        AtomN += 1
                elif currChar == 'O':
                    if nextChar.lower == nextChar:
                        AtomOther += 1
                        i += 1
                    else:
                        AtomO += 1
                elif currChar == 'P':
                    if nextChar.lower == nextChar:
                        AtomOther += 1
                        i += 1
                    else:
                        AtomP += 1
                elif currChar == 'S':
                    if nextChar.lower == nextChar:
                        AtomOther += 1
                        i += 1
                    else:
                        AtomS += 1
                elif currChar == 'F':
                    if nextChar.lower == nextChar:
                        AtomOther += 1
                        i += 1
                    else:
                        AtomF += 1
                elif currChar == 'I':
                    if nextChar.lower == nextChar:
                        AtomOther += 1
                        i += 1
                    else:
                        AtomI += 1
                elif currChar == 'H':
                    if nextChar.lower() == nextChar:
                        AtomOther += 1
                        i += 1
                    elif '1'<= nextChar and nextChar <= '9':
                        AtomH += nextChar - '0'
                    else:
                        AtomH += 1
                else:
                    AtomOther += 1
                    if nextChar.lower() == nextChar:
                        i += 1
    
    SMIfp = list(range(34))
    SMIfp[0] = BrkRndOp
    SMIfp[13] = BrkSqrOp
    SMIfp[19] = SB
    SMIfp[2] = DB
    SMIfp[9] = TB
    
    for ri,smi in zip([0,1,2,3,4,5,6,7,8],[8,7,10,12,18,22,23,31,32]):
        SMIfp[smi] = RingIdx[ri]
        if RingIdx[ri] != 0:
            SMIfp[smi] = RingIdx[ri]/2
            
    SMIfp[26] = PerCent
    SMIfp[21] = ChargePlus
    SMIfp[20] = ChargeMinus
    SMIfp[28] = AtomB
    SMIfp[5] = AtomC
    SMIfp[3] = AtomN
    SMIfp[1] = AtomO
    SMIfp[27] = AtomP
    SMIfp[11] = AtomS
    SMIfp[24] = AtomF
    SMIfp[17] = AtomCl
    SMIfp[25] = AtomBr
    SMIfp[30] = AtomI
    SMIfp[4] = Atomc
    SMIfp[14] = Atomo
    SMIfp[16] = Atoms
    SMIfp[6] = Atomn
    SMIfp[33] = Atomp
    SMIfp[15] = AtomH
    SMIfp[29] = AtomOther
    return SMIfp
                








 
