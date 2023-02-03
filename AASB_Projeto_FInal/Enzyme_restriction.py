class Enzyme_Restriction:
    def __init__(self, enzyme, seq):
        """Construtor da classe Enzyme_Restriction.

        Args:
            enzyme (str): Sequência da enzima
            seq (str): Sequência nucleótida
        """
        self.enzyme_seq = enzyme 
        self.seq = seq 

    def __str__ (self):
        """Função que permite representar as sequências como uma string.

        Returns:
            str: Representação das sequências
        """
        return "Enzyme sequence = " + str(self.enzyme_seq) + "\nNucleotide sequence = " + str(self.seq)

    def __repr__(self):
        """Função que permite representar as sequências como uma string na shell.

        Returns:
            str: Representação das sequências
        """
        return f"Enzyme_Restriction(\"{self.enzyme_seq}\",\"{self.seq}\")"

    def IUB(self):
        """Conversão da enzima de restrição na expressao regular correspondente.

        Returns:
            str: Expressão regular
        """
        try:
            rebase={"R":"[GA]", "Y":"[CT]","M":"[AC]", "K":"[GT]", "S":"[GC]", "W": "[AT]", "B":"[CTG]", "D":"[ATG]", "H":"[ACT]","V":"[ACG]", "N":"[ACGT]", "A":"A","C":"C","G":"G","T":"T"}
            cut = self.enzyme_seq.replace("^","")
            regexp = ""
            for x in cut:
                regexp += rebase[x]
            return regexp
        except:
            return "Nucleotide sequence or Enzyme invalid"

    def cut_positions(self): 
        """Função que permite determinar as posições onde a enzima de restrição "corta" a sequência.

        Returns:
            list: lista das posições de "corte"
        """
        from re import finditer
        cutposition = self.enzyme_seq.find("^")
        regexp = Enzyme_Restriction.IUB(self)
        matches = finditer(regexp, self.seq)
        locs = []
        for x in matches:
            locs.append(x.start() + cutposition)
        return locs

    def cut_positions_re(self): 
        """Função que permite determinar as posições onde a enzima de restrição "corta" a sequência complementar.

        Returns:
            list: lista das posições de "corte"
        """
        from re import finditer
        cutposition = self.enzyme_seq.find("^")
        seq_trans = self.seq.replace("A","t").replace("T","a").replace("C","g").replace("G","c")[::-1].upper()
        regexp = Enzyme_Restriction.IUB(self)
        matches = finditer(regexp, seq_trans)
        locs = []
        for x in matches:
            locs.append(x.start() + cutposition)
        return locs

    def cut_subsequences(self):
        """Função que determina as subsequências resultantes do corte da sequência (seq) usando os índices na lista locs.

        Returns:
            list: lista das subsequências
        """
        res =[]
        positions = Enzyme_Restriction.cut_positions(self)
        positions.insert(0,0)
        positions.append(len(self.seq))
        for i in range (len(positions)-1):
            res.append(self.seq[positions[i]:positions[i+1]])
        return res

    def cut_subsequences_re(self):
        """Função que determina as subsequências resultantes do corte da sequência completamentar (seq) usando os índices na lista locs.

        Returns:
            list: lista das subsequências
        """
        res =[]
        positions = Enzyme_Restriction(self.enzyme_seq,self.seq).cut_positions()
        positions.insert(0,0)
        seq_trans = self.seq.replace("A","t").replace("T","a").replace("C","g").replace("G","c")[::-1].upper()
        positions.append(len(self.seq))
        for i in range (len(positions)-1):
            res.append(seq_trans[positions[i]:positions[i+1]])
        return res

import unittest
class Test_Enzyme_Restriction(unittest.TestCase):

    def test_IUB(self):
        """Função usada para testar a função IUB() da classe Enzyme_Restriction.
        """
        self.assertEqual(Enzyme_Restriction("G^AMTV","TAAGACTGAAT").IUB(),"GA[AC]T[ACG]")
        self.assertEqual(Enzyme_Restriction("GPAMTV","TAAGACTGAAT").IUB(),"Nucleotide sequence or Enzyme invalid")
        self.assertEqual(Enzyme_Restriction("GPAMTV","TGCAACGU").IUB(),"Nucleotide sequence or Enzyme invalid")
        self.assertEqual(Enzyme_Restriction("G_AMTV","TAAGACTGAAT").IUB(),"Nucleotide sequence or Enzyme invalid")
    
    def test_cut_positions(self):
        """Função usada para testar a função cut_positions() da classe Enzyme_Restriction.
        """
        self.assertEqual(Enzyme_Restriction("G^AMTV","TAAGACTGAAT").cut_positions(), [4])
        self.assertEqual(Enzyme_Restriction("GAAACT","TAAGACTGAAT").cut_positions(), [])
        self.assertEqual(Enzyme_Restriction("GA    AACT","TAAGACTGAAT").cut_positions(), [])
        self.assertEqual(Enzyme_Restriction("GA*_*PAACT","TAAGACTGAAT").cut_positions(), [])

    def test_cut_positions_re(self):
        """Função usada para testar a função cut_positions_re() da classe Enzyme_Restriction.
        """
        self.assertEqual(Enzyme_Restriction("G^AMTV","ATTCAGTCTTA").cut_positions_re(), [4])
        self.assertEqual(Enzyme_Restriction("G^AMTV","ATAGTCTCAGTCTTA").cut_positions_re(), [4, 10])
        self.assertEqual(Enzyme_Restriction("GA    AACT","ATAGTCTCAGTCTTA").cut_positions_re(), [])
        self.assertEqual(Enzyme_Restriction("GA12AACT","TAAGACTGAAT").cut_positions_re(), [])

    def test_cut_subsequences(self):
        """Função usada para testar a função cut_subsequences() da classe Enzyme_Restriction.
        """
        self.assertEqual(Enzyme_Restriction("G^AMTV","TAAGACTGAAT").cut_subsequences(), ['TAAG', 'ACTGAAT'])
        self.assertEqual(Enzyme_Restriction('GA*_*PAACT', 'TAAGACTGAAT').cut_subsequences(), ["TAAGACTGAAT"])
        self.assertEqual(Enzyme_Restriction('GA   PAACT', 'TAAGACTGAAT').cut_subsequences(), ["TAAGACTGAAT"])
        self.assertEqual(Enzyme_Restriction('PPPAA12', 'TAAGACTGAAT').cut_subsequences(), ["TAAGACTGAAT"])

    def test_cut_subsequences_re(self):
        """Função usada para testar a função cut_subsequences_re() da classe Enzyme_Restriction.
        """
        self.assertEqual(Enzyme_Restriction("G^AMTV","GTAGAAGATTCT").cut_subsequences_re(), ['AGAATCTTCTAC'])
        self.assertEqual(Enzyme_Restriction("G^AMTV","TAAGACTGAAT").cut_subsequences_re(), ['ATTC', 'AGTCTTA'])
        self.assertEqual(Enzyme_Restriction("GA   PAACT","ATTCAGTCTTA").cut_subsequences_re(), ['TAAGACTGAAT'])
        self.assertEqual(Enzyme_Restriction("GA1223PAACT","ATTCAGTCTTA").cut_subsequences_re(), ['TAAGACTGAAT'])


unittest.main(argv=[""], exit=False)