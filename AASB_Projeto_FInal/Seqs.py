class Seq:
    """ 
    Classe usada para manipular e transformar os diferentes tipos de sequências através de diferentes processos.
    
    """
    def __init__(self, seq):
        """Construtor da classe Seq.

        Args:
            seq (str): Sequência biológica.
        """
        self.seq = str(seq.upper())

    def __len__ (self): 
        """ Função que permite implementar a função len().

        Returns:
            int: Tamanho da sequência.
        """   
        return len(self.seq)
    
    def __str__(self):
        """ Função que permite representar a sequência como uma string.

        Returns:
            str: Representação da sequência.      
        """       
        return self.seq

    def __getitem__(self,n):
        """Função que permite aplicar a indexação na classe.

        Args:
            n (int): Index

        Returns:
            str: Base ou aminoácido correspondente ao index.
        """
        return self.seq[n]
    
    def __getslice__(self,i,j):
        """Função que permite aplicar a indexação na classe.

        Args:
            i (int): index start
            j (int): index final

        Returns:
            str: Sequência, base ou aminoácido correspondente ao index.
        """
        return self.seq[a]    

    def get_type (self):
        """Função que a partir de uma sequência, classifica-a em DNA, RNA ou AMINO ACID SEQUENCE.

        Returns:
            str: Classificação da sequência.
        """
        if len([c for c in self.seq if c not in "ACGT"]) == 0:
            return "DNA"

        elif len([c for c in self.seq if c not in "ACGU"]) == 0:
            return "RNA"

        elif len([c for c in self.seq if c not in "ABCDEFGHIKLMNPQRSTVWYZ_"]) == 0:
            return "AMINO ACID SEQUENCE"
        
        else:
            return "Not valid"

    def count_nucleotides (self):
        """Função que a partir de uma sequência de DNA ou RNA e conta o número dos diferentes nucleótidos presentes.

        Returns: 
            list: Lista com a contagem dos diferentes núcleotidos.           
        """
        if Seq(self.seq).get_type() == "DNA":
            A = self.seq.count("A")
            C = self.seq.count("C")
            G = self.seq.count("G")
            T = self.seq.count("T") 
            return [A, C, G, T]
        
        elif Seq(self.seq).get_type() == "RNA":                    
            A = self.seq.count("A")
            C = self.seq.count("C")
            G = self.seq.count("G")
            U = self.seq.count("U") 
            return [A, C, G, U] 
        else:
            return "Not DNA or RNA"

    def freq_nucleotides (self):
        """Função que a partir de uma sequência de DNA ou RNA, estima a frequência dos diferentes nucleótidos presentes.

        Returns:
            list: Lista com a frequência dos diferentes nucleótidos. 
        """
        if Seq(self.seq).get_type() == "DNA":
            A = self.seq.count("A") / len(self.seq)
            C = self.seq.count("C") / len(self.seq)
            G = self.seq.count("G") / len(self.seq)
            T = self.seq.count("T") / len(self.seq)
            return [round(A, 3), round(C, 3), round(G, 3), round(T, 3)]             
        
        elif Seq(self.seq).get_type() == "RNA":
            A = self.seq.count("A") / len(self.seq)
            C = self.seq.count("C") / len(self.seq)
            G = self.seq.count("G") / len(self.seq)
            U = self.seq.count("U") / len(self.seq)
            return [round(A, 3), round(C, 3), round(G, 3), round(U, 3)]
        else:
            return "Not DNA or RNA"

    def transcricao(self):
        """Função que a partir de uma sequência de DNA ou RNA, devolve o transcrito da mesma.

        Returns:
            str: Sequência do transcrito.
        """
        if Seq(self.seq).get_type() == "DNA":  
           DNA = self.seq.lower() 
           seq_trans = DNA.replace("a","U").replace("u","A").replace("c","G").replace("g","C")[::-1].upper()
           return seq_trans  

        elif Seq(self.seq).get_type() == "RNA": 
           RNA = self.seq.lower() 
           seq_trans = RNA.replace("u","A").replace("a","T").replace("g","C").replace("c","G")[::-1].upper()
           return seq_trans  

        else:
            return "Not DNA or RNA" 

    def get_codons (self):
        """Função que recebe uma sequência de DNA ou RNA e transforma cada conjunto de 3 nucleótidos, presentes ao longo da sequência, no respetivo codão.  

        Returns:
            list: Lista de codões.
        """
        if Seq(self.seq).get_type() == "DNA" or Seq(self.seq).get_type() == "RNA":
            l = []
            codon_list =[]
            for c in range (0, len(self.seq), 3):                             
                codon_list.append(self.seq[c : c + 3])
            if len(codon_list[-1]) < 3:
                codon_list.pop()
            return codon_list
        else:
            return "Not DNA or RNA"   

    def traducao(self):
        """"Função que a partir de uma sequência de DNA, origina uma sequência de aminoácidos. 

        Returns:
            str: Sequência de aminoácidos
        """
        if Seq(self.seq).get_type() == "DNA":                                 
            gencode = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
        
            l = []
            for c in Seq(self.seq).get_codons():
                l.append(gencode.get(c))
                l_final = "".join(l)    
            return l_final  
        else:
            return "Not DNA"

    def get_prots (self):   
        """Função que a partir de uma sequência de aminoácidos gera as diferentes proteínas.

        Returns:
            list: Lista de proteínas
        """
        if Seq(self.seq).get_type() == "AMINO ACID SEQUENCE":
            lista_final = []
            lista_stop_tratados = self.seq.split("_") 
            for c in lista_stop_tratados:        
                for idx, amino_s in enumerate(c):       
                    if amino_s == "M":                  
                        lista_final.append(c[idx:]) 
            return lista_final 
        else:
            return "not AMINO ACID SEQUENCE"

    def get_orfs (self):    
        """Função que a partir de uma sequência de DNA ou RNA origina diferentes tipos de ORFS. 

        Returns:
            list: Lista das ORFS
        """
        if  Seq(self.seq).get_type() == "DNA": 
            dna = self.seq.lower()
            dna_f = self.seq
            lista_orfs =[]
            lista_orfs.append(dna_f[0:])   
            lista_orfs.append(dna_f[1:])   
            lista_orfs.append(dna_f[2:])
            dna_rev_comp = dna.replace("a","T").replace("t","A").replace("c","G").replace("g","C")[::-1] 
            lista_orfs.append(dna_rev_comp[0:])
            lista_orfs.append(dna_rev_comp[1:])  
            lista_orfs.append(dna_rev_comp[2:])
            return lista_orfs
        elif Seq(self.seq).get_type() == "RNA":
            rna = self.seq.lower()
            rna_f = self.seq
            lista_orfs =[]
            lista_orfs.append(rna_f[0:])      
            lista_orfs.append(rna_f[1:])   
            lista_orfs.append(rna_f[2:])
            rna_rev_comp = rna.replace("a","U").replace("u","A").replace("c","G").replace("g","C")[::-1] 
            lista_orfs.append(rna_rev_comp[0:])
            lista_orfs.append(rna_rev_comp[1:]) 
            lista_orfs.append(rna_rev_comp[2:])
            return lista_orfs
        else:
            return "Not DNA or RNA"

    def get_all_prots(self):
        """Função que a partir de uma sequência de DNA e permite obter todas as possíveis proteínas. 

        Returns:
            list: Lista das proteínas
        """
        if  Seq(self.seq).get_type() == "DNA": 
            orfs = Seq.get_orfs(self)
            list_codons = []
            list_aminos = []
            list_prots =  []
            
            for c  in orfs:
                a = Seq(c).get_codons()
                list_codons.append(a)
                for e in list_codons:
                    b = Seq(c).traducao()
                list_aminos.append((b))
                for r in list_aminos:
                    c = Seq(r).get_prots()
                list_prots.append(c)       
            return [c for a in list_prots for c in a]
        else:
            return "Not DNA"


import unittest
class Test_Seq(unittest.TestCase):
    
    def test_get_type(self):
        """Função usada para testar a função get_type() da classe Seq.
        """
        self.assertEqual(Seq("ATGC").get_type(),"DNA")
        self.assertEqual(Seq("augc").get_type(),"RNA")
        self.assertEqual(Seq("MNAAQTHCTA").get_type(),"AMINO ACID SEQUENCE")
        self.assertEqual(Seq("A-G-T-A-T").get_type(),"Not valid")
        self.assertEqual(Seq("ATGC1").get_type(),"Not valid")

    def test_count_nucleotides(self):
        """Função usada para testar a função count_nucleotides() da classe Seq.
        """
        self.assertEqual(Seq("ATTTGC").count_nucleotides(),[1, 1, 1, 3])
        self.assertEqual(Seq("auuugc").count_nucleotides(),[1, 1, 1, 3]) 
        self.assertEqual(Seq("AUTUUGC").count_nucleotides(),"Not DNA or RNA") 

    def test_freq_nucleotides(self):
        """Função usada para testar a função freq_nucleotides() da classe Seq.
        """
        self.assertEqual(Seq("ATTTGC").freq_nucleotides(),[0.167, 0.167, 0.167, 0.5])
        self.assertEqual(Seq("auuagc").freq_nucleotides(),[0.333, 0.167, 0.167, 0.333])
        self.assertEqual(Seq("AT TAC-AC").freq_nucleotides(),"Not DNA or RNA")
        
    def test_transcricao(self):
        """Função usada para testar a função transcricao() da classe Seq.
        """
        self.assertEqual(Seq("ATTTGC").transcricao(),"GCTTTU")
        self.assertEqual(Seq("auuagc").transcricao(),"GCTAAT")
        self.assertEqual(Seq("AUUHAGC").transcricao(),"Not DNA or RNA")
    
    def test_get_codons(self):
        """Função usada para testar a função get_codons() da classe Seq.
        """
        self.assertEqual(Seq("ATTTGC").get_codons(),['ATT', 'TGC'])
        self.assertEqual(Seq("auuagc").get_codons(),['AUU', 'AGC'])
        self.assertEqual(Seq("acacc_ca").get_codons(),"Not DNA or RNA")
        self.assertEqual(Seq("ATTTGC1").get_codons(),"Not DNA or RNA")
    
    def test_traducao(self):
        """Função usada para testar a função traducao() da classe Seq.
        """
        self.assertEqual(Seq("ATTTGC").traducao(),"IC")
        self.assertEqual(Seq("auuagc").traducao(),"Not DNA")

    def test_get_prots(self):
        """Função usada para testar a função get_prots() da classe Seq.
        """
        self.assertEqual(Seq("MA_KLCMNB").get_prots(),['MA', 'MNB'])
        self.assertEqual(Seq("ATTMTGC").get_prots(),['MTGC'])
        self.assertEqual(Seq("ATTTGC").get_prots(),"not AMINO ACID SEQUENCE")
        self.assertEqual(Seq("auuagc").get_prots(),"not AMINO ACID SEQUENCE")
        
    def test_get_orfs(self):
        """Função usada para testar a função get_orfs() da classe Seq.
        """
        self.assertEqual(Seq("ATTTGC").get_orfs(),['ATTTGC', 'TTTGC', 'TTGC', 'GCAAAT', 'CAAAT', 'AAAT'])
        self.assertEqual(Seq("auuagc").get_orfs(),['AUUAGC', 'UUAGC', 'UAGC', 'GCUAAU', 'CUAAU', 'UAAU'])
        self.assertEqual(Seq("MA_KLCMNB").get_orfs(),"Not DNA or RNA")
        self.assertEqual(Seq("ATTTGC1").get_orfs(),"Not DNA or RNA")

    def test_get_all_prots(self):
        """Função usada para testar a função get_all_prots() da classe Seq.
        """
        self.assertEqual(Seq("ATTTGC").get_all_prots(),[])
        self.assertEqual(Seq("ATTATGTGC").get_all_prots(),['MC'])
        self.assertEqual(Seq("auguagc").get_all_prots(),"Not DNA")
    
unittest.main(argv=[""], exit=False)