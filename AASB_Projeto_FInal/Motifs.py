class Motifs:
    """Classe que permite gerar Motifs probabilisticos 
    """
    def __init__(self, list_seqs):
        """Contrutor do Motifs onde criámos o nosso Alfabeto(letters) e onde gerámos a nossa Matriz vazia (mat)

        Args:
            list_seqs (list): Lista contendo diferentes sequências biológicas
        """
        from collections import Counter
        self.seqs = list_seqs
        self.letters = list(Counter([a for c in list_seqs for a in c]))
        self.mat = Mat(len(self.letters)+1,len(list_seqs[0])+1)
    
    def __repr__(self):
        return print(self.mat)

    def __str__(self) -> str:
        return print(self.mat)

    def contagens (self):
        """Função que permite devolver a frequência dos elementos das sequências

        Returns:
            list: Matriz das Contagens 
        """
        for c in self.seqs:
            for coluna in range(len(self.seqs[0])):
                linha = self.letters.index(c[coluna])
                self.mat[linha + 1][coluna + 1] += 1
                self.mat[0][coluna + 1] = coluna +1 
                self.mat[0][0] = "."                
                self.mat[linha +1][0] = self.letters[linha] 
        return print(self.mat)


    def PWM (self, pseudo_count = 0):
        """Função que permite gerar a Matriz PWM

        Args:
            pseudo_count (int, optional): pseudocontagem. Defaults to 0.

        Returns:
            list: Matriz PWM
        """
        self.mat = Mat(len(self.letters)+1,len(self.seqs[0])+1)
        for c in self.seqs:
            for coluna in range(len(self.seqs[0])):
                linha = self.letters.index(c[coluna])
                self.mat[linha + 1][coluna + 1] += 1
                self.mat[0][coluna + 1] = f"  {coluna +1} "   " "  
                self.mat[0][0] = "."                
                self.mat[linha +1][0] = self.letters[linha] 
        try:      
            for c in range(1,len(self.letters)+ 1):
                for a in range(1,len(self.seqs)+ 1): 
                    self.mat[c][a] = round(float(self.mat[c][a] + pseudo_count + 0.01) / float(len(self.seqs)),3)
            return print(self.mat)

        except:
            for c in range(1,len(self.letters)+ 1):
                for a in range(1,len(self.seqs)+ 2):
                    self.mat[c][a] = round(float(self.mat[c][a] + pseudo_count + 0.01) / float(len(self.seqs)),3)
            return print(self.mat)

    
    def PSSM (self, pseudo_count = 0):
        """Função que permite gerar a Matriz PSSM

        Args:
            pseudo_count (int, optional): pseudocontagem. Defaults to 0.

        Returns:
            list: Matriz PSSM
        """
        self.mat = Mat(len(self.letters)+1,len(self.seqs[0])+1)
        from math import log2
        for c in self.seqs:
            for coluna in range(len(self.seqs[0])):
                linha = self.letters.index(c[coluna])
                self.mat[linha + 1][coluna + 1] += 1
                self.mat[0][coluna + 1] = f"  {coluna +1} "   " "  
                self.mat[0][0] = "."               
                self.mat[linha +1][0] = self.letters[linha] 
        try:
            for c in range(1,len(self.letters)+ 1):
                for a in range(1,len(self.seqs)+ 1): 
                    self.mat[c][a] = round((log2(((self.mat[c][a]+ pseudo_count)/(len(self.seqs)+ len(self.letters)* pseudo_count)/(1/len(self.letters)) ))),3)
            return print(self.mat)
        except:
            for c in range(1,len(self.letters)+ 1):
                for a in range(1,len(self.seqs)+ 2):
                    self.mat[c][a] = round((log2(((self.mat[c][a]+ pseudo_count)/(len(self.seqs)+ len(self.letters* pseudo_count)/(1/len(self.letters)) )))),3)
            return print(self.mat)
    
    def prob_seq(self, seq):
        """Função que permite devolver a probabilidade de uma sequência

        Args:
            seq (str): Sequência 

        Returns:
            float: Probabilidade da sequência 
        """

        self.mat = Mat(len(self.letters)+1,len(self.seqs[0])+1)
        prob = 1 
        for c in self.seqs:
            for coluna in range(len(self.seqs[0])):
                linha = self.letters.index(c[coluna])
                self.mat[linha + 1][coluna + 1] += 1
                self.mat[0][coluna + 1] = f"  {coluna +1} "   " "  
                self.mat[0][0] = "."               
                self.mat[linha +1][0] = self.letters[linha] 
        for c in range(1,len(self.letters)+ 1):
            for a in range(1,len(self.seqs)+ 1):
                self.mat[c][a] = round(float(self.mat[c][a] + 0.01) / float(len(self.seqs)),3)
        
        for c in range(1,len(self.seqs[0])+ 1):
            a = self.letters.index(seq[c - 1])
            prob *= self.mat[a +1][c]
        return round(prob,4)


    def most_prob_seq (self,pseudo_count = 0):
        """Função que permite devolver a sequência mais provável de ocorrer

        Args:
            pseudo_count (int, optional): pseudocontagem. Defaults to 0.

        Returns:
            str: Sequência mais provável de ocorrer
        """

        self.mat = Mat(len(self.letters)+1,len(self.seqs[0])+1)
        resul = ""
        for c in self.seqs:
            for coluna in range(len(self.seqs[0])):
                linha = self.letters.index(c[coluna])
                self.mat[linha + 1][coluna + 1] += 1
                self.mat[0][coluna + 1] = f"  {coluna +1} "   " "  
                self.mat[0][0] = "."               
                self.mat[linha +1][0] = self.letters[linha]
        for c in range(1,len(self.letters)+ 1):
            for a in range(1,len(self.seqs)+ 1):
                self.mat[c][a] = round(float(self.mat[c][a] + pseudo_count + 0.01) / float(len(self.seqs)),3)
        list_temp= []
        list_final = []
        for c in range(1,len(self.seqs)+1):
            for a in range(1,len(self.seqs)+ 1):
                list_temp.append((self.mat[a][c],a))
            list_final.append(max(list_temp))
            list_temp =[]
        
        for c in list_final:
            for a in range(1,len(self.seqs[0])+1):
                if a == c[1]:
                    resul += self.mat[a][0]
        return resul

import unittest  
seqs = ['ATTG','ATCG','ATTC','ACTC']

class Test_Motifs (unittest.TestCase):
    """Conjunto de testes para a Classe Motifs
    """

    def test_prob_seq (self):
        """Testar a função prob_seq
        """

        self.assertEqual(Motifs(seqs).prob_seq("ATTG"), 0.2845)


    def test_most_prob_seq (self):
        """Testar a função most_prob_seq
        """
        
        self.assertEqual(Motifs(seqs).most_prob_seq(), "ATTC")

unittest.main(argv=[''], exit=False)

    