class Mat:
    """
    Classe que permite criar e manipular Matrizes
   
    """
    def __init__(self, rows, cols):
        """Construtor da Matriz 

        Args:
            rows (int): Número de linhas
            cols (int): Número de colunas 
        """
        self.mat = [[0 for c in range(cols)]
                    for r in range(rows)]
    

    def numRows (self): 
        """Função que permite obter o número de linhas da Matriz

        Returns:
            int : Valor de linhas da Matriz
        """
        return len(self.mat)


    def numCols (self): 
        """Função que permite obter o númeor de colunas da Matriz

        Returns:
            int : Valor de colunas da Matriz
        """
        return len(self.mat[0])


    def __str__(self):
        """ Função que permite devolver a Matriz como uma string

        Returns:
            str : String da Matriz
        """
        return '\n'.join(' '.join(str(val) for val in row)
                         for row in self.mat)

    
    def __getitem__ (self, n):
        return self.mat[n]

import io

def subst(x, y):
  import blosum as bl

  """Funçao que permite obter o valor da substituição do Aminoácido da Blossum 62

  Args:
      x (str): Aminoácido 1
      y (str): Aminoácido 2

  Returns:
      int: Valor de substituição do Aminoácido 
  """
  
  dic = dict(bl.BLOSUM(62))
  return int(dic[x+y])   

class NW:
  """Classe que permite efetuar o Algoritmo de Needleman Wunsch de forma a realizar o alinhamento Global de um par de sequências.
  """

  def __init__(self, s1, s2, g = -4):
    """Contrutor da classe NW

    Args:
        s1 (str): Sequência 1
        s2 (str): Sequência 2
        g (int, optional): Gap Penalty. Defaults to -4.
    """
    self.s1 = s1
    self.s2 = s2
    self.mat = Mat(len(s1) + 1, len(s2) + 1)  
    self.tr  = Mat(len(s1) + 1, len(s2) + 1)
    
    for L in range(len(s1)):                     
      self.mat[L + 1][0] = g * (L + 1)           
      self.tr[L + 1][0]  = 'C'                      

    for C in range(len(s2)):                     
      self.mat[0][C + 1] = g * (C + 1)             
      self.tr[0][C + 1]  = 'E'                        

    for L, x1 in enumerate(s1):
      for C, x2 in enumerate(s2):
        possiveis = [
            self.mat[L  ][C    ] + subst(x1, x2),   
            self.mat[L+1][C    ] + g,               
            self.mat[L  ][C + 1] + g,               
        ]
        dirs = "DEC"

        self.mat[L + 1][C + 1] = max(possiveis)   
        self.tr[L + 1][C + 1] = dirs[possiveis.index(self.mat[L + 1][C + 1])]  

  def rebuild(self):
    """Função que permite obter um possivel alinhamento das sequências

    Returns:
        str: Sequências alinhadas 
    """
    L = len(self.s1)
    C = len(self.s2)
    S1 = ""
    S2 = ""
    
    dirs = {
        'D' : (-1, -1),
        'E' : ( 0, -1),
        'C' : (-1,  0)
    }

    while L > 0 or C > 0:

      DL, DC = dirs[self.tr[L][C]]

      if self.tr[L][C] == "D":
        S1 = self.s1[L - 1] + S1    
        S2 = self.s2[C - 1] + S2
      elif self.tr[L][C] == "E":
        S1 = '-' + S1
        S2 = self.s2[C - 1] + S2
      else:
        S1 = self.s1[L - 1] + S1
        S2 = '-' + S2        

      L += DL
      C += DC

    return S1, S2

 
  def max_score (self):
    """Função que permite obter o max score do alinhamento 

    Returns:
        int: Max Score
    """
    return self.mat[len(self.s1)][len(self.s2)]

  def __repr__(self):
    """Função que permite obter a representação em String da própria classe 

    Returns:
        str: Matriz S (Score) e Matriz T (Trace)
    """
    cols = "-" + self.s2
    lins = "-" + self.s1
    with io.StringIO("") as S:
      print(' ', *cols, sep = '   ', file = S)           
      for L, linha in zip(lins, self.mat):
        print(L, *[f'{x:3d}' for x in linha], file = S)

      print(file = S)

      print(' ', *cols, file = S)
      for L, linha in zip(lins, self.tr):
        print(L, *linha, file = S)

      return S.getvalue()

  def score_mat(self):
    """Função que permite obter a Matriz S (Score)

    Returns:
        str: Matriz S
    """
    cols = "-" + self.s2
    lins = "-" + self.s1
    with io.StringIO("") as S:
      print(' ', *cols, sep = '   ', file = S)           
      for L, linha in zip(lins, self.mat):
        print(L, *[f'{x:3d}' for x in linha], file = S)

      print(file = S)
      return print(S.getvalue())

  def trace_mat(self):
    """Função que permite obter a Matriz T (Trace)

    Returns:
        str: Matriz T (Trace)
    """
    cols = "-" + self.s2
    lins = "-" + self.s1
    with io.StringIO("") as S:
      print(' ', *cols, file = S)
      for L, linha in zip(lins, self.tr):
        print(L, *linha, file = S)

      return print(S.getvalue())

import unittest

class Test_NW(unittest.TestCase):
    """Conjunto de testes para a Classe NW

    """
    
    
    def test_rebuild (self):
        """Testar a função rebuild
        """

        self.assertEqual(NW("ACTG", "ACTG").rebuild(), ('ACTG', 'ACTG'))
        self.assertEqual(NW("ATCAT", "ACAT").rebuild(), ('ATCAT', 'A-CAT'))


    def test_max_score (self):
        """Testar a função max_score
        """

        self.assertEqual(NW("ACTG", "ACTG").max_score(), 24)
        self.assertEqual(NW("ATCAT", "ACAT").max_score(), 18)



unittest.main(argv=[''], exit=False)

class SW:
  """Classe que permite efetuar o Algoritmo de Smith Waterman de forma a realizar o alinhamento Local de um par de sequências.
  """


  def __init__(self, s1, s2, g = -4):
    """Contrutor da classe SW

    Args:
        s1 (str): Sequência 1
        s2 (str): Sequência 2
        g (int, optional): Gap Penalty. Defaults to -4.
    """

    self.s1 = s1
    self.s2 = s2
    self.mat = Mat(len(s1) + 1, len(s2) + 1) 
    self.tr  = Mat(len(s1) + 1, len(s2) + 1)


    for L, x1 in enumerate(s1):
      for C, x2 in enumerate(s2):
        possiveis = [
            self.mat[L  ][C    ] + subst(x1, x2),   
            self.mat[L+1][C    ] + g,               
            self.mat[L  ][C + 1] + g,               
            0]
        dirs = "DEC."

        self.mat[L + 1][C + 1] = max(possiveis)   
        if self.mat[L + 1][C + 1] != 0:
          self.tr[L + 1][C + 1] = dirs[possiveis.index(self.mat[L + 1][C + 1])] 
  
  def max_score (self):
    """Função que permite obter o max score do alinhamento 

    Returns:
        int: Valor do Max Score
    """
    max_score = 0
    for L, x1 in enumerate(self.mat):
      for C, x2 in enumerate(self.mat):
        if max(x1) > max_score:
          max_score = max(x1)
    return max_score

  
  def rebuild(self):
    """Função que permite obter um possivel alinhamento das sequências

    Returns:
        str: Sequências alinhadas 
    """
    max_score = 0
    linha = 0
    coluna = 0
    for L, x1 in enumerate(self.mat):
      for C, x2 in enumerate(self.mat):
        if max(x1) > max_score:
          max_score = max(x1)
          linha = L
        coluna = C
    coluna -= 1
    L = linha
    C = coluna + 1
    S1 = ""
    S2 = ""
    
    dirs = {
        'D' : (-1, -1),
        'E' : ( 0, -1),
        'C' : (-1,  0),
        "." : (0, 0)
    }

    while L >= 0 or C >= 0:

      try:
        DL, DC = dirs[self.tr[L][C]]
        if self.tr[L][C] == "D":
            S1 = self.s1[L - 1] + S1    
            S2 = self.s2[C - 1] + S2
        elif self.tr[L][C] == "E":
            S1 = '-' + S1
            S2 = self.s2[C - 1] + S2
        elif self.tr[L][C] == "C":
            S1 = self.s1[L - 1] + S1
            S2 = '-' + S2 
        else:
          break
    
        L += DL
        C += DC
      except:
        break

    
    return S1, S2


  def __repr__(self):
    """Função que permite obter a representação em String da própria classe 

    Returns:
        str: Matriz S (Score) e Matriz T (Trace)
    """
    
    cols = "-" + self.s2
    lins = "-" + self.s1
    with io.StringIO("") as S:
      print(' ', *cols, sep = '   ', file = S)           
      for L, linha in zip(lins, self.mat):
        print(L, *[f'{x:3d}' for x in linha], file = S)

      print(file = S)

      print(' ', *cols, file = S)
      for L, linha in zip(lins, self.tr):
        print(L, *linha, file = S)

      return S.getvalue()

  def score_mat(self):
    """Função que permite obter a Matriz S (Score)

    Returns:
        str: Matriz S
    """
    cols = "-" + self.s2
    lins = "-" + self.s1
    with io.StringIO("") as S:
      print(' ', *cols, sep = '   ', file = S)           
      for L, linha in zip(lins, self.mat):
        print(L, *[f'{x:3d}' for x in linha], file = S)

      print(file = S)
      return print(S.getvalue())

  def trace_mat(self):
    """Função que permite obter a Matriz T (Trace)

    Returns:
        str: Matriz T (Trace)
    """
    cols = "-" + self.s2
    lins = "-" + self.s1
    with io.StringIO("") as S:
      print(' ', *cols, file = S)
      for L, linha in zip(lins, self.tr):
        print(L, *linha, file = S)

      return print(S.getvalue())

import unittest

class Test_SW(unittest.TestCase):
    """Conjunto de testes para a Classe SW

    """
    
    
    def test_rebuild (self):
        """Testar a função rebuild
        """

        self.assertEqual(SW("ACTG", "ACTG").rebuild(), ('ACTG', 'ACTG'))
        self.assertEqual(SW("ATCATAA", "ATTACAT").rebuild(), ('T-CAT', 'TACAT'))


    def test_max_score (self):
        """Testar a funçao max_score
        """

        self.assertEqual(SW("ACTG", "ACTG").max_score(), 24)
        self.assertEqual(SW("ATCATAA", "ATTACAT").max_score(), 19)

        
unittest.main(argv=[''], exit=False)

def consensus(s1, s2):
    """Função que permite obter o Consensus das 2 Sequências 

    Args:
        s1 (str): Sequencia 1
        s2 (str): Sequencia 2

    Returns:
        str: Consensus das Sequências
    """
    res = ""
    for x, y in zip(s1, s2):
        if x == y:
            res += x
        elif x == '-':
            res += y
        else:
            res += x
    return res
    
class MultipleAlign:
    """Classe que permite realizar o alinhamento múltiplo de várias sequências 
    """
    def __init__(self, lista_seqs):
        """Contrutor da classe MultpleAlign

        Args:
            lista_seqs (list): Lista de sequências 
        """
        self.lista_seqs = lista_seqs
        self.len = len(lista_seqs)



    def __len__ (self):
        """Função que permite implementar a função len()

        Returns:
            int: Numero de Sequências presentes na Classe
        """

        return len(self.lista_seqs)
    

    def __repr__(self):
        """Função que permite representar a Classe como String 

        Returns:
            str: Diferentes Sequências presentes na Classe
        """

        return str(self.lista_seqs)
    
    def __getitem__(self,a):
        """Função que permite implementar a Indexação na Classe

        Args:
            a (int): Index

        Returns:
            str: Sequência correspondente ao Index 
        """
        return self.lista_seqs[a]
    
    
    def align (self):
        """Função que permite realizar o Alinhamento Progressivo das Sequências 

        Returns:
            list: Lista com o alinhamento entre todas as sequências
        """
        lista_final_con = self.lista_seqs.copy()
        lista_final = []
        c = 0
        while c < len(self.lista_seqs)-1:
            lista_temp = lista_final_con[0:2].copy()
            a = NW(lista_temp[0],lista_temp[1]).rebuild()
            con = consensus(a[0],a[1])
            del lista_final_con[0:2]
            lista_final_con.insert(0,con)
            c += 1

        i = len(self.lista_seqs)

        for c in range(i):
            lista_final.append(NW(self.lista_seqs[c],lista_final_con[0]).rebuild()[0])
        

        return lista_final

class Test_MultipleAlign (unittest.TestCase):
    """Conjunto de testes para a Classe MultipleAlign

    """

    def test_align (self):
        """Testar função Align
        """
        self.assertEqual(MultipleAlign(["ATAGC", "AACC", "ATGAC","ATAGC", "AACC", "ATGAC"]).align(), ['ATAGC', 'A-ACC', 'ATGAC', 'ATAGC', 'A-ACC', 'ATGAC'])
        self.assertEqual(MultipleAlign(["TACCCGCT","AACCA"]).align(), ['TACCCGCT', '-A-AC-CA'])


unittest.main(argv=[''], exit=False)