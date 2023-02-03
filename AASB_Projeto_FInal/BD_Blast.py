import mysql.connector as SQLC  

DataBase = SQLC.connect(
   host ="127.0.0.1",
   user ="root",
   password ="login123",
   database = "Base Dados IAP",
   auth_plugin="mysql_native_password"
)
DataBase.autocommit= True
cur = DataBase.cursor()

def connect_database (host ="127.0.0.1", user ="root",password ="login123", database = "Base Dados IAP"):
    """Funçao que permite ao utilizador conectar à base de dados (Base Dados IAP)

    Args:
        host (str, optional):  Host . Defaults to "127.0.0.1".
        user (str, optional):  User . Defaults to "root".
        password (str, optional):  Password . Defaults to "login123".
        database (str, optional):  Database . Defaults to "Base Dados IAP".
    """

    DataBase = SQLC.connect(

        host = host,
        user = user,
        password = password,
        database = database,
        auth_plugin="mysql_native_password"
    )
    DataBase.autocommit= True
    cur = DataBase.cursor()

def get_all_organisms ():
    """Função que permite obter todos os organismos presentes na Base de Dados Predefinida

    Returns:
        lista: Contêm nome de Organimos 
    """
    connect_database()
    cur.execute("""select Organism from Genbank""")
    organisms_DB = cur.fetchall()

    return organisms_DB


def get_seq_from_organism (organism):
    """Função que permite obter a sequência do Organismo pretendido

    Args:
        organism (str): Organismo Presente na base de dados 

    Returns:
        str: Sequência do Organismo
    """
    connect_database()
    cur.execute(f""" select Sequence from Sequences
                    join Genbank on Sequences.ID_version_seq = Genbank.ID_version_seq
                    where Organism = '{organism}'
    """)

    sequence = cur.fetchall()
    
    return sequence[0][0]

def get_all_seqs():
    """Função que permite obter todas as sequências presentes na Base de Dados

    Returns:
        list: Lista com todas sequências presentes na Base de Dados
    """
    connect_database()
    cur.execute(f""" select Sequence from Sequences
    """)

    sequence = cur.fetchall()

    return sequence

import unittest

class Test_DataBase (unittest.TestCase):
    """Conjunto de testes para a DataBase
    """

    def test_database (self):
        """Testar a integridade da Base de Dados
        """
        self.assertEqual(len(get_all_organisms()),len(get_all_seqs()))


    def test_database_2 (self):
        """Testar a integridade da Base de Dados
        """
        self.assertEqual(type(get_seq_from_organism('Sulfolobus spindle-shaped virus 1')), str)


unittest.main(argv=[''], exit=False)

class Blast ():
    """Classe que permite realizar um Blast Simplificado 
    """
    def __init__(self, query, seq, w = 3):
        """Construtor do Blast Simplificado 

        Args:
            query (str): Query 
            seq (str): Seq
            w (int, optional): W . Defaults to 3.
            
        """
        self.seq = seq
        self.query = query
        self.w = w

    def __len__(self):
        return len(self.query)

    
    def query_map(query, w):
        """Função que permite devolver um dicionário em que as chaves são as sequências e os valores são uma lista dos índices

        Args:
            query (str): Query
            w (int): W

        Returns:
            dict: Dicionário em que as chaves são as sequências e os valores são uma lista dos índice
        """
        
        tam = len(query)
        res = {}
        for chave, offset in [(query[p : p + w], p) for p in range(0, tam - w +1)]:
            if chave not in res:
                res[chave] = []
            res[chave].append(offset)
        return(res)


    def get_all_offsets(s1, s2):  
        """Função Auxiliar que permite devolver todos os índices das ocorrências de uma substring numa string

        Args:
            s1 (str): Substring
            s2 (str): String

        Returns:
            dict: Dicionário com todas as ocorrèncias da substring na string
        """
        w = len(s1)

        return[p for p in range(0, len(s2) - w + 1) if s2[p : p + w] == s1]


    def hits(qm, seq):
        """Função que permite gerar uma lista de hits

        Args:
            qm (dict): Dicionário proveninete da Função query_map
            seq (str): Sequência Proveniente da base de dados
        
        Returns:
            list: Lista de hits em que cada elemento é um tuplo com os indices
        """
              
        return[(o_query, o_seq) for chave, offsets in qm.items() 
        for o_query in offsets 
        for o_seq in Blast.get_all_offsets(chave, seq)]


    def extend_hit_any_direction(query, seq, o1, o2, dir):
        """Função auxliar que permite extender o hit em ambas as direções

        """
        matches = 0
        count = 0
        while o1 >= 0 and o2 >= 0 and o1 < len(query) and o2 < len(seq):
            matches += 1 if query[o1] == seq[o2] else 0
            count += 1
            if 2 * matches < count:
                return o1, o2, matches, count
            o1 += dir
            o2 += dir
        return o1 - dir, o2 - dir, matches, count


    def extend_hit(query, seq, hit, w):
        """Função que permite estende um hit em cada direção se o nº de matches for de pelo menos metade do tamanho da extensão

        Args:
            query (str): Query
            seq (str): Seq
            hit (tuple): Hit
            w (int): w

        Returns:
            tuple: Tuplo com o índice do início do hit estendido na query, na sequência, o tamanho e o nºde matches
        """

        o1, o2 = hit
        left  = Blast.extend_hit_any_direction(query, seq, o1 - 1, o2 - 1, -1)
        right = Blast.extend_hit_any_direction(query, seq, o1 + w, o2 + w, +1)

        O1, O2, ML, SL = left
        _,   _, MR, SR = right
        return O1, O2, w + SL + SR, ML + w + MR
    

    def best_hit (self):
        """Função que permite devolver a extensão de maior score (no caso de empate, deverá devolver a de menor tamanho que aparece primeiro)

        Returns:
            tuple: Tuplo com a extensão que possui maior score
        """
        qm = Blast.query_map(self.query,self.w)
        lista_hits = Blast.hits(qm, self.seq)
        all_scores =[Blast.extend_hit(self.query, self.seq, c, self.w) for c in lista_hits]
        all_scores.sort(key=lambda c: c[2] and c[3], reverse=True)
        if all_scores == []:
            return print("Not Found")
        else:
            return all_scores [0]

import unittest

query = "AATATAT"
seq = "AATATGTTATATAATAATATTT"
w = 3

class Test_Blast (unittest.TestCase):
    """Estes Testes foram baseados nos exemplos fornecidos pelo Professor Rui (inputs/outputs) durantes
    as aulas teóricas 
    """

    def test_query_map (self):
        """Testar a função query_map
        """
        self.assertEqual(Blast.query_map(query,w),{'AAT': [0], 'ATA': [1, 3], 'TAT': [2, 4]})
    

    def test_get_all_offsets (self):
        """Testar a função get_all_offsets
        """
        self.assertEqual(Blast.get_all_offsets(query,seq), [])
    
    
    def test_hits (self):
        """Testar a função test_hits
        """
        self.assertEqual(Blast.hits(Blast.query_map(query,w),seq),[
            (0, 0), (0, 12), (0, 15), (1, 1), (1, 8),
            (1, 10), (1, 13), (1, 16), (3, 1), (3, 8), (3, 10),
            (3, 13), (3, 16), (2, 2), (2, 7), (2, 9), (2, 17),
            (4, 2), (4, 7), (4, 9),(4, 17)])


    def test_extend_hit_any_direction_left (self):
        """Testar a função extend_hit_any_direction
        """
        self.assertEqual(Blast.extend_hit_any_direction(query, seq, 1, 16, -1),(0, 15, 2, 2))

    
    def test_extend_hit_any_direction_right (self):
        """Testar a função extend_hit_any_direction
        """
        self.assertEqual(Blast.extend_hit_any_direction(query, seq, 1, 16, 1),(6, 21, 5, 6))

    def test_extend_hit (self):
        """Testar a função extend_hit
        """
        self.assertEqual(Blast.extend_hit(query, seq, (1, 16), w), (0, 15, 7, 6))

    
    def test_best_hit (self):
        """Testar a função best_hit
        """
        self.assertEqual(Blast(query, seq, w).best_hit(), (0, 0, 7, 6))

unittest.main(argv=[''], exit=False)

