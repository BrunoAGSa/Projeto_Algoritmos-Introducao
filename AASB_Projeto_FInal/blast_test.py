import mysql.connector as SQLC 
 
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

    cur.execute(f""" select Sequence from Sequences
                    join Genbank on Sequences.ID_version_seq = Genbank.ID_version_seq
                    where Organism = '{organism}'
    """)

    sequence = cur.fetchall()
    
    return sequence[0][0]

def get_all_seqs():
    """Função que permite obter todas as sequências presentes na Base de Dados

    Returns:
        list: Lista com as sequências presentes na Base de Dados
    """
    cur.execute(f""" select Sequence from Sequences
    """)

    sequence = cur.fetchall()

    return sequence



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


    def get_all_offsets(s1, s2):  # s1 --> Padrão    s2 --> seq_analisar
        """Função Auxiliar que permite devolver todos os índices das ocorrências de uma substring numa string

        Args:
            s1 (str): _description_
            s2 (str): _description_
        """
        w = len(s1)

        return[p for p in range(0, len(s2) - w + 1) if s2[p : p + w] == s1]


    def hits(qm, seq):      # qm --> query_map(query, w)   
        return[(o_query, o_seq) for chave, offsets in qm.items() 
        for o_query in offsets 
        for o_seq in Blast.get_all_offsets(chave, seq)]


    def extend_hit_any_direction(query, seq, o1, o2, dir):
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

        o1, o2 = hit
        left  = Blast.extend_hit_any_direction(query, seq, o1 - 1, o2 - 1, -1)
        right = Blast.extend_hit_any_direction(query, seq, o1 + w, o2 + w, +1)

        O1, O2, ML, SL = left
        _,   _, MR, SR = right
        return O1, O2, w + SL + SR, ML + w + MR
    

    def best_hit (self):
        qm = Blast.query_map(self.query,self.w)
        lista_hits = Blast.hits(qm, self.seq)
        all_scores =[Blast.extend_hit(self.query, self.seq, c, self.w) for c in lista_hits]
        all_scores.sort(key=lambda c: c[2] and c[3], reverse=True)
        if all_scores == []:
            return print("Not Found")
        else:
            return all_scores [0]