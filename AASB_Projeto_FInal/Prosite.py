def prosite(seq, pattern): 
                        
    """Função que permite converter uma expressão do padrão Prosite numa expressão regular

    Returns:
        str: Expressão regular
    """
    
    from re import search
    er = pattern.replace("-","")
    er = er.replace("x",".")
    er = er.replace("(","{")
    er = er.replace(")","}")
    motif = search(er, seq)

    if (motif != None):
        zona = motif.span()              
        return ("Motifs: {}" .format(seq[zona[0]:zona[1]]))
    else :
        return "Not found"

import unittest

class Test_Prosite(unittest.TestCase):
    """Conjunto de testes para a classe Prosite
    """

    def test_prosite(self):
        """Testar a função prosite
        """
        self.assertEqual(prosite("HKMMLASCKHLLCLKCIVKLG","C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]"),"Motifs: CKHLLCLKCI")
        self.assertEqual(prosite("HKMMLLCLKCIVKLG","C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]"),"Not found")



unittest.main(argv=[''], exit=False)