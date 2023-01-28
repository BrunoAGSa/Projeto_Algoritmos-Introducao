class Seq:
    """Esta Classe permite fazer X
    """

    def __init__(self, seq):
        """_summary_

        Args:
            seq (_type_): _description_
        """
        self.seq = str(seq.upper())

    def __len__ (self): 
        """
        
        """ 
        return len(self.seq)
    
    def __str__(self):
        """  
        """       
        return self.seq

    def get_type (self):
        """_summary_

        Returns:
            _type_: _description_
        """
        if len([c for c in self.seq if c not in "ACGT"]) == 0:
            return "DNA"

        elif len([c for c in self.seq if c not in "ACGU"]) == 0:
            return "RNA"

        elif len([c for c in self.seq if c not in "BDEFHIJKLMNOPQRSUVWXYZ"]) == 0:
            return "PROTEIN"
        
        else:
            return "Not valid"

    def count_nucleotides (self):
        """_summary_

        Returns:
            _type_: _description_
        """
       
        if Seq(self.seq).get_type() == "DNA":
            A = self.seq.count("A")
            C = self.seq.count("C")
            G = self.seq.count("G")
            T = self.seq.count("T") 
            print (" A   C   G   T ")
            print (" |   |   |   |")
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
        """_summary_

        Returns:
            _type_: _description_
        """

        if Seq(self.seq).get_type() == "DNA":
            A = self.seq.count("A") / len(self.seq)
            C = self.seq.count("C") / len(self.seq)
            G = self.seq.count("G") / len(self.seq)
            T = self.seq.count("T") / len(self.seq)
            return [round(A, 3), round(C, 3), round(G, 3), round(T, 3)]
        
        elif Seq(self.seq).get_type() == "DNA":
            A = self.seq.count("A") / len(self.seq)
            C = self.seq.count("C") / len(self.seq)
            G = self.seq.count("G") / len(self.seq)
            U = self.seq.count("U") / len(self.seq)
            return [round(A, 3), round(C, 3), round(G, 3), round(U, 3)]
        else:
            return "Not DNA or RNA"

