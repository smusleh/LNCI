# ################################################################
# Fickett Class to generate Fickett features for a given sequence
#

class Fickett:
    def __init__(self,seqid,seq):
        self.id = seqid
        self.sequence = seq
    
    def print_id(self):
        print("seq id = " + str(self.id))
    
    def print_seq(self):
        print(str(self.sequence))
    
    def print_len(self):
        print("Seq length = ", str(len((self.sequence))))
    
    def GetBasePositionScore(self, seq, base):
        '''compute base(A, C, G, T) position score based
            on Fickett's methods
        '''
        base_num = [0, 0, 0]
        for i in range(0, len(seq), 3):
            if seq[i] == base:
                base_num[i%3] += 1
            if i+1 == len(seq):
                break
            if seq[i+1] == base:
                base_num[(i+1)%3] += 1
            if i+2 == len(seq):
                break
            if seq[i+2] == base:
                base_num[(i+2)%3] += 1
        
        base_pos_score = max(base_num) * 1.0 / (min(base_num) + 1)

        return base_pos_score

           
    # Fickett feature
    def gen_fic_fea (self) :
        A_pos_fea = self.GetBasePositionScore(self.sequence, 'A')
        C_pos_fea = self.GetBasePositionScore(self.sequence, 'C')
        G_pos_fea = self.GetBasePositionScore(self.sequence, 'G')
        T_pos_fea = self.GetBasePositionScore(self.sequence, 'T')
        return [A_pos_fea,C_pos_fea,G_pos_fea,T_pos_fea]
