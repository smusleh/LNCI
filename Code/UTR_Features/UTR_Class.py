import numpy
class UTRkmer:
    def __init__(self, seqid, seq):
        self.id = seqid
        self.sequence = seq


    def GetORF_UTR(self,seq):
        '''Get ORF and UTR from sequence'''
        STP = {}
        STP[0] = []
        STP[1] = []
        STP[2] = []
        STP[0].append(0)
        STP[1].append(1)
        STP[2].append(2)
        AAnum = len(seq) // 3

        for i in range(0, 3):
            for j in range(0, AAnum):
                tmp = seq[(i+3*j):(i+3*(j+1))]
                if tmp == 'TAG' or tmp == 'TAA' or tmp == 'TGA':
                    STP[i].append(i+3*j)

        ORF = {}

        for i in range(0,3):
            if len(STP[i]) < 2:
                continue
            for j in range(1, len(STP[i])):
                tmpN = int( (STP[i][j] - STP[i][j-1])/3 )
                for k in range(0, tmpN):
                    tmpS = seq[ (STP[i][j-1] + 3*k):(STP[i][j-1] + 3*(k+1)) ]
                    if tmpS == 'ATG':
                        ORF[3*k + STP[i][j-1]] = STP[i][j]
                        break

        # longest ORF
        ORFseq = []
        ORFlen = []
        ORFstart = []
        ORFend = []
        for (k,v) in ORF.items():
            ORFseq.append(seq[k:(v-3)])
            ORFlen.append(v - k)
            ORFstart.append(k)
            ORFend.append(v)

        # ORF doesn't exist
        if ORF:
            ORF_l = ORFseq[numpy.argmax(ORFlen)]
            UTR5 = ''
            UTR3 = ''
            if len(seq[ 0:ORFstart[numpy.argmax(ORFlen)] ]) > 0:
                UTR5 = seq[ 0:ORFstart[numpy.argmax(ORFlen)] ]
            if len(seq[ORFend[numpy.argmax(ORFlen)]:]) > 0:
                UTR3 = seq[ORFend[numpy.argmax(ORFlen)]:]
            return [ORF_l, UTR5, UTR3]
        else:
            return ['', '', '']


    # def gen_ORF(self):
    #     # search ORF
    #     seq = self.sequence.strip()
    #     ORF, UTR5, UTR3 = self.GetORF_UTR(seq)
    #     transcript_len = len(self.sequence)

    def GetGC_Content(self, seq):
        '''calculate GC content of sequence'''
        A, C, G, T = 1e-9, 1e-9, 1e-9, 1e-9 
        for i in range(0, len(seq)):
            if seq[i] == 'A':
                A += 1
            elif seq[i] == 'C':
                C += 1
            elif seq[i] == 'G':
                G += 1
            elif seq[i] == 'T':
                T += 1

        GC = (G + C) / (A + C + G + T)

        return GC 

    def gen_utrs53(self):
        
        seq = self.sequence.strip()
        ORF, UTR5, UTR3 = self.GetORF_UTR(seq)
        transcript_len = len(self.sequence)
        
        UTR5_ratio = len(UTR5) * 1.0 / transcript_len
        UTR5_gc_content = self.GetGC_Content(UTR5)

        UTR3_ratio = len(UTR3) * 1.0 / transcript_len
        UTR3_gc_content = self.GetGC_Content(UTR3)

        return [UTR5_ratio, UTR5_gc_content,  UTR3_ratio, UTR3_gc_content]

