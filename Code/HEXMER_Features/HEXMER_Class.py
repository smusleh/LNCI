# #############################
# Hexamer Feature 
#
import numpy
class HexScore:
    def __init__(self, seqid, seq):
        self.id = seqid
        self.sequence = seq
        self.logscore_dict = self.ReadLogScore()

    def ReadLogScore(self):
        '''return a dict of logscore'''

        logscore_file = "./log_file/human.hexamer.logscore"

        logscore_dict = {}

        with open(logscore_file, "r") as fl:
            for line in fl.readlines():
                line = line.strip()
                logscore_dict[ line.split()[0] ] = float(line.split()[1])

            return logscore_dict
    
    def HexamerGenerate(self, seq):
        '''Generate hexamer'''

        hexamer_list = []
        hexamer_score_list = []
        if(len(seq) > 3):
            num = len(seq) // 3
            for i in range(0,num-1) :
                tmp =  seq[ i*3:(i+2)*3] 
                hexamer_list.append(tmp)
                #print("tmp = ", tmp)
                if  tmp in self.logscore_dict:
                    hexamer_score_list.append( self.logscore_dict[tmp] )
                else:
                    hexamer_score_list.append(0)

        #print("Hex score list = ", hexamer_score_list)
        return hexamer_score_list


    def HexamerScore2 (self, ORF_seq):
        '''compute the hexamer score given a sequence'''
        
        total = 0.0
        log_score = 0.0
        #print("ORF_seq in HexamerScore2", ORF_seq)
        for k in self.HexamerGenerate(ORF_seq):
            #print("k = ",k)
            log_score += k
            total += 1
            #print("sum log score = ", log_score)
            #print("total = ",total)
        try:
            return log_score/total
        except:
            return -1

    def GetORF_UTR(self):
        '''Get ORF and UTR from sequence'''
        STP = {}
        STP[0] = []
        STP[1] = []
        STP[2] = []
        STP[0].append(0)
        STP[1].append(1)
        STP[2].append(2)
        AAnum = len(self.sequence) // 3

        for i in range(0, 3):
            for j in range(0, AAnum):
                tmp = self.sequence[(i+3*j):(i+3*(j+1))]
                if tmp == 'TAG' or tmp == 'TAA' or tmp == 'TGA':
                    STP[i].append(i+3*j)

        ORF = {}

        for i in range(0,3):
            if len(STP[i]) < 2:
                continue
            for j in range(1, len(STP[i])):
                tmpN = int( (STP[i][j] - STP[i][j-1])/3 )
                for k in range(0, tmpN):
                    tmpS = self.sequence[ (STP[i][j-1] + 3*k):(STP[i][j-1] + 3*(k+1)) ]
                    if tmpS == 'ATG':
                        ORF[3*k + STP[i][j-1]] = STP[i][j]
                        break

        # longest ORF
        ORFseq = []
        ORFlen = []
        ORFstart = []
        ORFend = []
        for (k,v) in ORF.items():
            ORFseq.append(self.sequence[k:(v-3)])
            ORFlen.append(v - k)
            ORFstart.append(k)
            ORFend.append(v)

        # ORF doesn't exist
        if ORF:
            ORF_l = ORFseq[numpy.argmax(ORFlen)]
            UTR5 = ''
            UTR3 = ''
            if len(self.sequence[ 0:ORFstart[numpy.argmax(ORFlen)] ]) > 0:
                UTR5 = self.sequence[ 0:ORFstart[numpy.argmax(ORFlen)] ]
            if len(self.sequence[ORFend[numpy.argmax(ORFlen)]:]) > 0:
                UTR3 = self.sequence[ORFend[numpy.argmax(ORFlen)]:]
            return [ORF_l, UTR5, UTR3]
        else:
            return ['', '', '']

    def gen_hex_score (self):
        
        # search ORF
        seq = self.sequence.strip()
        ORF, UTR5, UTR3 = self.GetORF_UTR()
        ave_hexamer_score = self.HexamerScore2(ORF)
        return ave_hexamer_score