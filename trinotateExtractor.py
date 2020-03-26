

'''Extract proteins and nucleotide sequence from a transdecoder PEP 
or from Trinity assembly file based on trinotate annotation.  Software searches 
Trinotate blastx columns for keywords, extracts Trinity ids, then finds Trinity
ids in pep or fasta file.
'''
 

import re, pandas

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-k', '--keyword', help="Keyword contained in gene name")
parser.add_argument('-t', '--trinotate', help="Trinotate XLS file")
parser.add_argument('-f','--file', help="Your fasta or pep file")
parser.add_argument('-o', '--output', help="Name to ouputnew file")
args = parser.parse_args()

class Extractor(object):
    def __init__(self):
    
        self.keyword = args.keyword
        self.trinotate = pandas.read_table(args.trinotate)
        self.file = open(args.file,'r')
        self.newFile = open(args.output, 'w')
        self.trinLib = {}
        self.isoformList = []

        self.geneDescriptionRegex = r'\S+ Full=([\S\s]*?);'
        self.transcriptRegex = r"(TRINITY_DN[0-9]+_c[0-9]_g[0-9]_i[0-9])"

    def makeLists(self):
        for entry in self.trinotate.index:
            try:
                sprot = self.trinotate['sprot_Top_BLASTX_hit'][entry]
                geneName = re.search(self.geneDescriptionRegex, sprot)
                gene = geneName.groups()[0]
                if self.keyword in gene:
                    self.trinLib[self.trinotate['transcript_id'][entry]] = gene
                    self.isoformList.append(self.trinotate['transcript_id'][entry])
            except AttributeError:
                continue
    
    def writeNewFile(self):
        protWrite = False
        for line in self.file:
            if line.startswith('>') and protWrite == True:
                protWrite = False
            try:
                if line.startswith('>'):
                    found = 0
                    output = re.search(self.transcriptRegex, line)
                    result = output.groups()
                    if result[0] in self.isoformList:
                        self.newFile.write(line)
                        found = 1
                elif found == 1:
                    self.newFile.write(line)
                    protWrite = True
            except AttributeError:
                continue 

        self.file.close()
        self.newFile.close()    
    
      
    def makeItHappen(self):
        self.makeLists()
        self.writeNewFile()
            

if __name__ == "__main__":
    obj = Extractor()
    obj.makeItHappen()