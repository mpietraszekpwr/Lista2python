'''
Aby przeprowdzić transkrypcję mRNA na białko potrzeby jest słownik który tłumaczy kod : 3 zasady -> aminokwas.Sktrypt, który generuje ten słownik jest na stronie:
https://www.petercollingridge.co.uk/tutorials/bioinformatics/codon-table/
kod został napisany
'''
def get_codon_table():
    bases = "UCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    return dict(zip(codons, amino_acids))


class DNASequence:
    valid_chars = set("ATGC")
    complement_trans = str.maketrans("ATGC", 'TACG')
    transcribe_trans = str.maketrans("T", "U")

    def __init__(self, id, data):
        self.id = id
        self.data = data
        return

    def __len__(self):
        return len(self.data)

    def __str__(self):
        return ''.join(('>', self.id, '\n', self.data))

    def mutate(self, pos, value):
        self.data = ''.join((self.data[0:pos], value, self.data[pos + 1:]))

    def find_motif(self, motif):
        return self.data.index(motif)

    def complement(self):
        return self.data.translate(self.complement_trans)[::-1]

    def transcribe(self):
        return RNASequence(self.id, self.data.translate(self.transcribe_trans))

    def is_valid(self):
        for k in self.data:
            if k not in self.valid_chars:
                break
        else:
            return True
        return False


class RNASequence:
    valid_chars = set("AUGC")
    codon_table = get_codon_table()

    def __init__(self, id, data):
        self.id = id
        self.data = data
        return

    def __len__(self):
        return len(self.data)

    def __str__(self):
        return ''.join(('>', self.id, '\n', self.data))

    def mutate(self, pos, value):
        self.data = ''.join((self.data[0:pos], value, self.data[pos + 1:]))

    def find_motif(self, motif):
        return self.data.index(motif)

    def transcribe(self, from_start=True, to_stop=True):
        """
        from_start:
            True  transkrypcja od pierwszego "AUG"
            False transkrypcja od początku
        to_stop:
            True  transkrypcja do pierwszego stopu (UAA,UAG,UGA)
            False transkrypcja do końca. Stop kodowany jako *
        """
        beg = self.data.index("AUG") if from_start else 0
        temp = []
        for k in range(beg, len(self.data), 3):
            kodon = self.data[k:k + 3]
            if len(kodon) < 3:
                break
            a = self.codon_table[kodon]
            if (a == '*') and to_stop:
                break
            temp.append(a)
        seq = ''.join(temp)
        return ProteinSequence(self.id, seq)

    def is_valid(self):
        for k in self.data:
            if k not in self.valid_chars:
                break
        else:
            return True
        return False


class ProteinSequence:
    valid_chars = set('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')

    def __init__(self, id, data):
        self.id = id
        self.data = data
        return

    def __len__(self):
        return len(self.data)

    def __str__(self):
        return ''.join(('>', self.id, '\n', self.data))

    def mutate(self, pos, value):
        self.data = ''.join((self.data[0:pos], value, self.data[pos + 1:]))

    def find_motif(self, motif):
        return self.data.index(motif)

    def is_valid(self):
        for k in self.data:
            if k not in self.valid_chars:
                break
        else:
            return True
        return False

d=DNASequence('pierwsza','GCT') ## powinna kodować alaninę oznaczenie A
print(d)
print(len(d))
print(d.is_valid())

r=d.transcribe() ## transkrypcja do RNA
print(r)

p=r.transcribe(False,False) ## transkrypcja do białka
print(p)

p=r.transcribe(False,False)
print(p)

print(d)
d.mutate(2,'s')
print(d)
print(len(d))
print(d.is_valid())

# Przykładowa sekwencja z internetu pobrana
with open("C:\\Users\\pietr\\Desktop\\Pwr\\Języki programowania\\NC_005816.fna") as f:
    id=f.readline().strip()[1:]
    s=f.readlines()
seq=''.join(l.strip() for l in s)
yes=DNASequence(id,seq)

r = yes.transcribe()
p = r.transcribe(True, True)
print(p)
