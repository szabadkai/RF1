from Bio import SeqIO 
import re

# A program megallapitja, hogy van e transzlacios readthrouhg indukalo stop 
# kontextusa a betaplalt mRNS-eknek. Ha mar ezzel megvan, megnezi, hol a kovetkezo
# in Frame stop kodon.

table = 11
min_pro_len = 11


'''
	a funkcio Nyilt leovasasi kereteket (ORF) keres minden Frameben, 
	es ezek kozul mindig a leghosszabbakat adja vissza egy listaban.
									BioPython Tutorial page nyoman 
			(http://biopython.org/DIST/docs/tutorial/Tutorial.html)
'''
def find_orfs(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3                        
                    answer.append((start, end, strand,
                                   trans[aa_start:aa_end]))
                aa_start = aa_end+1
    answer.sort()
    return answer

def trace_stop(orf_list):
	'''
	a find_orf altal talalt ORFek vegeit keresi ki es tovabbadja a leghosszabb
	ORF stop, start es srtand(melyik szal) ertekeit.
	'''
	max_len = 0	
	for start, end, strand, pro in orf_list: 
		if len(pro)>max_len: #borzaszto szofisztikalt algoritmus a legnagyobb ertek kivalasztasahoz.
			max_len = len(pro)
			(max_start, max_stop,  max_strand) = (start, end, strand)
		else:
			pass
	max_len = 0	#erre lehet, hogy nincs szukseg, mondjuk artani se art...
	return(max_start,max_stop,max_strand)
	

def trace_stop_context(record):
	'''
	ez itt a trace_stop altal kivalasztott stop kontextust kikeresi az
	RNS szekvenciabol, majd tovabbadja a szekvenciat (stp+6nt formatumban)
	'''
	orf_list = find_orfs(record.seq, table, min_pro_len)
	(max_start , max_stop, max_strand) = trace_stop(orf_list)

	if max_strand == -1:
		return(record.seq[max_start-6:max_start+3].reverse_complement(),'-1', max_stop, max_start)
	else:
		return(record.seq[max_stop-3:max_stop+6],'+1', max_stop, max_start)


handle = open('longest.fa', "rU")

for record in SeqIO.parse(handle, "fasta") :
	try:
		(sequence, strand, max_stop, max_start) = trace_stop_context(record)
		if strand == '-1':
			if 10<len(record.seq[:max_start].reverse_complement()) :
				print('>' + record.id + '\n' + record.seq[:max_start].reverse_complement())
		else:
			if 10<len(record.seq[max_stop:]):
				print('>'+record.id +'\n'+ record.seq[max_stop:])
	except:
		pass
