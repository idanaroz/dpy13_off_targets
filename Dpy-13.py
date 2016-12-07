from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

input_fasta = "E:\google drive\year 3\BioInformatics project\F30B5.1spliced_+_UTR.fasta"
output_fasta = "E:\google drive\year 3\BioInformatics project\multiFasta"


def create_input_multi_fasta_file(input_file ):
    fasta_record = SeqIO.read(open(input_file), 'fasta')
    dpy_13_frags = []
    n= len(fasta_record.seq)
    for i in range (0,n-25):
        start = i
        end = i+25
        frag = fasta_record.seq[start:end]
        record = SeqRecord(frag, 'Dpy-13:%i' % (i + 1)+"_"+str(i+25), '', '')
        dpy_13_frags.append(record)

    output_handle = open("E:\google drive\year 3\BioInformatics project\Dpy-13_multiFasta.fasta", "w")
    SeqIO.write(dpy_13_frags, output_handle, "fasta")
    output_handle.close()

result =  create_input_multi_fasta_file(input_fasta )
