import json

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def blast_to_multifasta_input(input_fasta_file, output_file_xml):
    """"""
    """
    :type input_fasta_file: str
    :type output_file_xml: str
    :param: input_fasta_file: string ,  Should be the multifasta file (with multi queries)
    :param: output_file_xml: string ,  output xml file (would need to be analyzed later)
    :return: string: output xml file (would need to be analyzed later)
    """
    save_file = open(output_file_xml, "w")
    for record in SeqIO.parse(input_fasta_file, format="fasta"):
        result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta") \
                                       , results_file="E:\google drive\year 3\BioInformatics project\myTrial!!!!.fasta" \
                                       , hitlist_size=50, entrez_query="txid6239",
                                       expect=40)

        save_file.write(result_handle.read())
        result_handle.close()
    save_file.close()
    return output_file_xml


def analyze_blast_xml_results(my_blast_xml, outputfile):
    """
    Analyze the xml file with blast results\n
    :type my_blast_xml: str
    :type outputfile: str
    :param: my_blast_xml: the output from the previous function (blast_to_multifasta_input)
    :param: outputfile: file that will include all the blast results
    :return: None
    """
    f = open(outputfile, "w")  # Use "a" instead of "w" to append to file
    result_handle = open(my_blast_xml)  # we can just open the saved file for input
    blast_records = NCBIXML.parse(result_handle)
    blast_records = list(blast_records)
    for blast_record in blast_records:
        f.write('#####Query= ' + blast_record.query + "\n")
        n = len(blast_record.alignments)
        for i in range(n):
            # if( ("Caenorhabditis elegans" not in blast_record.descriptions[i].title) and ("C.elegans" not in blast_record.descriptions[i].title )  ):
            #     continue;
            f.write('***Alignment_' + str(i) + "\n")
            f.write("title= " + blast_record.descriptions[i].title + "\n")

            f.write("accession= " + blast_record.alignments[i].accession + "\n")
            f.write("hit_def= " + blast_record.alignments[i].hit_def + "\n")
            f.write("hit_id= " + blast_record.alignments[i].hit_id + "\n")

            f.write("-hsps\n")
            if (len(blast_record.alignments[0].hsps) > 1):
                f.write('$$$$$$$$$$$$$$$$$$$ hsp include more than 1 segment pair $$$$$$$$$$$$$$$$$$$' + "\n")
            hsp = blast_record.alignments[i].hsps[0]
            f.write('align_length=' + str(hsp.align_length) + "\n")
            f.write('score=' + str(hsp.score) + "\n")
            f.write('bits=' + str(hsp.bits) + "\n")
            f.write('expect=' + str(hsp.expect) + "\n")
            f.write('frame=' + str(hsp.frame) + "\n")
            f.write('gaps=' + str(hsp.gaps) + "\n")
            f.write('identities=' + str(hsp.identities) + "\n")
            f.write('query=' + str(hsp.query) + "\n")
            f.write('match=' + str(hsp.match) + "\n")
            f.write('sbjct=' + str(hsp.sbjct) + "\n")
            f.write('sbjct_start=' + str(hsp.sbjct_start) + "\n")
            f.write('sbjct_end=' + str(hsp.sbjct_end) + "\n")
    result_handle.close()
    f.close()


def hsps_str_repr(hsps):
    """C:\Python34\Lib\site-packages\Bio\Blast\Record.py"""
    lines = ["Score %i (%i bits), expectation %0.1e, alignment length %i"
             % (hsps.score, hsps.bits, hsps.expect, hsps.align_length)]
    if hsps.align_length < 50:
        lines.append("Query:%s %s %s" % (str(hsps.query_start).rjust(8),
                                         str(hsps.query),
                                         str(hsps.query_end)))
        lines.append("               %s"
                     % (str(hsps.match)))
        lines.append("Sbjct:%s %s %s" % (str(hsps.sbjct_start).rjust(8),
                                         str(hsps.sbjct),
                                         str(hsps.sbjct_end)))
    else:
        lines.append("Query:%s %s...%s %s"
                     % (str(hsps.query_start).rjust(8),
                        str(hsps.query)[:45],
                        str(hsps.query)[-3:],
                        str(hsps.query_end)))
        lines.append("               %s...%s"
                     % (str(hsps.match)[:45],
                        str(hsps.match)[-3:]))
        lines.append("Sbjct:%s %s...%s %s"
                     % (str(hsps.sbjct_start).rjust(8),
                        str(hsps.sbjct)[:45],
                        str(hsps.sbjct)[-3:],
                        str(hsps.sbjct_end)))
    return "\n".join(lines)


def processed_dic(txt_file_json):
    """take 'full_dic' format in format : [{"val": [], "key": ["NM_001026753", "2", "-", "7765795", "7771965"]}, {"val": [], "key": ...
        and create 2 files DPY-13-homologs.txt with clear results of blast
     """
    with open(txt_file_json, 'r') as content_file:
        with open('FINAL RESULTS!.txt', 'w') as outfile:
            with open('FINAL RESULTS short.txt', 'w') as outfile_short:
                content = content_file.read()
                data = json.loads(content)
                # outfile.write("key" + 13*"\t" + "#hits" +"\t\t" +"hsps" +"\n")
                for elem in data:
                    key = str(elem["key"]).strip('[]')
                    value = elem["val"]
                    length = len(value)
                    if length > 0:
                        outfile.write("$KEY:" + "\t" + key + "\n")
                        outfile.write("#HITS: " + str(length) + "\n")

                        outfile_short.write("$KEY:" + "\t" + key + "\n")
                        outfile_short.write("#HITS: " + str(length) + "\n")

                        for j in range(length):
                            outfile.write("-hsps " + str(j) + ":" + "\n")
                            outfile.write(str(value[j]).strip('[]') + "\n")

                        outfile.write("------------------------------------------------------------------------\n")
                        # outfile.write(key + "\t\t" + str(length) + "\t\t")
                        # outfile.write(key + "\t\t" +str(length) +"\t\t")
                        # outfile.write(value + "\n")
    content_file.close()
    outfile.close()
    outfile_short.close()


def parse_refseq_to_dic(refseq_file):
    """
    :param refseq_file: path to refseq file (str)
    :return: dictionary
     E.g. {(name, chromosome_num,strand , txstart , txend ) : [] }
    """
    chromosome_dict = {"chrI": "1", "chrII": "2", "chrIII": "3", "chrIV": "4", "chrV": "5", "chrX": "10"}
    d = {}
    with open(refseq_file) as f:
        next(f)
        for line in f:
            list = [s.strip() for s in line.split("\t")]
            chrom_number = chromosome_dict[list[2]]
            tup = (list[1], chrom_number, list[3], list[4], list[5])
            d[tup] = []

    return d;


def add_blastRecord_to_refseq_dic(refseq_dic, alignment, query_name):
    """
    :param refseq_dic:
    :param alignment:
    :param query_name:
    fill the refseq_dic with values. the values are strings - detailed alignments to the keys.
    """
    chrom_dict = {"CHROMOSOME_I": "1", "CHROMOSOME_II": "2", "CHROMOSOME_III": "3", "CHROMOSOME_IV": "4",
                  "CHROMOSOME_V": "5", "CHROMOSOME_X": "10"}
    if alignment.hit_def != "CHROMOSOME_MtDNA":
        chrom_string_num = chrom_dict[alignment.hit_def]  # to fix format
    else:
        chrom_string_num = "666666"
    # handle each hsps in the alignment
    for k in range(len(alignment.hsps)):
        hsp = alignment.hsps[k]
        strand = str(hsp.frame[1])  # check it
        sbjct_start = hsp.sbjct_start
        sbjct_end = hsp.sbjct_end

        for key in refseq_dic:
            dic_chrom_num = key[1];
            dic_strand = key[2];
            dic_txstart = int(key[3]);
            dic_txend = int(key[4])

            if chrom_string_num != dic_chrom_num:
                continue
            if (dic_strand == "+" and strand == "1"):
                if dic_txstart <= sbjct_start and sbjct_start <= dic_txend:
                    # BINGO !
                    refseq_dic[key] += ["query_name:" + query_name + "\n" + hsps_str_repr(hsp)]
            elif (dic_strand == "-" and strand == "-1"):
                if dic_txstart <= sbjct_start and sbjct_start <= dic_txend:
                    # BINGO !
                    refseq_dic[key] += ["query_name:" + query_name + "\n" + hsps_str_repr(hsp)]


def remap_keys(dict_with_tuple_key):
    """return list that each contain dictionary with one key and val
        E.g: input d ata ={(5, 6, 7): 'qers', (1, 2, 3): 'asdfdf'}
             become : [{"key": [5, 6, 7], "val": "qers"}, {"key": [1, 2, 3], "val": "asdfdf"}]"""
    return [{'key': k, "val": v} for k, v in dict_with_tuple_key.items()]


def create_full_dic(refseq_file, dpy13_blastn_wordsize7_xml):
    """
    :param refseq_file: path to refseq file (str)
    :param dpy13_blastn_wordsize7_xml: path to xml blast result xml file (E.G. "/past used files dpy13_blastn_wordsize7.xml")
    :return:
    """
    refseq_dic = parse_refseq_to_dic(refseq_file)
    result_handle = open(dpy13_blastn_wordsize7_xml)  # we can just open the saved file for input
    blast_records = NCBIXML.parse(result_handle)
    # blast_records=list(blast_records)
    for blast_record in blast_records:
        n = len(blast_record.alignments)
        for i in range(n):
            alignment = blast_record.alignments[i]

            # if alignment.hit_def == "CHROMOSOME_III":
            add_blastRecord_to_refseq_dic(refseq_dic, alignment, blast_record.query)

    with open('full_dic.txt', 'w') as outfile:
        json.dump(remap_keys(refseq_dic), outfile)
        outfile.close()
    return refseq_dic

# def create_10length_dic(refseq_file , dpy13_blastn_wordsize7_xml):
#     """
#     tmp function to check myself (can be deleted)
#     :param refseq_file:
#     :param dpy13_blastn_wordsize7_xml:
#     :return:
#     """
#     refseq_dic=parse_refseq_to_dic(refseq_file)
#     result_handle=open(dpy13_blastn_wordsize7_xml)  # we can just open the saved file for input
#     blast_records=NCBIXML.parse(result_handle)
#     j=0
#     for blast_record in blast_records:
#         j+=1;
#         if (j > 30):
#             break;
#         # f.write('#####Query= ' + blast_record.query + "\n")
#         n=len(blast_record.alignments)
#         for i in range(n):
#             # if (len(blast_record.alignments[i].hsps) > 1):
#             # print('$$$$$$$$$$$$$$$$$$$ hsp include more than 1 segment pair $$$$$$$$$$$$$$$$$$$' + "\n")
#             alignment = blast_record.alignments[i]
#             if alignment.hit_def == "CHROMOSOME_III":
#                 add_blastRecord_to_refseq_dic(refseq_dic , alignment , blast_record.query)
#     print("be ready")
#     with open('10length_dic.txt' , 'w') as outfile:
#         json.dump(remap_keys(refseq_dic) , outfile)
#         outfile.close()
#     return refseq_dic


############################# BLAST 6/9 #############################
# dpy13_blastn_wordsize7_xml= "dpy13_blastn_wordsize7.xml"
# parse_refseq("WS220_refSeq.txt")
# create_full_dic("WS220_refSeq.txt", "dpy13_blastn_wordsize7.xml")
# processed_dic("full_dic.txt")
# create_10length_dic("WS220_refSeq.txt", "dpy13_blastn_wordsize7.xml")
# processed_dic("10length_dic.txt")
# analyze_blast_xml_results(dpy13_blastn_wordsize7_xml , "total_blast_results.fasta") #save blast results to "total_blast_results.fasta"


# We need to  run the folowing :

############################ BLAST (for all dpy13 gene - in 25 nt queries)#############################
# input_multifasta_file="E:\google drive\year 3\BioInformatics project\Dpy-13_multiFasta.fasta" #point to multifasta file that include the DPY-13+UTR (MY FILE) !
# output_multifasta_xml= "multifasta_blast_results_xml.xml"
# output_multifasta_xml = blast_to_multifasta_input(input_multifasta_file,output_multifasta_xml)
# analyze_blast_xml_results(output_multifasta_xml , "total_blast_results.fasta") #save blast results to "total_blast_results.fasta"
