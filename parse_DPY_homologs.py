import re


def count_gaps(line):
    """pre : The line that is given include '-' only to symbol gaps
            return  the number of gaps in a sequence (by counting '-') """
    gaps_number=0
    for e in line:
        if e == "-":
            gaps_number+=1
    return gaps_number


def count_matches(line):
    """pre : The line that is given include '|' only to symbol matches
            return  the number of matches in a sequence (by counting '|') """
    match_number=0
    for e in line:
        if e == "|":
            match_number+=1
    return match_number


def extract_gene_from_key(key):
    return key.split()[6][5:].replace("\n" , "")


def transript_synonyms_file(WS220_refSeq , synonyms_file):
    """create file with the gene origins of the transcript"""
    with open(WS220_refSeq) as f:
        next(f)
        with open(synonyms_file , 'w') as o:
            for line in f:
                l=[s.strip() for s in line.split("\t")]
                o.write(l[1] + "\t" + l[12] + "\n")
    f.close()
    o.close()
    return


def transript_synonyms_dic(WS220_refSeq):
    """return dic with transcript as key and gene origin as value"""
    d={}
    with open(WS220_refSeq) as f:
        next(f)
        for line in f:
            l=[s.strip() for s in line.split("\t")]
            d[l[1]]=l[12]
    f.close()
    return d


# transript_synonyms_dic("WS220_refSeq.txt" )
# transript_synonyms_file("WS220_refSeq.txt" , "synonyms.txt")

##############################################################################################################
# PARSE BLAST RESULT INTO DIC
##############################################################################################################

def parse_by_transcript2dic(file , score_limit , expect_upper_limit , alignment_min_limit , min_identity , min_matches):
    """Parse BLAST result into dic
     e.g: {'NR_070178', '4', '+', '3674859', '4526957':set(hit1, hit2, hit3, hit4 ...) }
     """

    synonyms_dic=transript_synonyms_dic("WS220_refSeq.txt")
    dic={}
    # prepare dic with all genes as keys e.g  : 'NR_070178', '4', '+', '3674859', '4526957'
    with open(file , 'r') as f:
        for line in f:
            if line[0] == "$":  # New hit
                # Use set to guarantee uniqueness - each DPY-13 fragment can hit specific gene one time at most
                queries_set=set()
                key=line.strip("\n")
                transcript=key.split()[1][
                           1:-2]  # extract the transcript from whole key (<class 'list'>: ['$KEY:', "'NM_001025813',", "'1',", "'+',", "'4035894',", "'4078129'"])
                if (transcript in synonyms_dic):
                    key+="\t" + "gene:" + synonyms_dic[transcript]
            if line[0:18] == "query_name:Dpy-13:":
                query_name=line[11:].strip("\n")
            if line[0:5] == "Score":  # score line
                l=re.sub('[!@#$,]' , '' , line)
                arr=l.split()
                score=float(arr[1])
                expect=float(arr[5])
                gaps_num=0;
                matches=0;
                mismatches=0  # tmp
                alignment_len=int(arr[8])
                # count Query gaps
                line=f.readline();
                assert (line[:5] == "Query")
                gaps_num=count_gaps(line)
                # count matches
                line=f.readline();
                assert (line.count('|') > 2)
                matches=count_matches(line)
                # count Sbjct gaps
                line=f.readline();
                assert (line[:5] == "Sbjct")
                gaps_num+=count_gaps(line)
                mismatches=alignment_len - (gaps_num + matches)
                identity_precentage=(alignment_len - (gaps_num + mismatches)) / alignment_len * 100  # e.g 33.4

                # check that the required conditions are met and add
                if (score >= score_limit) and (expect <= expect_upper_limit) and (alignment_len >= alignment_min_limit)\
                        and (identity_precentage >= min_identity) and (matches>=min_matches):
                    queries_set.add(query_name)
            if line[0:10] == "----------":
                if (len(queries_set) > 0):
                    dic[key]=queries_set
    f.close()
    return dic


def convert2gene_dic(dic):
    """return dictionary with converted keys. e.g: '$KEY:	'NM_069126', '4', '+', '8824800', '8825959'	gene:col-3'
        into: 'col-3'
     """
    genes_dic={}
    for k in dic.keys():
        genes_dic[extract_gene_from_key(k)]=dic[k]
    return genes_dic


def parse_by_transcript2file_genes_concise(file , out_file_short , out_file_long , score_limit , expect_upper_limit ,
                                           alignment_min_limit , min_identity , min_matches):
    """E.g: $KEY:	'NR_070178', '4', '+', '3674859', '4526957'
            #hits:964                                          """

    dic=parse_by_transcript2dic(file , score_limit , expect_upper_limit , alignment_min_limit , min_identity , min_matches)
    gene_dic=convert2gene_dic(dic)
    # Write dic to files
    with open(out_file_short , 'w') as outfile_s:
        with open(out_file_long , 'w') as outfile_l:

            for k in sorted(gene_dic , key=lambda k: len(gene_dic[k]) , reverse=True):  # sort by value (higher first)
                val=gene_dic[k]
                # outfile_s.write("gene: " + k + "\n")
                # outfile_s.write("#hits: " + str(len(val)) + "\n")
                outfile_s.write(k +"\t"+ str(len(val))+ "\n")
                outfile_l.write("gene: " + k + "\n")
                outfile_l.write("#hits: " + str(len(val)) + "\n")
                for h in val:
                    outfile_l.write(h + "\n")
    outfile_s.close()
    outfile_l.close()
    return


def parse_by_transcript2files_basic(file , out_file_short , out_file_long , score_limit , expect_upper_limit ,
                                    alignment_min_limit , min_identity , min_matches):
    """E.g: $KEY:	'NR_070178', '4', '+', '3674859', '4526957'
            #hits:964                                          """

    dic=parse_by_transcript2dic(file , out_file_short , out_file_long , score_limit , expect_upper_limit ,
                                alignment_min_limit , min_identity , min_matches)
    # Write dic to files
    with open(out_file_short , 'w') as outfile_s:
        with open(out_file_long , 'w') as outfile_l:

            for k in sorted(dic , key=lambda k: len(dic[k]) , reverse=True):  # sort by value (higher first)
                val=dic[k]
                outfile_s.write(k + "\n")
                outfile_s.write("#hits: " + str(len(val)) + "\n")

                outfile_l.write(k + "\n")
                outfile_l.write("#hits: " + str(len(val)) + "\n")
                for h in val:
                    outfile_l.write(h + "\n")
    outfile_s.close()
    outfile_l.close()
    return


def make_report(full_result_f , score_limit , expect_upper_limit , alignment_min_limit , min_identity , min_matches , m_delta):
    """
    @pre: excell file name : excel_genes.txt
    :param full_result_f:
    :param score_limit:
    :param expect_upper_limit:
    :param alignment_min_limit:
    :param min_identity:
    :param min_matches:
    :param m_delta:
    """
    dic=parse_by_transcript2dic(full_result_f , score_limit , expect_upper_limit , alignment_min_limit , min_identity , min_matches)
    gene_dic=convert2gene_dic(dic)
    excel_genes_set=set()
    cnt_caught_genes=0
    count_off_off_targets=0  # genes that got alignment even thought it was not supposed by excel
    p="," ; q = ", "
    output_location="midgam/" + str(score_limit) + p + str(expect_upper_limit) + p + \
                    str(alignment_min_limit) + p + str(min_identity)+p + str(min_matches) \
                    +p+str(m_delta) + ".txt"
    with open("excel_genes.txt" , 'r') as excel_genes_f:
        with open(output_location , "w") as out_f:
            first_line="params: " + "score_limit=" + str(score_limit) + ", eVal_upper_limit=" + str(
                expect_upper_limit) + ", align_min_limit=" + str(alignment_min_limit) + ", min_identity="\
                       +str(min_identity)+", matches=" +str(min_matches)+", m_delta=" + str(m_delta)+"\n"
            second_line = "exl_b" + 2*"\t" + "exl_c" + "\t" + "exl_d" + "\t" + "hits" +"\n"
            out_f.write(first_line); out_f.write(second_line)
            for line in excel_genes_f:
                line=line.replace("\n" , "")
                arr=line.split()
                gene=arr[0]
                delta=int(arr[1]) - int(arr[2])
                excel_genes_set.add(gene)
                newline=line
                if gene in gene_dic:
                    hits4gene = len(gene_dic[gene])
                if gene in gene_dic:
                    newline+=("\t\t" + str(hits4gene) + "\n")
                    if delta >= m_delta:
                        cnt_caught_genes+=1
                    if delta <= 0: #=
                        count_off_off_targets+=1
                else:
                    newline+="\t\t0\n"
                out_f.write(newline)
            # compute score
            ###
            out_f.write("-----------------------------------------------------------\n")
            out_f.write("The other genes that we found\n")
            out_f.write("-----------------------------------------------------------\n")
            extra_genes_number=0
            for gene in gene_dic.keys():
                if gene not in excel_genes_set:
                    extra_genes_number+=1
                    out_f.write(gene + "\t" + str(len(gene_dic[gene])) + "\n")
            out_f.write("-----------------------------------------------------------\n")
            out_f.write("Summary:\n")
            out_f.write("-----------------------------------------------------------\n")
            out_f.write("* " + str(cnt_caught_genes) + " genes from excel that D>=" +str(m_delta) + "\n")
            out_f.write("** " + str(count_off_off_targets) + " Got off off targets (hits for a gene that his D<=0)\n")
            out_f.write("*** " + str(extra_genes_number) + " extra genes")
    excel_genes_f.close()
    out_f.close()
# make_report("FINAL RESULTS!.txt" , 10 , 3 , 10 , 60 , 10)
def make_all_table_reports():
    score_limit_set={14,16,18,20}
    expect_upper_limit_set={1e-03 , 1e-01 , 3 , 4}
    alignment_min_limit={14,16,18,20}
    min_identity={60,80,100}
    min_matches ={16,18,20,22,25}
    m_delta = {7,19,59,99}
    for a in score_limit_set:
        for b in expect_upper_limit_set:
            for c in alignment_min_limit:
                for d in min_identity:
                    for e in min_matches:
                        for delta in m_delta:
                            make_report("FINAL RESULTS!.txt" , a , b , c , d ,e , delta)
# make_report("FINAL RESULTS!.txt" , 10 , 4 , 10 , 10)
# make_all_table_reports()


# out_file0_short = '16.10 parse_by_transcript2file_genes_concise_short_ 14 , 4 , 14 , 80 , 17.txt'
# out_file0_long = '16.10 parse_by_transcript2file_genes_concise_long_ 14 , 4 , 14 , 80 , 17.txt'
# parse_by_transcript2file_genes_concise("FINAL RESULTS!.txt" , out_file0_short , out_file0_long , 14 , 4 , 14 , 80 , 17)

##############################################################################################################
###################################################################################################################

def parse_by_DPY_13(file , out_file_short , out_file_long , score_limit , expect_upper_limit , alignment_min_limit , min_identity,min_matches):
    """
    Dpy-13:95_119
       #hits:9
    :param file:
    :param out_file_short:
    :param out_file_long:
    :param score_limit:
    :param expect_upper_limit:
    :param alignment_min_limit:
    :param min_identity:
    :param min_matches:
    example outputs : 4.11 parsed_by_DPY_13_short score22.txt' , out_file1_long = '4.11 parsed_by_DPY_13_long score22.txt'
    """

    d={}
    startIndex=1;
    # prepare dic with all DPY-13 fragments as keys
    for i in range(startIndex , 965):  # in multifasta 964 is the last fragment
        end_index=i + 24
        # Use set to guarantee uniqueness -  each DPY-13 fragment can hit specific gene one time at most
        d["Dpy-13:" + str(i) + "_" + str(end_index)]=set()
    with open(file , 'r') as f:  #'final results final!' file to parse from
        for line in f:
            if line[0] == "$":  # New hit
                ref_seq_key=line
            if line[0:18] == "query_name:Dpy-13:":
                query_name=line[11:].strip("\n")
            if line[0:5] == "Score":  # score line
                l=re.sub('[!@#$,]' , '' , line)
                arr=l.split()
                score=float(arr[1])
                expect=float(arr[5])
                alignment_len=float(arr[8])

                # count Query gaps
                line = f.readline();
                assert (line[:5] == "Query")
                gaps_num = count_gaps(line)
                # count matches
                line = f.readline();
                assert (line.count('|') > 2)
                matches = count_matches(line)
                # count Sbjct gaps
                line = f.readline();
                assert (line[:5] == "Sbjct")
                gaps_num += count_gaps(line)
                mismatches = alignment_len - (gaps_num + matches)
                identity_precentage = (alignment_len - (gaps_num + mismatches)) / alignment_len * 100  # e.g 33.4

                # check that the required conditions are met and add
                if (score >= score_limit) and (expect <= expect_upper_limit) and (alignment_len >= alignment_min_limit) \
                        and (identity_precentage >= min_identity) and (matches >= min_matches):
                    k = ref_seq_key.replace('$KEY:\t', 'Trascript_details:')
                    d[query_name].add(k)
        f.close()
        #Write the dictionary to files
        with open(out_file_short , 'w') as outfile_short:  # watch # hits only
            with open(out_file_long , 'w') as outfile:  # Detailed
                for k in sorted(d , key=lambda k: len(d[k]) , reverse=True):  # sort by value (higher first)
                    val=d[k]
                    outfile_short.write(k + "\n")
                    outfile_short.write("#hits:" + str(len(val)) + "\n")
                    if len(val) > 0:
                        outfile.write(str(k) + "\n")
                        outfile.write("#hits:" + str(len(val)) + "\n")
                        for ref_key in val:
                            outfile.write(ref_key)
        outfile_short.close()
        outfile.close()
    return

out_file1_short = '4.11 parsed_by_DPY_13_short score22.txt'
out_file1_long = '4.11 parsed_by_DPY_13_long score22.txt'
parse_by_DPY_13("FINAL RESULTS!.txt" , out_file1_short , out_file1_long , 16 , 3 , 16 , 80 , 22)


##############################################################################################################
# calculate the identity table
##############################################################################################################
