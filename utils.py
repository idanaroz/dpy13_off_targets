import json
import os
import collections
import statistics

MAX_OFF_OFF = "62"
T="\t"; TT=2 * "\t";TTT=3 * "\t";TTTT = 4*"\t";

from math import floor
def floored_percentage(val, digits):
    val *= 10 ** (digits + 2)
    return '{1:.{0}f}%'.format(digits, floor(val) / 10 ** digits)


def extract_summary_from_report(report):
    """input report file
       extract the last 3 lines of summary from there
       return tuple of : (D>0 genes, off_off  , extra genes) specifically
    """
    positive_delta_genes =0;
    off_off=0
    extra_genes=0
    with open(report,'r') as f:
        for line in f:
            if line[0]=='*':
                positive_delta_genes = line.split()[1]
                line = f.readline(); assert(line[:2]=='**')
                off_off = line.split()[1]
                line = f.readline(); assert(line[:3]=='***')
                extra_genes = line.split()[1]
    f.close()
    return (positive_delta_genes, off_off, extra_genes)


def count_max_exl_genes4delta(delta):
    """need execel_genes file"""
    cnt=0
    with open ("excel_genes.txt") as f:
        for line in f:
            b = int(line.split()[1])
            c = int(line.split()[2])
            if b-c >= delta:
                cnt+=1
    f.close()
    return cnt

def count_max_exl_off_off(delta):
    """ need execel_genes file """
    cnt=0
    with open ("excel_genes.txt") as f:
        for line in f:
            b = int(line.split()[1])
            c = int(line.split()[2])
            if b-c <= 0:
                cnt+=1
    f.close()
    return cnt



def parameters_overview_report():
    """
    create overview report from midgam files. for Further understanding go to :'midgam\summary.txt'
    """
    d={}
    max_delta_dic={}
    max_off_off_dic={}
    with open("midgam\summary.txt" , 'w') as summary_f:
        for fn in os.listdir('midgam'):
            if fn == "summary.txt":
                continue
            f_locat="midgam\\" + fn
            tup=extract_summary_from_report(f_locat)
            delta = int (fn.replace(',',' ').replace('.',' ').split()[-2])
            #find max delta
            if delta not in max_delta_dic:
                max_d_genes = count_max_exl_genes4delta(delta)
                max_delta_dic[delta] = max_d_genes
            else:
                max_d_genes  = max_delta_dic[delta]
            D_genes = str(tup[0])
            ratio_d_genes =  floored_percentage((int(D_genes) / max_d_genes),2) #precentage
            off_off = str(tup[1])
            max_off_off = MAX_OFF_OFF #constant '62'
            ratio_off_off =floored_percentage(int(off_off) / int(max_off_off),2)
            extra_genes = str(tup[2])
            line =  "\t" + D_genes + TTT + str(max_d_genes)+ TTT +  str(ratio_d_genes) + TTT + off_off + TTT + max_off_off + TTT \
                   +str(ratio_off_off) + TTT + extra_genes + 6 * "\t" + fn + "\n"
            if delta in d:
                d[delta] += [line]
            else:
                d[delta] = [line]
        #sort d by key from down to up
        od =  collections.OrderedDict(sorted(d.items()))
        for key in od:
            key_s = str(key)
            s=''
            if len(key_s)==1:
                s = '  '
            summary_f.write('---------------------------------------------------------------------------------------------------------------------------\n')
            summary_f.write("D>="+ key_s+ "genes"+ s + T +"maxD>=" + key_s+ 2* "\t" + "ratioD" + TT  + "  off_off" +T +" max_ofof"\
            +"\t"+ "ratio_off_off" +T+ "Extra_genes\t\tGrade?\t\tparams(score,e-val,length,identity,matches,D)\n")
            summary_f.write('---------------------------------------------------------------------------------------------------------------------------\n')
            value_list = od[key]
            for line in value_list:
                summary_f.write(line)
    summary_f.close()
    return


def coordinate_avg(coord_str):
    """
    :param coord_str:
    helper function for graph analysis
    """
    arr = coord_str.split("_")
    return str((int(arr[0])+int(arr[1]))//2)
def prepare_to_graph(parsed_by_DPY_13_file):
    global d
    d={}
    sorted_l=[]
    with open(parsed_by_DPY_13_file , 'r') as f:
        for line in f:
            if line[:6] == "Dpy-13":
                key=line[7:]  # 964_988
            if line[:5] == "#hits":
                val=line[6:]  # 296
                assert (key not in d)
                d[key]=val
        f.close()
        for k in sorted(d , key=lambda k: int(d[k]) , reverse=True):  # sort by value (higher first)
            sorted_l.append([k , d[k]])
        with open("sorted.txt" , 'w') as sorted_f:
            for t in sorted_l:
                sorted_f.write(coordinate_avg(t[0].replace("\n" , "")) + "\t" + t[1])
        sorted_f.close()
    return




######################################################################################################################
### Mission1:For VISUALIZATION graphs from midgam_summary
######################################################################################################################

def extract_params_from_midgam_summary_file(line , requsted_param):
    """line order [D>=7genes , maxD>=7 , ratioD , off_off , max_ofof , ratio_off_off , Extra_genes , params(score,e-val,length,identity,matches,D)]"""
    arr = line.split()
    assert(len(arr)==8)# make sure it's params line
    if requsted_param == "ratioD":
        return arr[2]
    elif requsted_param == "ratio_off_off" :
        return arr[5]
    params = arr[-1] # 14,0.001,14,100,16,7.txt
    params = params.replace(".txt","").replace(',', ' ').split()#[score,e-val,length,identity,matches,D]
    param_dic = {"score":0,"e-val":1,"length":2,"identity":3,"matches":4,"D":5}
    if requsted_param in param_dic:
        return params[param_dic[requsted_param]]
    else:
        print("you inserted Not valid requsted_param: " + requsted_param)
        raise exec()

# score_limit_list = ['14', '16', '18', '20']
def build_SCORE_red_gray_graph_data():
    """The output of this function is in :
    '/mission1 - paramater examination (to each parm) - data'  folder  """
    summary_file = "midgam\midgam_summary.txt"
    sum_f = open(summary_file, 'r')
    score_limit_list = ['14', '16', '18', '20']
    score_limit_dic = {'14':0, '16':1, '18':2, '20':3}
    score_file_d_list = []
    for score in score_limit_list:
        score_file_d =  open("mission1\score\score" +score + ".txt", 'w'); #"mission1\score\\14.txt"
        score_file_d.write("ratioD\tratio_off_off\n")
        score_file_d_list.append(score_file_d)
    for line in sum_f:
        if line[0] == "-" or line[0] == "D" :
            continue
        ratioD = extract_params_from_midgam_summary_file(line, "ratioD")[:-1]
        ratio_off_off = extract_params_from_midgam_summary_file(line, "ratio_off_off")[:-1]
        score = extract_params_from_midgam_summary_file(line, "score")
        line_to_write = ratioD + "\t" + ratio_off_off + "\n"
        assert (score=="14" or score=="16" or score =="18" or score == "20")
        index = score_limit_dic[score]
        score_file_d_list[index].write(line_to_write)
    for descriptor_file in score_file_d_list:
        descriptor_file.close()
    sum_f.close()

# expect_upper_limit_set={1e-03 , 1e-01 , 3 , 4}
def build_eVAL_red_gray_graph_data():
    summary_file = "midgam\midgam_summary.txt"
    sum_f = open(summary_file, 'r')
    e_val_limit_list = ['0.001' , '0.1' , '3', '4']
    e_val_limit_dic = {'0.001':0, '0.1':1, '3':2, '4':3}
    e_val_file_d_list = []
    for e_val in e_val_limit_list:
        score_file_d =  open("mission1\e-val\eval" +e_val + ".txt", 'w'); #"mission1\score\\14.txt"
        score_file_d.write("ratioD\tratio_off_off\n")
        e_val_file_d_list.append(score_file_d)

    for line in sum_f:
        if line[0] == "-" or line[0] == "D" :
            continue
        ratioD = extract_params_from_midgam_summary_file(line, "ratioD")[:-1]
        ratio_off_off = extract_params_from_midgam_summary_file(line, "ratio_off_off")[:-1]
        eVal = extract_params_from_midgam_summary_file(line, "e-val")
        line_to_write = ratioD + "\t" + ratio_off_off + "\n"
        assert (eVal in e_val_limit_list)
        index = e_val_limit_dic[eVal]
        e_val_file_d_list[index].write(line_to_write)
    for descriptor_file in e_val_file_d_list:
        descriptor_file.close()
    sum_f.close()

# alignment_min_length_limit={14,16,18,20}
def build_LENGTH_red_gray_graph_data():
    """The output of this function is in :
    '/mission1 - paramater examination (to each parm) - data'  folder  """
    summary_file = "midgam\midgam_summary.txt"
    sum_f = open(summary_file, 'r')
    length_limit_list = ['14','16','18','20']
    length_limit_dic = {'14':0, '16':1, '18':2, '20':3}
    length_file_d_list = []
    for len in length_limit_list:
        length_file_d =  open("mission1\length\length" +len + ".txt", 'w'); #"mission1\score\\14.txt"
        length_file_d.write("ratioD\tratio_off_off\n")
        length_file_d_list.append(length_file_d)
    for line in sum_f:
        if line[0] == "-" or line[0] == "D" :
            continue
        ratioD = extract_params_from_midgam_summary_file(line, "ratioD")[:-1]
        ratio_off_off = extract_params_from_midgam_summary_file(line, "ratio_off_off")[:-1]
        length = extract_params_from_midgam_summary_file(line, "length")
        line_to_write = ratioD + "\t" + ratio_off_off + "\n"
        assert (length in length_limit_list)
        index = length_limit_dic[length]
        length_file_d_list[index].write(line_to_write)
    for descriptor_file in length_file_d_list:
        descriptor_file.close()
    sum_f.close()

# min_identity={60,80,100}
def build_IDENTITY_red_gray_graph_data():
    """The output of this function is in :
    '/mission1 - paramater examination (to each parm) - data'  folder  """
    summary_file = "midgam\midgam_summary.txt"
    sum_f = open(summary_file, 'r')
    identity_limit_list = ['60','80','100']
    identity_limit_dic = {'60':0, '80':1, '100':2}
    identity_file_d_list = []
    for len in identity_limit_list:
        id_file_d =  open("mission1\identity\identity" +len + ".txt", 'w'); #"mission1\score\\14.txt"
        id_file_d.write("ratioD\tratio_off_off\n")
        identity_file_d_list.append(id_file_d)
    for line in sum_f:
        if line[0] == "-" or line[0] == "D" :
            continue
        ratioD = extract_params_from_midgam_summary_file(line, "ratioD")[:-1]
        ratio_off_off = extract_params_from_midgam_summary_file(line, "ratio_off_off")[:-1]
        identity = extract_params_from_midgam_summary_file(line, "identity")
        line_to_write = ratioD + "\t" + ratio_off_off + "\n"
        assert (identity in identity_limit_list)
        index = identity_limit_dic[identity]
        identity_file_d_list[index].write(line_to_write)
    for descriptor_file in identity_file_d_list:
        descriptor_file.close()
    sum_f.close()

# min_matches ={16,18,20,22,25}
def build_MATCHES_red_gray_graph_data():
    """The output of this function is in :
    '/mission1 - paramater examination (to each parm) - data'  folder  """
    summary_file = "midgam\midgam_summary.txt"
    sum_f = open(summary_file, 'r')
    matches_limit_list = ['16','18','20','22','25']
    matches_limit_dic = {'16':0, '18':1, '20':2 ,'22':3,'25':4}
    matches_file_d_list = []
    for m in matches_limit_list:
        id_file_d =  open("mission1\matches\match" + m + ".txt", 'w'); #"mission1\score\\14.txt"
        id_file_d.write("ratioD\tratio_off_off\n")
        matches_file_d_list.append(id_file_d)
    for line in sum_f:
        if line[0] == "-" or line[0] == "D" :
            continue
        ratioD = extract_params_from_midgam_summary_file(line, "ratioD")[:-1]
        ratio_off_off = extract_params_from_midgam_summary_file(line, "ratio_off_off")[:-1]
        matches = extract_params_from_midgam_summary_file(line, "matches")
        line_to_write = ratioD + "\t" + ratio_off_off + "\n"
        assert (matches in matches_limit_list)
        index = matches_limit_dic[matches]
        matches_file_d_list[index].write(line_to_write)
    for descriptor_file in matches_file_d_list:
        descriptor_file.close()
    sum_f.close()


# m_delta = {7,19,59,99}
def build_DELTA_red_gray_graph_data():
    """The output of this function is in :
    '/mission1 - paramater examination (to each parm) - data'  folder  """
    summary_file = "midgam\midgam_summary.txt"
    sum_f = open(summary_file, 'r')
    delta_limit_list = ['7','19','59','99']
    delta_limit_dic = {'7':0, '19':1, '59':2 ,'99':3}
    delta_file_d_list = []
    for D in delta_limit_list:
        id_file_d =  open("mission1\D\delta" + D + ".txt", 'w'); #"mission1\score\\14.txt"
        id_file_d.write("ratioD\tratio_off_off\n")
        delta_file_d_list.append(id_file_d)
    for line in sum_f:
        if line[0] == "-" or line[0] == "D" :
            continue
        ratioD = extract_params_from_midgam_summary_file(line, "ratioD")[:-1]
        ratio_off_off = extract_params_from_midgam_summary_file(line, "ratio_off_off")[:-1]
        delta = extract_params_from_midgam_summary_file(line, "D")
        line_to_write = ratioD + "\t" + ratio_off_off + "\n"
        assert (delta in delta_limit_list)
        index = delta_limit_dic[delta]
        delta_file_d_list[index].write(line_to_write)
    for descriptor_file in delta_file_d_list:
        descriptor_file.close()
    sum_f.close()


######################################################################################################################
### Mission 2:For VISUALIZATION graphs from midgam_summary
######################################################################################################################


def missio2_median_combination():
    """
    output of this function is in: 'mission2 - medians visualization DATA' folder
    :return:
    """
    length_set = {'14', '16', '18', '20'}
    identity_set = {'60', '80', '100'}
    matches_set = {'16', '18', '20', '22', '25'}
    combinations_dic = {}
    for l in length_set:
        for id in identity_set:
            for m in matches_set:
                combinations_dic[(l, id, m)] = {"ratioD": [], "off_off_ratios": []}
    summary_file = "midgam\midgam_summary.txt"
    sum_f = open(summary_file, 'r')
    for line in sum_f:
        if line[0] == "-" or line[0] == "D":
            continue
        ratioD = extract_params_from_midgam_summary_file(line, "ratioD")
        ratio_off_off = extract_params_from_midgam_summary_file(line, "ratio_off_off")
        length = extract_params_from_midgam_summary_file(line, "length")
        identity = extract_params_from_midgam_summary_file(line, "identity")
        matches = extract_params_from_midgam_summary_file(line, "matches")
        tup = (length, identity, matches)
        assert (tup in combinations_dic)
        combinations_dic[tup]["ratioD"] += [float(ratioD.strip('%'))]
        combinations_dic[tup]["off_off_ratios"] += [float(ratio_off_off.strip("%"))]
    sum_f.close()
    with open("mission2\medians.txt", 'w') as f:
        f.write("(len,id,matches)\t ratioD_median \t off_off_ratios_median\n")
        for tup_key in combinations_dic:
            ratioD = combinations_dic[tup_key]["ratioD"]
            off_off_ratios = combinations_dic[tup_key]["off_off_ratios"]
            ratioD_median = "{0:.2f}".format(statistics.median(ratioD))
            ratio_off_off_median = "{0:.2f}".format(statistics.median(off_off_ratios))
            params = tup_key[0]+"," +tup_key[1] + ","+ tup_key[2]
            f.write(params + T + ratioD_median + T + ratio_off_off_median + "\n")
    f.close()
    print("finish")

######################################################################################################################
### Mission 3 explore by "square" - to understand better go to the paper page 10 at the folder 'Paper+excell' (Computational biology project)
######################################################################################################################

def vectors_params_dictionary(lower_range_off,upper_range_off,lower_range_ratioD,upper_range_ratioD ):
    score_d = {'14':0, '16':1, '18':2, '20':3}
    eval_d = {'0.001':0 , '0.1':1 , '3':2 , '4':3}
    length_d = {'14':0,'16':1,'18':2,'20':3}
    id_d={'60':0,'80':1,'100':2}
    matches_d = {'16':0, '18':1, '20':2, '22':3, '25':4}
    D_d = {'7':0, '19':1, '59':2, '99':3}
    summary_file = "midgam\midgam_summary.txt"
    sum_f = open(summary_file, 'r')
    hits_cnt=0
    vectors_params_dic = {"score":[0,0,0,0],"e-val":[0,0,0,0] , "length":[0,0,0,0], \
                          "id":[0,0,0],"matches":[0,0,0,0,0],"D":[0,0,0,0] }
    for line in sum_f:
        if line[0] == "-" or line[0] == "D" :
            continue
        ratioD = float (extract_params_from_midgam_summary_file(line, "ratioD")[:-1])
        ratio_off_off = float (extract_params_from_midgam_summary_file(line, "ratio_off_off")[:-1])
        if  (lower_range_ratioD < ratioD) and(ratioD<=upper_range_ratioD) and \
        (lower_range_off < ratio_off_off) and(ratio_off_off<= upper_range_off):
            hits_cnt+=1

            file_score = extract_params_from_midgam_summary_file(line, "score")
            file_eval = extract_params_from_midgam_summary_file(line, "e-val")
            file_length = extract_params_from_midgam_summary_file(line, "length")
            file_identity = extract_params_from_midgam_summary_file(line, "identity")
            file_matches = extract_params_from_midgam_summary_file(line, "matches")
            file_D = extract_params_from_midgam_summary_file(line, "D")

            vectors_params_dic["score"][score_d[file_score]] +=1
            vectors_params_dic["e-val"][eval_d[file_eval]] +=1
            vectors_params_dic["length"][length_d[file_length]] +=1
            vectors_params_dic["id"][id_d[file_identity]] +=1
            vectors_params_dic["matches"][matches_d[file_matches]] +=1
            vectors_params_dic["D"][D_d[file_D]] +=1

    assert(sum(vectors_params_dic["score"]) == hits_cnt)
    assert(sum(vectors_params_dic["e-val"]) == hits_cnt)
    assert(sum(vectors_params_dic["length"]) == hits_cnt)
    assert(sum(vectors_params_dic["id"]) == hits_cnt)
    assert(sum(vectors_params_dic["matches"]) == hits_cnt)
    assert(sum(vectors_params_dic["D"]) == hits_cnt)

    sum_f.close()
    return str(hits_cnt), vectors_params_dic


def vectors_params_dictionary_with_constraints(lower_range_off, upper_range_off, lower_range_ratioD, upper_range_ratioD, constraint_name,
                              constraint_values_dic):
    assert (
    constraint_name == "score" or constraint_name == "e-val" or constraint_name == "length" or constraint_name == "identity" or constraint_name == "matches" or constraint_name == "D")
    score_d = {'14': 0, '16': 1, '18': 2, '20': 3}
    eval_d = {'0.001': 0, '0.1': 1, '3': 2, '4': 3}
    length_d = {'14': 0, '16': 1, '18': 2, '20': 3}
    id_d = {'60': 0, '80': 1, '100': 2}
    matches_d = {'16': 0, '18': 1, '20': 2, '22': 3, '25': 4}
    D_d = {'7': 0, '19': 1, '59': 2, '99': 3}
    summary_file = "midgam\midgam_summary.txt"
    sum_f = open(summary_file, 'r')
    hits_cnt = 0
    vectors_params_dic = {"score": [0, 0, 0, 0], "e-val": [0, 0, 0, 0], "length": [0, 0, 0, 0], \
                          "id": [0, 0, 0], "matches": [0, 0, 0, 0, 0], "D": [0, 0, 0, 0]}
    for line in sum_f:
        if line[0] == "-" or line[0] == "D":
            continue
        ratioD = float(extract_params_from_midgam_summary_file(line, "ratioD")[:-1])
        ratio_off_off = float(extract_params_from_midgam_summary_file(line, "ratio_off_off")[:-1])
        if (lower_range_ratioD < ratioD) and (ratioD <= upper_range_ratioD) and \
                (lower_range_off < ratio_off_off) and (ratio_off_off <= upper_range_off):

            file_score = extract_params_from_midgam_summary_file(line, "score")
            file_eval = extract_params_from_midgam_summary_file(line, "e-val")
            file_length = extract_params_from_midgam_summary_file(line, "length")
            file_identity = extract_params_from_midgam_summary_file(line, "identity")
            file_matches = extract_params_from_midgam_summary_file(line, "matches")
            file_D = extract_params_from_midgam_summary_file(line, "D")

            if constraint_name == "score":
                if file_score in constraint_values_dic:
                    hits_cnt += 1
                    vectors_params_dic["score"][score_d[file_score]] += 1
                    vectors_params_dic["e-val"][eval_d[file_eval]] += 1
                    vectors_params_dic["length"][length_d[file_length]] += 1
                    vectors_params_dic["id"][id_d[file_identity]] += 1
                    vectors_params_dic["matches"][matches_d[file_matches]] += 1
                    vectors_params_dic["D"][D_d[file_D]] += 1

            elif constraint_name == "e-val":
                if file_eval in constraint_values_dic:
                    hits_cnt += 1
                    vectors_params_dic["score"][score_d[file_score]] += 1
                    vectors_params_dic["e-val"][eval_d[file_eval]] += 1
                    vectors_params_dic["length"][length_d[file_length]] += 1
                    vectors_params_dic["id"][id_d[file_identity]] += 1
                    vectors_params_dic["matches"][matches_d[file_matches]] += 1
                    vectors_params_dic["D"][D_d[file_D]] += 1

            elif constraint_name == "length":
                if file_length in constraint_values_dic:
                    hits_cnt += 1
                    vectors_params_dic["score"][score_d[file_score]] += 1
                    vectors_params_dic["e-val"][eval_d[file_eval]] += 1
                    vectors_params_dic["length"][length_d[file_length]] += 1
                    vectors_params_dic["id"][id_d[file_identity]] += 1
                    vectors_params_dic["matches"][matches_d[file_matches]] += 1
                    vectors_params_dic["D"][D_d[file_D]] += 1

            elif constraint_name=="identity":
                if file_identity in constraint_values_dic:
                    hits_cnt += 1
                    vectors_params_dic["score"][score_d[file_score]] += 1
                    vectors_params_dic["e-val"][eval_d[file_eval]] += 1
                    vectors_params_dic["length"][length_d[file_length]] += 1
                    vectors_params_dic["id"][id_d[file_identity]] += 1
                    vectors_params_dic["matches"][matches_d[file_matches]] += 1
                    vectors_params_dic["D"][D_d[file_D]] += 1

            elif constraint_name=="matches":
                if file_matches in constraint_values_dic:
                    hits_cnt += 1
                    vectors_params_dic["score"][score_d[file_score]] += 1
                    vectors_params_dic["e-val"][eval_d[file_eval]] += 1
                    vectors_params_dic["length"][length_d[file_length]] += 1
                    vectors_params_dic["id"][id_d[file_identity]] += 1
                    vectors_params_dic["matches"][matches_d[file_matches]] += 1
                    vectors_params_dic["D"][D_d[file_D]] += 1

            elif constraint_name=="D":
                if file_D in constraint_values_dic:
                    hits_cnt += 1
                    vectors_params_dic["score"][score_d[file_score]] += 1
                    vectors_params_dic["e-val"][eval_d[file_eval]] += 1
                    vectors_params_dic["length"][length_d[file_length]] += 1
                    vectors_params_dic["id"][id_d[file_identity]] += 1
                    vectors_params_dic["matches"][matches_d[file_matches]] += 1
                    vectors_params_dic["D"][D_d[file_D]] += 1
    # assert(sum(vectors_params_dic["score"]) == hits_cnt)
    # assert(sum(vectors_params_dic["e-val"]) == hits_cnt)
    # assert(sum(vectors_params_dic["length"]) == hits_cnt)
    # assert(sum(vectors_params_dic["id"]) == hits_cnt)
    # assert(sum(vectors_params_dic["matches"]) == hits_cnt)
    # assert(sum(vectors_params_dic["D"]) == hits_cnt)

    sum_f.close()
    return str(hits_cnt), vectors_params_dic


scores = {'14', '16', '18', '20'}
evals = {'0.001' , '0.1' , '3' , '4'}
lengths = {'14','16','18','20'}
identitys={'60','80','100'}
matches = {'16', '18', '20', '22', '25'}
deltas = {'7','19','59','99' }

# print(vectors_params_dictionary(50,60,90,100))
# print(vectors_params_dictionary_with_constraints(20,30,90,100 ,"D", {'59'}))
# print(vectors_params_dictionary_with_constraints(20,30,90,100 ,"D", {'99'}))





#######################################
####### cake
#######################################

def cakes(list_limits):
    """
    This function provide data for creating 'cake' visualization
    :param list_limits: len(list_limits)=6
    [0] = represents : ["score"]
    [1] = represents : ["e-val"]
    [2] = represents : ["length"]
    [3] = represents : ["id"]
    [4] = represents : ["matches"]
    [5] = represents : ["D"]

    :return:
    """


    table_cakes_file = "cake.txt"
    r_f = open(table_cakes_file, 'w')

    for limits in list_limits:
        lower_range_off, upper_range_off, lower_range_ratioD, upper_range_ratioD = limits[0],limits[1],limits[2],limits[3]
        (hits_cnt, vectors_params_dic) = vectors_params_dictionary(lower_range_off, upper_range_off, lower_range_ratioD, upper_range_ratioD)
        # [score ,e-val , length ,id , matches, d ]
        l = ['0','1','2','3','4','5']
        l[0] = str(vectors_params_dic["score"])
        l[1] = str(vectors_params_dic["e-val"])
        l[2] = str(vectors_params_dic["length"])
        l[3] = str(vectors_params_dic["id"])
        l[4] = str(vectors_params_dic["matches"])
        l[5] = str(vectors_params_dic["D"])
        r_f.write("off_off_ratio:"+str(lower_range_off) +"-" +str(upper_range_off)+"; ")
        r_f.write("ratioD: "+ str(lower_range_ratioD) +"-" + str(upper_range_ratioD)+"\n")
        r_f.write("(#hits: "+hits_cnt+", { score: " +l[0] +", ")
        r_f.write("e-val: " +l[1] +", ")
        r_f.write("length: " +l[2] +", ")
        r_f.write("id: " +l[3] +", ")
        r_f.write("matches: " +l[4] +", ")
        r_f.write("D: " +l[5] +"})\n")

    r_f.close()

# cakes([[20,30,90,100],[20,25,90,95],[25,30,90,95], [20,25,95,100] ,[25,30,95,100], \
#        [20, 30, 80, 90], [20,25,80,85], [25,30,80,85],[20,25,85,90],[25,30,85,90],\
#        [10,20,80,90],[10,15,80,85],[15,20,80,85],[10,15,85,90],[15,20,85,90],\
#        [10,20,70,80],[10,15,70,75],[15,20,70,75], [10,15,75,80],[15,20,75,80],\
#        [10,20,50,60],[10,15,50,55],[15,20,50,55], [10,15,55,60],[15,20,55,60],\
#        [0,10,60,70],[0,10,60,65] , [0,10,65,70],\
#        [0,10,50,60], [0,10,50,55], [0,10,55,60],\
#        [0,10,40,50], [0,10,40,45], [0,10,45,50],\
#        [40,50,90,100], [40,45,90,95],[45,50,90,95], [40,45,95,100],[45,50,95,100],\
#        [50,60,90,100],\
#        [60,70,90,100], [70,80,90,100]])

###################################
###last % match visualization
###################################

def build_density_function(match_value):
    """
    This function provide data for creating 'cake' visualization
    """
    matches_d = {'16': 0, '18': 1, '20': 2, '22': 3, '25': 4}
    summary_file = "midgam\midgam_summary.txt"
    sum_f = open(summary_file, 'r')
    off_off_vector=[0,0,0,0,0,0,0,0,0,0] #len=10
    ratioD_vector= [0,0,0,0,0,0,0,0,0,0]#len=10
    all_hits =0

    density_off_file = open("off ratio density "+match_value+".txt",'w')
    density_ratioD_file = open("ratioD density " +match_value+".txt",'w')
    density_off_file.write("0-10\t10-20\t20-30\t30-40\t40-50\t50-60\t60-70\t70-80\t80-90\t90-100\n")
    density_ratioD_file.write("0-10\t10-20\t20-30\t30-40\t40-50\t50-60\t60-70\t70-80\t80-90\t90-100\n")
    for line in sum_f:
        if line[0] == "-" or line[0] == "D":
            continue
        ratioD = float(extract_params_from_midgam_summary_file(line, "ratioD")[:-1])
        ratio_off_off = float(extract_params_from_midgam_summary_file(line, "ratio_off_off")[:-1])
        match_from_file = extract_params_from_midgam_summary_file(line, "matches")
        if(match_from_file==match_value):
            all_hits+=1
            if ratioD==100.0:
                ratioD_vector[9]+=1
            else:
                ratioD_vector[int(ratioD//10)] += 1
            off_off_vector[int(ratio_off_off//10)]+=1
    for i in range(10):
        ratioD_vector[i] /= (all_hits / 100) #percent
        density_ratioD_file.write(str(ratioD_vector[i])+"\t"  )

        off_off_vector[i] /= (all_hits / 100)
        density_off_file.write(str(off_off_vector[i])+"\t"  )

    density_off_file.close()
    density_ratioD_file.close()
    sum_f.close()
    print("ratioD_vector:" ,ratioD_vector)
    # print("off_off_vector:" ,off_off_vector)

# build_density_function("16")
# build_density_function("18")
# build_density_function("20")
# build_density_function("22")
# build_density_function("25")

