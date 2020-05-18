#!/usr/bin/env python
import argparse, os, sys, operator

def file_check(the_file) : # it checks if the input file and the search files exist
    if os.path.isfile(the_file) == False:
        print('ERROR No such file: ' + the_file + ' check the file name or the path!' + '\n' + 'ABORT')
        sys.exit()

def directory_check(out_file) : # it checks if the output directory exists
    if os.path.realpath(out_file) == False :
        print('ERROR No such output directory: ' + out_file + ' check the directory path!' + '\n' + 'ABORT')
        sys.exit()

def temporary(the_file, out_file) : # it creates the a temporary file where all the input sequences are stored
    if the_file.endswith('.fasta') :
        file_format = 0
        print('The input file is in the fasta format')
    elif the_file.endswith('.fastq') :
        file_format = 1
        print('The input file is in the fastq format')
    else :
        print('ERROR this is not a fasta/fastq file. Check the file format or file name')
    file_name = os.path.splitext(os.path.basename(out_file))[0]
    tmp_directory = os.path.dirname(os.path.realpath(out_file))
    tmp_locus = os.path.join(tmp_directory, file_name + '.tmp') # it creates the output file in the directory where the sequencing file is positioned
    with open(the_file) as f :
        read_name_list = []
        general_count = 0
        tmp_file = open(tmp_locus, 'w')
        if file_format == 0 :
            idx = 0
            for line in f : # it stores the names of the input reads (fasta format)
                if idx == 0 :
                    read_name_list.append(line)
                if idx == 1 :
                    tmp_file.write(line)
                    idx = -1
                idx += 1
        if file_format == 1 : # it stores the names of the input reads (fastq format)
            idx = 1
            for line in f :
                idx += 1
                if idx == 2 :
                    read_name_list.append(line)
                if idx == 3 :
                    tmp_file.write(line)
                    idx = -1
                general_count += 1
        tmp_file.close()
        
        read_name_list = [i.replace('\n', '') for i in read_name_list]
        print(read_name_list)
        return tmp_locus, read_name_list, tmp_directory, file_name

def input_search(search_file) : # it loads the sequences of the probes and makes the reverse complement of them
    with open(search_file, 'r') as f :
        sites = 0
        idx = 0
        probes_names = []
        fw_probes = []
        rv_probes = []
        for line in f :
            if idx%2 == 0 and line.find('>',) != -1:
                probes_names.append(line)
            elif idx%2 != 0 :
                fw_probes.append(line)
                sites += 1
            else :
                print('ERROR The search file is not in FASTA format!' + '\n' + 'ABORT')
                sys.exit()
            idx += 1

        fw_probes = [i.replace('\n','') for i in fw_probes]
        probes_names = [i.replace('\n', '\t') for i in probes_names]
        probes_number = len(probes_names)

        
        for seq in fw_probes : # it creates the reverse complement of each probe
            replace1 = seq.replace("A","t")
            replace2 = replace1.replace("T","a")
            replace3 = replace2.replace("G","c")
            replace4 = replace3.replace("C","g")
            replace5 = replace4.upper()
            reverse = replace5[::-1]
            rv_probes.append(reverse)
    
    return probes_names, fw_probes, rv_probes, probes_number
    
def find_function(target_sequence, mismatch_allowed, in_probes, out, read_name, read_eval) : # identifies the position of each probe in the input sequence
    check_probe = 0
    out.write(read_name[read_eval]+ '\t')
    round_cycle = 0
    tmp_list = []
    for x in range(0,2) : # by default, all the probes, fw and rv are searched in the read
        total_index = 0
        round_cycle += 1
        tmp_sublist = []
        cycle = 0
        for probe in in_probes[check_probe] :
            cycle += 1 
            found_index = 0
            nucleotide1 = 0
            tmp_sublist_inner = []
            for y in range(len(target_sequence) - (len(probe))) : # the overall identity between the probe and the sequence is evaluated
                nucleotide2 = 0
                match_index = 0
                for z in range(len(probe)) :
                    if target_sequence[nucleotide1 + nucleotide2] == probe[nucleotide2] :
                        match_index += 1
                    nucleotide2 += 1
                nucleotide1 += 1
                if match_index >= (len(probe) - mismatch_allowed): # if the overall identity between the probe and the sequence reaches the minimum threshold, the position is registered
                    if check_probe == 1 :
                        reverse_nucleotide = len(target_sequence) - nucleotide1 - len(probe) # it reverses the matching position
                        tmp_sublist_inner.append(reverse_nucleotide)
                    else :
                        tmp_sublist_inner.append(nucleotide1)
                    found_index += 1
                    total_index += 1
            if found_index == 0 :
                tmp_sublist_inner.append(str(0))
            tmp_sublist.append(tmp_sublist_inner)
        tmp_list.append(tmp_sublist)
        if total_index > 0 : # if the forward probes are found in the sequence, the cycle is broken and the reverse probes are not searched             
            break
        else : # if the forward probes are not found in the sequence, the reverse probes are seached
            check_probe = 1
    if len(tmp_list) == 1 : # the matching positions are printed in a new file
        for element in tmp_list[0] :
            count = 0
            for x in element :
                count += 1
                out.write(str(x))
                if count < len(element) :
                    out.write(',')
            out.write('\t')
    else :
        for element in tmp_list[1] :
            count = 0
            for x in element :
                count += 1
                out.write(str(x))
                if count < len(element) :
                    out.write(',')
            out.write('\t')

    out.write('\n')   

def distance_finder(out, out_directory, out_distance_name, read_name, probes_number) : # it is able to calculate the distance between two consecutive probes 
    with open(out, 'r') as f :
        count = 0
        count_2 = -1
        distance_file = os.path.join(out_directory, out_distance_name + '_distance' + '.txt') #it creates the output file where the distances are annotated
        distance_file_write = open(distance_file, 'w')
        for line in f: #every sequence contained in the previous output file is evaluated
            order_list = []
            count += 1
            if count > 1 : #the first line in the previous output file is not evaluated because it is the header
                list1 = line.split('\t') 
                list1 = list1[1:-1]
                list2 = []
                list3 = []
                for element in list1 :
                    list2 = (element.split(','))
                    list3.append(list2) 
                probe_eval = 1
                for probe in list3 :
                    for position in probe :
                        single_pos_list = []
                        single_pos_list.append(int(position)) #it is necessary to declare that the position is a number for the sorting function
                        single_pos_list.append(int(probe_eval))
                        order_list.append(single_pos_list)
                    probe_eval += 1
            new_ordered_list = sorted(order_list, key=operator.itemgetter(0)) #the list (order_list) is ordered on the basis of the first element (position) of each sublist
            
            ref_point1 = 0
            break_idx = 0
            new_ordered_list_1 = []
            
            for x in range (0, len(new_ordered_list) - 1) : # it trims the annotated distances to the first complete repeat
                if new_ordered_list[ref_point1][1] == 1 and new_ordered_list[ref_point1 + 1][1] == 2 :
                    for x in range (0, len(new_ordered_list) - ref_point1) :
                        new_ordered_list_1.append(new_ordered_list[ref_point1])
                        ref_point1 += 1
                    break_idx += 1
                else :
                    ref_point1 += 1
                if break_idx == 1 :
                    break
            
            ref_point1 = 1
            break_idx = 0
            new_ordered_list_2 = []
            new_order_list_3 = []

            for x in range (0, len(new_ordered_list_1) - 1) : # it trims the annotated distances to the last ccomplete repeat
                if new_ordered_list_1[-(ref_point1)][1] == 1 and new_ordered_list[-(ref_point1) - 1][1] == probes_number :
                    reverse_idx = 0
                    for x in range (0, len(new_ordered_list_1) - ref_point1 + 1) :
                        new_ordered_list_2.append(new_ordered_list_1[-(ref_point1)])
                        ref_point1 += 1
                    new_ordered_list_2.reverse()
                    break_idx += 1
                else :
                    ref_point1 += 1
                if break_idx == 1 :
                    break

            ref_point1 = 0
            ref_point2 = 1
            distance_list = []
            for x in range (0, len(new_ordered_list_2) - 1) :
                if new_ordered_list_2[ref_point1][1] == (new_ordered_list_2[ref_point1 + 1][1] - 1) : # it calculates the distance between the two probes. The distance is calculated only if the succession of the probes is correct (e.g. probe1 preceeds always probe 2 etc)
                    if new_ordered_list_2[ref_point1][0] != 0 :
                        distance = new_ordered_list_2[ref_point1 + 1][0] - new_ordered_list_2[ref_point1][0]
                        distance_list.append(distance)
                    else :
                        distance = 0
                        distance_list.append(distance)
                elif new_ordered_list_2[ref_point1][1] == probes_number and (new_ordered_list_2[ref_point1][1]) - (probes_number -1) ==  (new_ordered_list_2[ref_point1 + 1][1]): #if the very last probe preceeds the first one, then the distance is calculated
                    if new_ordered_list_2[ref_point1][0] != 0 :
                        distance = new_ordered_list_2[ref_point1 + 1][0] - new_ordered_list_2[ref_point1][0]
                        distance_list.append(distance)
                    else :
                        distance = 0
                        distance_list.append(distance)
                else : #if the two cases above are not true, the distance '0' is assigned (e.g. if probe 1 preeceds probe 3 because probe 2 was not found, then a default distance equal to 0 is assigned to P1P3)
                    distance_list.append(int(0))
                ref_point1 += 1
                ref_point2 += 1
            ref_point1 = 0
            if count_2 >= 0 : # in this way it is possible to elude the first line of the output document which is the header
                distance_file_write.write(read_name[count_2] + '\t') # the name of the read evaluated is printed
                for x in range(len(distance_list)) :
                    if int(distance_list[ref_point1]) > 0 : # Only distances > 0 are printed in the distance document. In this way it is possible to get rid of the distances of non consecutive probes
                        distance_file_write.write('P' + str(new_ordered_list_2[ref_point1][1]) + 'P' + str(new_ordered_list_2[ref_point1 + 1][1]) + ' ' + str(distance_list[ref_point1]) + '\t')
                        ref_point1 += 1
                    else :
                        ref_point1 += 1
                distance_file_write.write('\n')
            count_2 += 1
    distance_file_write.close()
    return distance_file

def overlap_detective(in_file, out_directory, out_overlap_name, read_name, tolerance, min_ovl, max_ovl) : # it is able to find the overlapping bacs using the distances annotated in the distance file
    overlap_file = os.path.join(out_directory, out_overlap_name + '_overlap' + '.txt') #it creates the output file where the overlaps are annotated
    overlap_file_write = open(overlap_file, 'w')
    main_list = []
    nested_list_1 = []
    with open(in_file) as f :
        for line in f :
            probe_distance = line.split('\t')
            main_list.append(probe_distance)
        for element in main_list :
            del element[0]
            del element[-1]
            nested_list_2 = []
            for inner_element in element :
                inner_element = inner_element.split(' ')
                nested_list_2.append(inner_element)
            nested_list_1.append(nested_list_2)
    compare_reads = len(nested_list_1) - 1
    index = 0
    for x in range (0,2) : # the overlaps are evaluated by checking the begin and the end of each read
        seq_round = 0
        for seq in nested_list_1 :
            read_count = 1 + seq_round
            for x in range(0, compare_reads - seq_round) :
                try :
                    current_ovl = max_ovl
                    for z in range (0, max_ovl - min_ovl) :
                        score = 0
                        if index == 0 :
                            overlap_to_find = -(current_ovl)
                            overlap_to_find_begin = 0
                        elif index == 1 :
                            overlap_to_find = 0
                            overlap_to_find_begin = -(current_ovl)
                        if current_ovl <= len(nested_list_1[seq_round]) and current_ovl <= len(nested_list_1[read_count]):
                            for x in range (0, current_ovl) :
                                if seq[overlap_to_find][0] == nested_list_1[read_count][(overlap_to_find_begin)][0] : # a first check is performed looking at the right succession of the probes
                                    if int(seq[overlap_to_find][1]) >= ((int(nested_list_1[read_count][(overlap_to_find_begin)][1])) - tolerance) and int(seq[overlap_to_find][1]) <= (int((nested_list_1[read_count][(overlap_to_find_begin)][1])) + tolerance) : # if the distances annotated are compatible, then the overlap is called
                                        score += 1
                                overlap_to_find += 1
                                overlap_to_find_begin += 1
                        else :
                            print('TOO SHORT')
                        if score == current_ovl :
                            overlap_file_write.write(read_name[seq_round]+ ' and ' + read_name[read_count] + ' overlap. ' + str(current_ovl) + ' features (distances) overlap' + '\n')
                            break
                        current_ovl -= 1   
                except IndexError:
                    pass
                read_count += 1
            seq_round += 1
        overlap_file_write.write('End evaluation change!' + '\n' ) # it makes clear what ends are evaluted
        index += 1



parser = argparse.ArgumentParser(description='Prova')

parser.add_argument('-i', '--input', help = 'Input fastq file', type = str, required = True )
parser.add_argument('-s', '--search', help = 'Input search file in fasta format', type = str, required = True)
parser.add_argument('-m', '--mismatches', help = 'Input mismatches allowed. Default = 0', type = int, default = 0)
parser.add_argument('-o', '--output', help = 'Output file', type = str, required = True)
parser.add_argument('-t', '--overlap_tolerance', help = 'Distance tolerance allowed to identify an overlap. Default = 5', type = int, default = 5)
parser.add_argument('-min_n', '--min_overlaps', help = 'Minimum number of matching distances to identify an overlap', type = int, default = 3)
parser.add_argument('-max_n', '--max_overlaps', help = 'Max number of matching distances to identify an overlap', type = int, default = 3)

args = parser.parse_args()
input_file = args.input
output_file = args.output
search_input = args.search
mismatch_input = args.mismatches
ovl_tolerance = args.overlap_tolerance
min_nmb_overlaps = args.min_overlaps
max_nmb_overlaps = args.max_overlaps
#print(args.output)

file_to_check = []
file_to_check.append(input_file and args.search)

for x in file_to_check :
    file_check(x)

directory_check(output_file)

tmp_locus, read_name_list, tmp_directory, file_name = temporary(input_file, output_file)
probes_names, fw_probes, rv_probes, probes_number = input_search(search_input)

print(tmp_directory)
#print(probes_names)
#print(rv_probes)

tmp_file_read = open(tmp_locus, 'r')

reads = []
reads.append(fw_probes)
reads.append(rv_probes)
#print(reads)

output_file_write = open(output_file, 'w')
output_file_write.write('read name' + '\t')
for name in probes_names :
    output_file_write.write(name)
output_file_write.write('\n')

read_evaluated = -1

for line in tmp_file_read :
    read_evaluated += 1
    find_function(target_sequence = line, mismatch_allowed = mismatch_input, in_probes = reads, out = output_file_write, read_name = read_name_list, read_eval = read_evaluated)

output_file_write.close()

distance_file = distance_finder(out = output_file, out_directory = tmp_directory, out_distance_name = file_name, read_name = read_name_list, probes_number = probes_number)
print(distance_file)

print(read_name_list)
overlap_detective(in_file = distance_file, out_directory = tmp_directory, out_overlap_name = file_name, read_name = read_name_list, tolerance = ovl_tolerance, min_ovl = min_nmb_overlaps, max_ovl = max_nmb_overlaps)

os.remove(tmp_locus)

  
