import sys
with open(sys.argv[1], "r") as f:
    psl_lines = f.readlines()

#FIGURE OUT HOW MANY LINES ARE HEADERS. SAY N
#psl_headers = psl_lines[:N]
#psl_lines = psl_lines[N:]

output = open(sys.argv[2], "w")
#output.writelines(psl_headers)

i = 0
psl_new=[]
while i < len(psl_lines):
    x = psl_lines[i]
    t = x.rstrip().split("\t")
    exon_num_x = t[-4]
    start_positions_x = t[-1]
    lengths_x = t[-3]
    start_position_list_x = start_positions_x.split(",")[:-1]
    length_list_x = lengths_x.split(",")[:-1]
    #merged_status=0
    j = 0
    while j < len(psl_lines):
        y = psl_lines[j]
        k = y.rstrip().split("\t")
        exon_num_y = k[-4]
        print("K[0]:" + str(k[0]) + str(j))
        start_positions_y = k[-1]
        lengths_y = k[-3]
        start_position_list_y = start_positions_y.split(",")[:-1]
        length_list_y = lengths_y.split(",")[:-1]
        if (exon_num_x == exon_num_y) and (x != y):   # no need to merge the same line with itself
            exon1_start_x = int(start_position_list_x[0])
            exon1_end_x = int(start_position_list_x[0]) +  int(length_list_x[0]) - 1
            exon1_start_y = int(start_position_list_y[0])
            exon1_end_y = int(start_position_list_y[0]) +  int(length_list_y[0]) - 1
            exon2_start_x = int(start_position_list_x[-1])
            exon2_end_x = int(start_position_list_x[-1]) +  int(length_list_x[-1]) - 1
            exon2_start_y = int(start_position_list_y[-1])
            exon2_end_y = int(start_position_list_y[-1]) +  int(length_list_y[-1]) - 1
            print("Testing new")
            print("Exon Number:" + str(exon_num_x))
            print("start_position_list_x:" + str(start_position_list_x[1:-1]) + "start position list y:" + str(start_position_list_y[1:-1]))
            print("length list x:" + str(length_list_x[1:-1]) + "length list y:" + str(length_list_y[1:-1]))
            print("exon1_end_x:" + str(exon1_end_x) + " exon1_end_y:" + str(exon1_end_y))
            print("exon2_start_x:" + str(exon2_start_x) + " exon2_start_y:" + str(exon2_start_y))
            if int(exon_num_x) == 2:
                if (exon1_end_x == exon1_end_y) and (exon2_start_x == exon2_start_y):
                    #MERGE the two lines
                    if int(start_position_list_x[0]) >= int(start_position_list_y[0]):
                      start_position_list_x[0] = start_position_list_y[0]
                      length_list_x[0] = length_list_y[0]
                    if int(start_position_list_y[-1]) >= int(start_position_list_x[-1]):
                      start_position_list_x[-1] = start_position_list_y[-1]
                      length_list_x[-1] = length_list_y[-1]
                    
                    merged_status=1
                    t[0]=t[0] + "|"+ k[0]
                    t[-1] = ",".join(start_position_list_x) + ","
                    t[16] = str(start_position_list_x[0])#first element of start position x
                    t[17] = str((int(start_position_list_x[-1]) + int(length_list_x[-1]))-1)
                    t[15] = str((int(t[17])-int(t[16]))+ 1)
                    t[-3] = ",".join(length_list_x) + ","
                    new_line = "\t".join(t) + "\n"
                    #psl_lines[i] = new_line
                    psl_lines.remove(psl_lines[j])
                    #psl_lines.remove(y)   # REMOVE y from the looping list
                    #merged_status=1
                    t[0]=t[0] + "|"+ k[0]
                    j += 0
                else:
                    j += 1    
            if int(exon_num_x) >= 3:
                if (start_position_list_x[1:-1] == start_position_list_y[1:-1]) and (length_list_x[1:-1] == length_list_y[1:-1]) and (exon1_end_x == exon1_end_y) and (exon2_start_x == exon2_start_y):
                    #MERGE the two lines
                    print("I printed X.")
                    print(x)
                    #print(len(psl_lines))
                    if int(start_position_list_x[0]) >= int(start_position_list_y[0]):
                      start_position_list_x[0] = start_position_list_y[0]
                      length_list_x[0] = length_list_y[0]
                    if int(start_position_list_y[-1]) >= int(start_position_list_x[-1]):
                      start_position_list_x[-1] = start_position_list_y[-1]
                      length_list_x[-1] = length_list_y[-1]
                    #merged_status=1
                    t[0]=t[0] + "|"+ k[0]
                    t[-1] = ",".join(start_position_list_x) + ","
                    t[16] = str(start_position_list_x[0])#first element of start position x
                    t[17] = str((int(start_position_list_x[-1]) + int(length_list_x[-1]))-1)
                    t[15] = str((int(t[17])-int(t[16]))+ 1)
                    t[-3] = ",".join(length_list_x) + ","
                    new_line = "\t".join(t) + "\n"
                    #psl_lines[j]=new_line
                    #print("Printing I:" + str(i))
                    #psl_lines.remove(psl_lines[i])
                    #psl_lines.remove(psl_lines[j-1])
                    #psl_lines[i]=new_line
                    psl_lines.remove(psl_lines[j])
                    print("PSL LINES length: " + str(len(psl_lines)))
                    print(new_line)
                    j += 0
                else:
                    j += 1
        else:
            j += 1
         
    #NOW we have loop through the whole list for the current line x, we can update x
    t[-1] = ",".join(start_position_list_x) + ","
    t[16] = str(start_position_list_x[0])#first element of start position x
    t[17] = str((int(start_position_list_x[-1]) + int(length_list_x[-1]))-1)
    t[15] = str((int(t[17])-int(t[16]))+ 1)
    t[-3] = ",".join(length_list_x) + ","
    print(t[14:22])
    new_line = "\t".join(t) + "\n"
    output.write(new_line)
    i+=1
output.close()