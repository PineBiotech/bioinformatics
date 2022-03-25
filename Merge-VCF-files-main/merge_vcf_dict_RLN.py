#!/usr/bin/env python
import argparse
import gzip
import re
import string

des = """
A python tool to find intersection of VCF calling from two different sources into a single VCF file, 
so that overlapping variants are present only once in the output file.
A New tag in the INFO field is added named ‘calledBy’, 
that will have name of the tool which called the variant.
Any common INFO and FORMAT tags annotated by both tools are renamed 
by prefixing the tool name in the name of the tag e.g., Freebayes_DP, VarScan_DP
Raghavendran Lakshminarayanan, 23 March 2022 
Usage:
python merge_vcf_dict_RLN.py -i1 diploid_samples.vcf -i2 ClinVar_tp53.vcf -o common_dip_clvar_tp53.vcf"""


parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i1', '--vcf1', type=str,help='Freebayes_raw.vcf file \t[None]')
parser.add_argument('-i2', '--vcf2', type=str,help='Varscan_raw.vcf file \t[None]')
parser.add_argument('-o', '--out', type=str,help='vcf output file \t[None]')
args = parser.parse_args()

if args.vcf1 is not None:
    freebayes_path = args.vcf1
else:
    print('vcf1 path not given!')
    print ('Usage: \n python merge_vcf_dict_RLN.py -i1 diploid_samples.vcf -i2 ClinVar_tp53.vcf -o common_dip_clvar_tp53.vcf')
    raise IOError
if args.vcf2 is not None:
    varscan_path = args.vcf2
else:
    print('vcf2 path not given!')
    print ('Usage: \n python merge_vcf_dict_RLN.py -i1 diploid_samples.vcf -i2 ClinVar_tp53.vcf -o common_dip_clvar_tp53.vcf')
    raise IOError
if args.out is not None:
    out_file = args.out
else:
    print('out_file not given!')
    print ('Usage: \n python merge_vcf_dict_RLN.py -i1 diploid_samples.vcf -i2 ClinVar_tp53.vcf -o common_dip_clvar_tp53.vcf')
    raise IOError


### Read vcf file, dictionary is defined with [chr, pos, ref, alt] as indices.

def read_vcf(path):
    header, row, raw, last = [], [], [], ""
    R = {}
    #R = OrderedDict()
    if not path.endswith('.gz'):
        with open(path, 'r') as f:
            raw = f.readlines()
            for line in raw:
                if line.startswith('#'):
                    header += [line.replace('\n', '').rsplit('\t')]
                else:
                    if last != line:
                        splitedlines = line.split()
                        result= splitedlines[4].split(",")
                        match = re.search(r',', splitedlines[4])
                        if match:
                            lst=[]
                            for i in range(len(result)):
                                lst.append(splitedlines[:])

                            for i in lst:
                                i[4] = result.pop(0)
                                line1 = ' '.join(i)
                                #data += [line.replace('\n', '').rsplit('\t')]
                                row = [line1.replace('\n', '').rsplit(' ')]
                                #print (row[0][5])
                                C = row[0][0].replace('chr','')
                                C = C.replace('X',str(23))
                                C = C.replace('Y',str(24))
                                C = C.replace('M',str(25))
                                #q = row[0][0].replace('chr',''),int(row[0][1]),''.join(row[0][3]),''.join(row[0][4])
                                q = int(C),int(row[0][1]),''.join(row[0][3]),''.join(row[0][4])
                                t = tuple(q)
                                R[t] = row                
                        else:
                            #data += [line.replace('\n', '').rsplit('\t')]
                            row = [line.replace('\n', '').rsplit('\t')]
                            #print (row[0][5])
                            C = row[0][0].replace('chr','')
                            C = C.replace('X',str(23))
                            C = C.replace('Y',str(24))
                            C = C.replace('M',str(25))
                            #q = row[0][0].replace('chr',''),int(row[0][1]),''.join(row[0][3]),''.join(row[0][4])
                            q = int(C),int(row[0][1]),''.join(row[0][3]),''.join(row[0][4])
                            t = tuple(q)
                        R[t] = row
                    last = line 
                        
    else:
        with gzip.GzipFile(path,'rb') as f:
            raw = f.readlines()
            for line in raw:
                if line.startswith('#'):
                    header += [line.replace('\n', '').rsplit('\t')]
                else:
                    if last != line:
                        splitedlines = line.split()
                        result= splitedlines[4].split(",")
                        match = re.search(r',', splitedlines[4])
                        # If-statement after search() tests if it succeeded
                        if match:
                            lst=[]
                            for i in range(len(result)):
                                lst.append(splitedlines[:])

                            for i in lst:
                                i[4] = result.pop(0)
                                line1 = ' '.join(i)
                                #data += [line.replace('\n', '').rsplit('\t')]
                                row = [line1.replace('\n', '').rsplit(' ')]
                                #print (row[0][5])
                                C = row[0][0].replace('chr','')
                                C = C.replace('X',str(23))
                                C = C.replace('Y',str(24))
                                C = C.replace('M',str(25))
                                #q = row[0][0].replace('chr',''),int(row[0][1]),''.join(row[0][3]),''.join(row[0][4])
                                q = int(C),int(row[0][1]),''.join(row[0][3]),''.join(row[0][4])
                                t = tuple(q)
                                R[t] = row                
                        else:
                            #data += [line.replace('\n', '').rsplit('\t')]
                            row = [line.replace('\n', '').rsplit('\t')]
                            #print (row[0][5])
                            C = row[0][0].replace('chr','')
                            C = C.replace('X',str(23))
                            C = C.replace('Y',str(24))
                            C = C.replace('M',str(25))
                            #q = row[0][0].replace('chr',''),int(row[0][1]),''.join(row[0][3]),''.join(row[0][4])
                            q = int(C),int(row[0][1]),''.join(row[0][3]),''.join(row[0][4])
                            t = tuple(q)
                            R[t] = row
                    last = line 
    return header,R

print('Reading %s'%freebayes_path)
h1,d1 = read_vcf(freebayes_path)
print('Reading %s'%varscan_path)
h2,d2 = read_vcf(varscan_path)

# Determine the type of information in the header line.

def getheaderline(header):
    file_format, source,headerTitlestring, additional, info, fformat = "", "", "", [], [], []
    for j in header:
        if j[0].startswith("##fileformat"): 
            file_format = j[0].split("=",1)
        elif j[0].startswith("##source"): 
            source = j[0].split("=",1)
        elif j[0].startswith("##INFO"):
            info.append(j[0])
        elif j[0].startswith("##FORMAT"):
            fformat.append(j[0])
        elif j[0].startswith("#CHROM"):
            headerTitlestring = j
        else: additional.append(j[0])
    if not bool(file_format):
        file_format = "##fileformat=VCFv4.1".split("=",1)
    if not bool(source):
        source = "##source=freeBayes v1.0.2-29-g41c1313".split("=",1)
    return file_format, source, headerTitlestring, additional, info, fformat


def headerinfo(if_list):
    HeaderTags = {}
    HeaderString = {}
     
    for line in if_list:
        tag = line.split("=",1)
        tagID = (tag[1].split("ID=",1))[1].split(",",1)

        tagNumber = (tagID[1].split("Number=",1))[1].split(",",1)
        tagType = (tagNumber[1].split("Type=",1))[1].split(",",1)
        try: tagDescription = ( ( (tagType[1].split("Description=\"",1))[1] ).split("\">") )[0]
        except IndexError: tagDescription = ""
        tagID = tagID[0]; tagNumber = tagNumber[0]; tagType = tagType[0]

        HeaderTags[tagID] = tagNumber, tagType, tagDescription
        HeaderString[tagID] = line
    
    return HeaderTags, HeaderString

### Separate header line into different sections
h1_file, h1_source, h1_title, h1_additional, h1_info, h1_fformat = getheaderline(h1)
h2_file, h2_source, h2_title, h2_additional, h2_info, h2_fformat = getheaderline(h2)

### INFO and FORMAT separated to tags, descriptions 
h1_info_tags, h1_info_string = headerinfo(h1_info)
h1_format_tags, h1_format_string = headerinfo(h1_fformat)
h2_info_tags, h2_info_string = headerinfo(h2_info)
h2_format_tags, h2_format_string = headerinfo(h2_fformat)

### Tool source name derived, to be used to annotate common INFO and FORMAT field names
h1_source_name = h1_source[1].split(" ",1)[0]
h2_source_name = h2_source[1].split(" ",1)[0]

## Merge header lines for the output file
print('Merging header information')
h3_file_format = sorted(list(set(h1_file + h2_file)))
h3_additional = sorted(list(set(h1_additional + h2_additional)),reverse=True)
h3_source = sorted(list(set(h1_source + h2_source)))
h3_title = h1_title

### Identify common INFO and FORMAT tags to rename them with their respective tool name as a prefix 
comm_tags = h1_info_tags.keys() & h2_info_tags.keys()
for i in comm_tags:
    if i:
        if bool(h1_info_string.get('i')):
            h1_info_string[i] = h1_info_string[i].replace("##INFO=<ID=", ("##INFO=<ID="+h1_source_name+"_"))
        if bool(h2_info_string.get('i')):
            h2_info_string[i] = h2_info_string[i].replace("##INFO=<ID=", ("##INFO=<ID="+h2_source_name+"_"))

comm_tags = h1_format_tags.keys() & h2_format_tags.keys()
for i in comm_tags:
    if i:
        if bool(h1_format_string.get('i')):
            h1_format_string[i] = h1_format_string[i].replace("##FORMAT=<ID=", ("##FORMAT=<ID="+h1_source_name+"_"))
        if bool(h1_format_string.get('i')):
            h2_format_string[i] = h2_format_string[i].replace("##FORMAT=<ID=", ("##FORMAT=<ID="+h2_source_name+"_"))

### Use set operations on dictinoary keys to merge the vcf dictionaries, by finding common keys 
### concatenate them with unique keys from both vcf dictionaries
print('Merging variant calls')
d2_keys_not_in_d1 = d2.keys() - d1.keys()
d1_keys_not_in_d2 = d1.keys() - d2.keys()
common_keys = d1.keys() & d2.keys()

d = {}

### The INFO field in the common variants are edited to reflect they are called_by both the tools. 
### The common FORMAT fields are edited to differentiate between their source tool 
for i in common_keys:
    #d[i]=d1[i]+d2[i]
    #print(i)
    d[i]=d1[i]
    d[i][0][6] = d2[i][0][6]
#    d[i][0][7] = d1[i][0][7] + ";calledBy=Freebayes+VarScan;" + d2[i][0][7]
    d[i][0][7] = d1[i][0][7] + ";calledBy=" + h1_source_name + '+' + h2_source_name + ";" + d2[i][0][7]
    if len(d1[i][0])>8:
        split_formd1 = (d1[i][0][8].split(':'))
    else:
        split_formd1 = (d1[i][0][7].split(':'))
    if len(d2[i][0])>8:
        split_formd2 = (d2[i][0][8].split(':'))
    else:
        split_formd2 = (d2[i][0][7].split(':'))
    for j in split_formd1:
        if (j in split_formd2):
            idx = split_formd1.index(j)
            idy = split_formd2.index(j)
#            split_formd1[idx] = "_".join(["Freebayes", split_formd1[idx]])   
#            split_formd2[idy] = "_".join(["Varscan", split_formd2[idy]])   
            split_formd1[idx] = "_".join([h1_source_name, split_formd1[idx]])   
            split_formd2[idy] = "_".join([h2_source_name, split_formd2[idy]])   
    d[i][0][8] = ":".join(split_formd1 + split_formd2)
    if len(d1[i][0])>9:
        if len(d2[i][0])>9:
            d[i][0][9] = d1[i][0][9] + ":" + d2[i][0][9]
        else:
            d[i][0][9] = d1[i][0][9]
    else:
        if len(d2[i][0])>9:
            d[i][0][9] = d2[i][0][9]

### Unique variants in the respective dictionaries, not used as the intended use is to
### get intersection between two vcf files  
#for i in d1_keys_not_in_d2:
#    d[i]=d1[i]
#    d[i][0][7] = d1[i][0][7] + ";calledBy=" + h1_source_name + ";"

#for i in d2_keys_not_in_d1:
#    d[i]=d2[i]
#    d[i][0][7] = d2[i][0][7] + ";calledBy=" + h2_source_name + ";"

### Print out the merged VCF file

### Print out ##fileformat
with open(out_file, 'w') as f:
    for l in range(len(h3_file_format)):
        print (h3_file_format[l],end ="", file=f)
        if l == 0:
            print ("=",end ="", file=f)
        elif l != len(h3_file_format)-1 :
            print (",",end ="", file=f)
    print ("", file=f)        

### Print out ##source
with open(out_file, 'a') as f:
    for l in range(len(h3_source)):
        print (h3_source[l],end ="", file=f)
        if l == 0:
            print ("=",end ="", file=f)
        elif l != len(h3_source)-1 :
            print (",",end ="", file=f)
    print ("", file=f)        

### Print out phasing, filedate, contigs, commandline, filter etc    
    for l in range(len(h3_additional)):
        print (h3_additional[l], file=f)

### Print out INFO strings from file vcf1
    for l in h1_info_string:
        print (h1_info_string[l], file=f)

### Print out the new INFO field
    temp_str = '##INFO=<ID=calledBy,Number=A,Type=String,Description="Name of the tool which called the variant. If a variant is called by both tools, value is set as calledBy=' + h1_source_name + '+' + h2_source_name + '">'
    print (temp_str, file=f)
    
### Print out INFO strings from file vcf2
    for l in h2_info_string:
        print (h2_info_string[l], file=f)

### Print out Format strings from file vcf1
    for l in h1_format_string:
        print (h1_format_string[l], file=f)

### Print out Format strings from file vcf2
    for l in h2_format_string:
        print (h2_format_string[l], file=f)

### Print out Title
    for l in range(len(h3_title)):
        print (h3_title[l],"\t",end ="", file=f)
    print ("", file=f)
    
### Set operations on Dictionary messes up with the order, hence need to sort by chr, pos fields
### tried ordered dictionaries still the problem persists, hence the need for the sorting operation.
    com_key = sorted(list(d.keys()),key = lambda x: (int(x[0]), int(x[1])))

### Print out the merged variant calls 
    for value in com_key:
        for l in range(len(d[value][0])):
            print (d[value][0][l],"\t",end =" ", file=f)
        print ("", file=f)
f.close
print('Done !!! The output is in the file %s'%out_file)
