#!/usr/bin/env python2.7
import sys, string
import argparse
from Bio import SeqIO
from collections import defaultdict

parser = argparse.ArgumentParser(description="""

Given a fasta/fastq file and a file of entry names, return from the fastx file only those entries (or only other entries).

Assumes each given name only occurs once in the file by default.

For very large fastx files, if you are looking to get >50% of the entries with supplied names,
it might be more efficient to supply the names of entries you do not want with "-e" option.

Similarly, for very large fastx files, if you are looking to exclude >50% of the entries with -e,
it might be more efficient to supply the names of entries you do want (without -e).

It more realistically depends on your workflow though.

"very large fastx files" are likely those with tens of millions of entries or more.
On a macbook pro, for a fasta file with 87509 scaffolds:
-- the approach of names to keep for 87256 names took 20 seconds
-- the approach of names to exclude for the complementary 253 entry names took 17 seconds
Not a big difference.

    """, formatter_class= argparse.RawTextHelpFormatter)


inputtype = parser.add_mutually_exclusive_group()
inputtype.add_argument('--fastx', '-f',
                   type= str,
                   help='''Path to fasta/fastq file.''')
inputtype.add_argument('--stdin', action='store_true',
                   help='''Fastx is coming from stdin stream.''')


filetype = parser.add_mutually_exclusive_group()
filetype.add_argument('--fa', action='store_true',
                   help='''Explicitly define as fasta input.''')
filetype.add_argument('--fq', action='store_true',
                   help='''Explicitly define as fastq input.''')


parser.add_argument('--verbose', '-v',
                   action='store_true',
                   help='''Provides jobname and pct complete messages while running (only applicable to extracting entries based on names).''', default=False)

names = parser.add_mutually_exclusive_group()

names.add_argument('--namesfile', '-n', type=str,
                   help='''Path to file with fasta/fastq entry names -- one name per line''')

names.add_argument('--names', '-c', type=str, help='Enter comma-separated names at command line with this flag')


names.add_argument('--annotatednamesfile', '-a', type=str,
                   help='''Path to file with fasta/fastq entry names and annotations.
There needs to be one name per line in first column.
The annotation in 2nd column (tab-sep) is added to the fasta header of the extracted sequence.
Does NOT have any effect when --exclude is being used -- works same as --namesfile.
Does NOT have effect when --indexes, --head, --minlen, --maxlen used...''')

names.add_argument('--regex', '-r', type=str, help='''Provide a regular expression. If it is in the record.description, the record will be printed.''')


parser.add_argument('--ignore_case', '-I', action='store_true', default=False, help='''Only has an effect in conjunction with --regex.
This simplifies writing the regex when you want to search for the presence of a word, but do not care how it appears.
For example: word, Word, wOrd, WOrD....
Normal regex: [Ww][Oo][Rr][Dd]
Regex with --ignore_case: word''')


parser.add_argument('--annotation2id', '-A', action='store_true', default=False,
                   help='''When using --annotatednamesfile, append annotation to record.id (first string before white space in fasta header).
This overrides default behavior of adding the annotation to the record.description after a tab.
With this flag, the annotation is appended with a '_' character to the record.id.
Any white space in the annotation is replaced with '-'.''')

parser.add_argument('--exclude', '-e', action='store_true', default=False,
                   help='''All sequences in fastx file EXCEPT names given will be returned.''')

parser.add_argument('--multiple', '-m', action='store_true', default=False,
                    help=''' Only use this if given names possibly occur multiple times in the file. This is unusual.
By default, it is assumed that a given name only occurs once, so it is removed from the list when it is encountered.
This flag will turn that off.''')

parser.add_argument("--minlen",
                    type=int, default=0,
                    help='''Use this to extract reads >= int given. Default: 0.''')

parser.add_argument("--maxlen",
                    type=int, default=30000000000,
                    help='''Use this to extract reads <= int given. Default: 30 billion.''')

parser.add_argument('--separate', '-s',
                   action='store_true',
                   help='''Works with -n and -c. Puts each fasta record in its own fasta file called name.fasta (with given name). This overrides default behavior of printing to stdout.''', default=False)

parser.add_argument('-i', '--indexes', type=str, default=False,
                    help=''' Use this if you just want to extract specific entries according to their order of appearance.
For example, you may want only the first entry (use -i 0), or first 2 (-i 0:2), or the first 5 and last 5 (-i 0:5,N-5:N).
Or you can take every even numbered entry by -i 0:N:2, for example.
In the latter 2 cases, you need to know the values of N and N-5 for now.
In general, you can give comma-separated indexes and can give "slices" as colon-separated integers.
Indexes are 0-based like Python. Slices will go up to but not include end of range.''')

parser.add_argument('--head', type=int, default=False,
                    help=''' Use this if you just want to extract the head (first N bases) of all sequences.''')

parser.add_argument('--subseqs', type=str, default=False,
                    help=''' Use this if you just want to extract only reads that contain sub-sequences/kmers in given file (sub-sequences assumed to be in first column of file regardless of how many other columns).
This can handle kmers (subseqs all of length k) as well as variable-length subseqs. This does not currently work with --exclude. (TODO)''')

parser.add_argument('--revsubseqs', action='store_true', default=False,
                    help='''Use this with --subseqs if the sub-sequence set should also contain the reverse complements. This forms a set, so there is no concern about representing a sub-sequence/kmer more than once.''')

parser.add_argument('--Nsubseqs', type=int, default=1,
                    help='''Use this with --subseqs to instruct how many times a subseq needs to be found to extract read.
                            This is NOT requiring them to be N distinct subseqs - the same subseq can occur N times.
                            Default: 1.''')

##
##parser.add_argument('-U', '--touppercase', type=int, default=False,
##                    help=''' Will ensure all letters in sequence returned are uppercase.''')
##parser.add_argument('-L', '--tolowercase', type=int, default=False,
##                    help=''' Will ensure all letters in sequence returned are lowercase.''')

args = parser.parse_args()
if args.regex: import re
if args.ignore_case: import string

############################################
'''           functions                '''
############################################

def find(string, char):
    return [i for i, ltr in enumerate(string) if ltr == char]

############################################
'''           check options              '''
############################################

## ensure that one of these are present
## quit if no -n or -c
    ##if not args.namesfile and not args.names:
    ##    print "-n or -c is required. Use -h for help."
    ##    quit()

## ensure one of these is used, and if stdin then also explicity specify filetype    
if not args.fastx and not args.stdin:
    print "Specify --stdin or -f/--fastx file.fx"
    quit()
if args.stdin and not (args.fa or args.fq):
    print "When using --stdin, filetype should be explicitly defined --fa or --fq."
    quit()

# file or stdin?
if args.stdin:
    fastxFile = sys.stdin
else:
    fastxFile = open(args.fastx)



# FASTA or FASTQ?
if args.fa:
    fastx = "fasta"
elif args.fq:
    fastx = "fastq"
elif args.fastx: ## can figure out from file (right now not from stdin due to non-seekability)
    line1 = fastxFile.next()[0]
    if line1[0] == ">":
        fastx = "fasta"
    elif line1[0] == "@":
        fastx = "fastq"
    fastxFile.seek(0)
else:
    print "Expected fasta or fastq. File given needs to be reformatted if user thinks it is."
    quit()



# job name
if args.verbose:
    pathInd = []
    if args.fastx:
        pathInd = find(args.fastx, "/")
    if len(pathInd) == 0:
        if args.stdin:
            jobname = "stdin"
        else:
            jobname = args.fastx
    else:
        jobname = args.fastx[max(pathInd)+1:]
    jobname += "+"
    if args.namesfile:
        pathInd = find(args.namesfile, "/")
        if len(pathInd) == 0:
            jobname += args.namesfile
        else:
            jobname += args.namesfile[max(pathInd)+1:]
    elif args.names:
        jobname += "commandline"




############################################
'''           execute                '''
############################################    

out = sys.stdout
msg = sys.stderr

if args.namesfile or args.names or args.annotatednamesfile:
    ## Make set of record IDs (names)
    names = set()
    if args.namesfile or args.annotatednamesfile:
        if args.namesfile:
            fname = args.namesfile
        else:
            fname = args.annotatednamesfile
            annotations = defaultdict(str)
        if fname == "-" or fname == "stdin":
            nfile = sys.stdin
        else:
            nfile = open(fname)
        for line in nfile:
            line = line.rstrip().split("\t")
            names.add(line[0])
            if args.annotatednamesfile:
                annotations[line[0]] = line[1]
            
    elif args.names:
        for name in args.names.split(","):
            names.add(name)
    setsize = len(names)
    pctdone = 0

    # Go through fasta file and return records that have IDs that match an element in set of names
    gatepct=10
    if not args.exclude:
        for record in SeqIO.parse(fastxFile, fastx):
            if record.id in names:
                if args.annotatednamesfile:
                    if args.annotation2id:
                        record.id = record.id + "_" + ("-").join(annotations[record.id].split())
                    else:
                        record.description = record.description + "\t" + annotations[record.id]
                if args.separate:
                    with open(record.id+".fasta", 'w') as f:
                        SeqIO.write(record, f, fastx)
                else:
                    SeqIO.write(record, out, fastx)
                if not args.multiple:
                    if args.annotatednamesfile and args.annotation2id:
                        record.id = ("_").join(record.id.split("_")[:-1])
                    names.remove(record.id)
                pctdone += 100*1.0/setsize
                if pctdone >= gatepct and args.verbose:
                    msg.write(str(pctdone)+"% complete...."+jobname+"\n")
                    gatepct += 10
    else:
        for record in SeqIO.parse(fastxFile, fastx):
            if record.id not in names:
                if args.separate:
                    with open(record.id+".fasta", 'w') as f:
                        SeqIO.write(record, f, fastx)
                else:
                    SeqIO.write(record, out, fastx)
                pctdone += 100*1.0/setsize
                if pctdone >= gatepct and args.verbose:
                    msg.write(str(pctdone)+"% complete...."+jobname+"\n")
                    gatepct += 10
            else: #it is in names, so no need to check for this one anymore
                if not args.multiple:
                    names.remove(record.id)

elif args.indexes:
    extract = []
    indexes = args.indexes.split(",")
    for idx in indexes:
        idxrange = idx.split(":")
        if len(idxrange) == 2:
            extract += range(int(idxrange[0]), int(idxrange[1]))
        elif len(idxrange) == 3:
            extract += range(int(idxrange[0]), int(idxrange[1]), int(idxrange[2]))
        elif len(idxrange) == 1:
            extract.append(int(idxrange[0]))
    extract.sort()
    i = 0
    j = 0
    nidx = len(extract)
    for record in SeqIO.parse(fastxFile, fastx):
        if j < nidx and i == extract[j]:
            SeqIO.write(record, out, fastx)
            j += 1
        i += 1

elif args.head:
    for record in SeqIO.parse(fastxFile, fastx):
        record = record[:args.head]
        SeqIO.write(record, out, fastx)


elif args.regex:
    if args.ignore_case:
        new_re = ''
        for ltr in args.regex.lower():
            if ltr in string.ascii_lowercase:
                new_re += '[' + ltr.upper() + ltr + ']'
            else:
                new_re += ltr
        args.regex = new_re
    regex = re.compile(args.regex)
    for record in SeqIO.parse(fastxFile, fastx):
        regex_present = len( re.findall( regex, record.description ) )
        if regex_present:
            SeqIO.write(record, out, fastx)

elif args.subseqs:
    kmers = [line.strip().split()[0].upper() for line in open(args.subseqs).readlines()]
    k = len(kmers[0])
    constant = sum([len(kmer)==k for kmer in kmers]) == len(kmers)
    if args.revsubseqs:
        intab = 'ACGT'
        outtab = 'TGCA'
        trantab = string.maketrans(intab, outtab)
        revkmers = [e.translate(trantab)[-1::-1] for e in kmers]
        kmers = set(kmers+revkmers)
    else:
        kmers = set(kmers)
    if constant:
        for record in SeqIO.parse(fastxFile, fastx):
            nfound = 0
            reclen = len(record.seq)
            recseq = str(record.seq)
            for i in range(0,reclen-k+1):
                if recseq[i:i+k] in kmers:
                    nfound += 1
                    if nfound >= args.Nsubseqs:
                        SeqIO.write(record, out, fastx)
                        break
    else:
        subseqs = list(kmers) ## can be diff lengths
        for record in SeqIO.parse(fastxFile, fastx):
            nfound = 0
            for subseq in subseqs:
                subseqlen = len(subseq)
                reclen = len(record.seq)
                recseq = str(record.seq)
                for i in range(0,reclen-subseqlen+1):
                    if recseq[i:i+subseqlen] == subseq:
                        nfound += 1
                        if nfound >= args.Nsubseqs:
                            SeqIO.write(record, out, fastx)
                            break
        
    


## NEW ELIFs SHOULD BE ABOVE THIS LINE    
elif args.minlen or args.maxlen:
    for record in SeqIO.parse(fastxFile, fastx):
        if len(record) >= args.minlen and len(record) <= args.maxlen:
            SeqIO.write(record, out, fastx)
                   

fastxFile.close()
out.close()


