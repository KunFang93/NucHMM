'''Commonly used the data store here'''
hg19 = {'chr1': 249250621, 'chr2':243199373, 'chr3':198022430,'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067,
        'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
        'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248,
        'chr19':59128983, 'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566}

hg38 = {'chr1': 249698942, 'chr2':242508799, 'chr3':198450956,'chr4': 190424264, 'chr5': 181630948, 'chr6': 170805979,
        'chr7': 159345973, 'chr8': 145138636, 'chr9': 138688728, 'chr10': 133797422, 'chr11': 135186938, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 108136338, 'chr15': 102439437, 'chr16': 92211104, 'chr17': 83836422, 'chr18': 80373285,
        'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983, 'chr22': 51857515, 'chrX': 156040895, 'chrY': 57264655}

val2chr = {1: 'chr1', 2: 'chr2', 3: 'chr3', 4: 'chr4', 5: 'chr5', 6: 'chr6', 7: 'chr7', 8: 'chr8', 9: 'chr9',
           10: 'chr10', 11: 'chr11', 12: 'chr12', 13: 'chr13', 14: 'chr14', 15: 'chr15', 16: 'chr16', 17: 'chr17',
           18: 'chr18', 19: 'chr19', 20: 'chr20', 21: 'chr21', 22: 'chr22', 23: 'chrX', 24: 'chrY'}

spe_colors = {'1':'gold','2':'blue', '3':'orange', '4':'green', '5':'red', '6':'purple', '7':'brown', '8':'pink',
              '9':'gray', '10':'olive','11':'cyan','12':'lightblue','13':'peru', '14':'navy','15':'limegreen',
              '16':'darkcyan','17':'orchid','18':'saddlebrown','19':'silver','20':'black'}

info2strand = {'-':'back', '+':'for'}

# H3K4me1/3, H3K27ac, H3K79me2/3, H3K36me3 Location information come from
# https://www.cell.com/trends/biochemical-sciences/pdf/S0968-0004(17)30189-5.pdf
# H3K9me3, H3K27me3 location information come from paper
# High-Resolution Profiling of Histone Methylations in the Human Genome
# H3K9me3 actually has no specific location but has some functional meanings in Distal (10kb) according to paper

Histone_location = {'H3K4me1':['Distal','Proximal'],'H3K4me2':['Promoter','5\'-Genebody'],'H3K4me3':['Promoter','5\'-Genebody'],
                    'H3K27ac':['Distal','Proximal','Promoter'],'H3K79me2':['5\'-Genebody'],
                    'H3K79me3':['5\'-Genebody'],'H3K36me3':['3\'-Genebody'],
                    'H3K9me3':['Distal','Proximal','Promoter','Genebody'],
                    'H3K27me3':['Distal','Proximal','Promoter','Genebody']}




