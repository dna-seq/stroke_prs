import time
vcf_file_path = 'C:/dev/python/openCravatPlugin/vcf/antonkulaga.vcf'

chr_name_ind = 0
chr_position_ind = 1
effect_allele_ind = 2
effect_weight_ind = 4
hm_chr_ind = 7
hm_pos_ind = 8
ref_ind = 10

CHROM_ind = 0
POS_ind = 1
REF_ind = 3
ALT_ind = 4
default_ind = 9
FILTER_ind = 6

referens_count = 0
total_snps = 0

def parse_prs(path):
    global referens_count
    global total_snps
    res = dict()
    with open(path) as f:
        skip = True
        for row in f:
            row = row[0:-1]
            if row.startswith("chr_name"):
                skip = False
                continue
            if skip:
                continue
            parts = row.split("\t")
            key = parts[hm_chr_ind]+parts[hm_pos_ind]
            if key == "":
                continue

            ref = parts[ref_ind]
            isref = ref == parts[effect_allele_ind]
            if isref:
                referens_count += 1
            total_snps += 1
            res[key] = [parts[effect_allele_ind], float(parts[effect_weight_ind]), isref, True]
    return res

def get_alleles(ref, alt, info):
    if info == "1/1":
        return alt, alt
    if info == "0/1":
        return alt, ref
    if info == "0/2":
        return alt.split(",")[1], ref
    if info == "1/2":
        p = alt.split(",")
        return p[0], p[1]

start = time.time()
prs = parse_prs("PGS002724_hmPOS_GRCh38_plus_ref.txt")
sum = 0.0
count = 0
row_count = 0
with open(vcf_file_path) as f:
    skip = True
    for row in f:
        row = row[0:-1]
        row_count += 1
        if row_count % 100000 == 0:
            print("Processed:", row_count)
        if row.startswith("#CHROM"):
            skip = False
            continue
        if skip:
            continue
        parts = row.split("\t")
        if parts[FILTER_ind] == "RefCall" or parts[default_ind][0:3] == "0/0":
            continue
        key = parts[CHROM_ind]+parts[POS_ind]
        alelle_weight = prs.get(key)
        if alelle_weight is None:
            continue
        alelle_weight[3] = False
        a, b = get_alleles(parts[REF_ind], parts[ALT_ind], parts[default_ind][0:3])
        flag = False
        if a == alelle_weight[0]:
            sum += alelle_weight[1]
            flag = True
        if b == alelle_weight[0]:
            sum += alelle_weight[1]
            flag = True
        if flag:
            count += 1

for key in prs:
    alelle_weight = prs.get(key)
    if alelle_weight[2] and alelle_weight[3]:
        sum += alelle_weight[1]*2
        count += 1

avg = 0
if count > 0:
    avg = sum/count
print("sum:", sum, "count:", count, "avg:", avg)
print("time:", time.time() - start)
print("Overall referense snps:", referens_count, "/", total_snps)







