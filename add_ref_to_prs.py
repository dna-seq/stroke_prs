import sqlite3

chr_name_ind = 0
chr_position_ind = 1
effect_allele_ind = 2
effect_weight_ind = 4
hm_chr_ind = 7
hm_pos_ind = 8

dbsnp_mismatch_count = 0
dbsnp_total_count = 0
referens_count = 0
total_snps = 0

conn = sqlite3.connect("D:/dev/oakVar/modules/annotators/dbsnp/data/dbsnp.sqlite")
cursor = conn.cursor()

def getRef(chrom, pos):
    global cursor
    global dbsnp_mismatch_count
    global dbsnp_total_count
    chrom = "chr"+chrom
    cursor.execute(f"SELECT ref FROM {chrom} WHERE pos = {pos}")
    row = cursor.fetchone()
    dbsnp_total_count += 1
    if row is None or len(row) == 0:
        dbsnp_mismatch_count += 1
        return ""

    return row[0]


def parse_prs(path, new_path):
    global referens_count
    global total_snps
    fwrite = open(new_path, "w")
    with open(path) as f:
        skip = True
        for row in f:
            row = row[0:-1]
            if row.startswith("chr_name"):
                fwrite.write(row+"\tref\n")
                skip = False
                continue
            if skip:
                fwrite.write(row+"\n")
                continue
            parts = row.split("\t")
            key = parts[hm_chr_ind]+parts[hm_pos_ind]
            if key == "":
                fwrite.write(row+"\t\n")
                continue
            ref = getRef(parts[hm_chr_ind], parts[hm_pos_ind])
            isref = (ref == parts[effect_allele_ind])
            if isref:
                referens_count += 1
            total_snps += 1
            if total_snps % 10000 == 0:
                print(total_snps)
            fwrite.write(row + "\t"+ref+"\n")
    fwrite.close()

parse_prs("PGS002724_hmPOS_GRCh38.txt", "PGS002724_hmPOS_GRCh38_plus_ref.txt")

cursor.close()
conn.close()

print("Finish")