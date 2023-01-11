import polars as pl
import time

start = time.time()
prs = pl.read_csv("PGS002724_hmPOS_GRCh38_plus_ref.txt", sep="\t", comment_char="#", ignore_errors=True)
prs = prs.with_columns([pl.col("hm_chr").cast(pl.Utf8), pl.col("hm_pos").cast(pl.Utf8)])
prs = prs.with_columns([(pl.col("hm_chr") + pl.col("hm_pos")).alias("key_right")])
vcf = pl.read_csv("C:/dev/python/openCravatPlugin/vcf/antonkulaga.vcf", sep="\t",
                  comment_char="#", ignore_errors=True, has_header=False,
                  new_columns=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","default"])

vcf = vcf.filter((pl.col("FILTER") == "PASS") & (pl.col("default").str.starts_with("0/0")).is_not())
vcf = vcf.with_columns([pl.col("CHROM").cast(pl.Utf8), pl.col("POS").cast(pl.Utf8)])
vcf = vcf.with_columns([(pl.col("CHROM") + pl.col("POS")).alias("key_left")])

vcf_1_1 = vcf.filter(pl.col("default").str.starts_with("1/1"))
vcf_1_1 = vcf_1_1.with_columns([pl.col("ALT").alias("A"), pl.col("ALT").alias("B")])

vcf_0_1 = vcf.filter(pl.col("default").str.starts_with("0/1"))
vcf_0_1 = vcf_0_1.with_columns([pl.col("ALT").alias("A"), pl.col("REF").alias("B")])

vcf_0_2 = vcf.filter(pl.col("default").str.starts_with("0/2"))
vcf_0_2 = vcf_0_2.with_columns([pl.col("ALT").str.split(by=",").arr.last().alias("A"), pl.col("REF").alias("B")])

vcf_1_2 = vcf.filter(pl.col("default").str.starts_with("1/2"))
vcf_1_2 = vcf_1_2.with_columns([pl.col("ALT").str.split(by=",").arr.first().alias("A"), pl.col("ALT").str.split(by=",").arr.last().alias("B")])

vcf = vcf_1_1.vstack(vcf_0_1).vstack(vcf_0_2).vstack(vcf_1_2)
unite = vcf.join(prs, left_on="key_left", right_on="key_right")
res1 = unite.filter(pl.col("A") == pl.col("effect_allele")).select(pl.col("effect_weight")).sum()
res2 = unite.filter(pl.col("B") == pl.col("effect_allele")).select(pl.col("effect_weight")).sum()

print(res1 + res2)
print("time:", time.time() - start)