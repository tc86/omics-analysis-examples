#!/usr/bin/env python3
import argparse, gzip, sys
from pathlib import Path
from bisect import bisect_right
from collections import defaultdict
import pandas as pd

def open_maybe_gzip(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def parse_sample_from_header(vcf_path):
    with open_maybe_gzip(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 10:
                    return cols[9]
                break
    return Path(vcf_path).stem

def parse_variant_counts(format_str, sample_str):
    keys = format_str.split(":")
    vals = sample_str.split(":")
    fmap = dict(zip(keys, vals))
    if "AD" in fmap and fmap["AD"] not in (".", ""):
        parts = [p for p in fmap["AD"].split(",") if p not in (".","")]
        if len(parts) >= 2:
            try:
                return int(parts[0]), int(parts[1])
            except ValueError: pass
    if "RO" in fmap and "AO" in fmap and fmap["RO"] not in (".","") and fmap["AO"] not in (".",""):
        try:
            return int(fmap["RO"]), int(fmap["AO"].split(",")[0])
        except ValueError: pass
    if "AF" in fmap and "DP" in fmap:
        try:
            dp = int(fmap["DP"]); af = float(fmap["AF"].split(",")[0])
            alt = int(round(af*dp)); ref = dp - alt
            return ref, alt
        except ValueError: pass
    return None

def load_cnv_info(cnv_path):
    df = pd.read_csv(cnv_path, sep="\t")
    required = {"sample_id","chrom","start","end","major_cn","minor_cn","normal_cn"}
    if not required.issubset(df.columns):
        sys.exit("ERROR: CNV table missing required cols")
    def norm_chrom(x): return x[3:] if str(x).lower().startswith("chr") else str(x)
    df["chrom"] = df["chrom"].map(norm_chrom)
    cnv_index = defaultdict(dict); tc_map={}
    if "tumour_content" in df.columns:
        for sid,g in df.groupby("sample_id"):
            vals = pd.to_numeric(g["tumour_content"],errors="coerce").dropna()
            tc_map[sid]=float(vals.median()) if len(vals) else None
    for sid,gsid in df.groupby("sample_id"):
        for chrom,gc in gsid.groupby("chrom"):
            gc=gc.sort_values("start"); starts=gc["start"].tolist()
            intervals=list(zip(gc["start"],gc["end"],gc["major_cn"],gc["minor_cn"],gc["normal_cn"]))
            cnv_index[sid][chrom]={"starts":starts,"intervals":intervals}
    return cnv_index,tc_map

def lookup_cn(cnv_index,sid,chrom,pos):
    chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
    table=cnv_index.get(sid,{}).get(chrom)
    if not table: return (None,None,None)
    i=bisect_right(table["starts"],pos)-1
    if i<0: return (None,None,None)
    s,e,maj,mino,norm=table["intervals"][i]
    if s<=pos<=e: return int(maj),int(mino),int(norm)
    return (None,None,None)

def iter_vcf_records(vcf_path):
    with open_maybe_gzip(vcf_path) as f:
        for line in f:
            if line.startswith("#") or not line: continue
            cols=line.rstrip().split("\t")
            if len(cols)<10: continue
            chrom,pos,_,ref,alt,_,_,_,fmt,sample=cols[:10]
            try: pos=int(pos)
            except: continue
            yield chrom,pos,ref,alt.split(",")[0],fmt,sample

def build_pyclone_input(vcf_dir,cnv_path,out_tsv,require_all=False,default_tc=1.0):
    cnv_index,tc_map=load_cnv_info(cnv_path); rows=[]; samples=[]
    for vcf in sorted(Path(vcf_dir).glob("*.vcf*")):
        sid=parse_sample_from_header(vcf); samples.append(sid)
        for chrom,pos,ref,alt,fmt,sample_field in iter_vcf_records(vcf):
            counts=parse_variant_counts(fmt,sample_field)
            if counts is None: continue
            ref_c,alt_c=counts; mut=f"{chrom}:{pos}:{ref}>{alt}"
            maj,mino,norm=lookup_cn(cnv_index,sid,chrom,pos)
            if maj is None: continue
            tc=tc_map.get(sid,default_tc)
            rows.append({"mutation_id":mut,"sample_id":sid,"ref_counts":ref_c,
                         "alt_counts":alt_c,"major_cn":maj,"minor_cn":mino,
                         "normal_cn":norm,"tumour_content":tc})
    if not rows: sys.exit("No rows produced")
    df=pd.DataFrame(rows)
    if require_all:
        keep=df.groupby("mutation_id")["sample_id"].nunique().eq(len(set(samples)))
        keep_ids=df.groupby("mutation_id").filter(lambda x: len(set(x["sample_id"]))==len(set(samples)))["mutation_id"]
        df=df[df["mutation_id"].isin(keep_ids)]
    cols=["mutation_id","sample_id","ref_counts","alt_counts","major_cn","minor_cn","normal_cn","tumour_content"]
    df[cols].to_csv(out_tsv,sep="\t",index=False)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--vcf-dir",required=True)
    ap.add_argument("--cnv-info",required=True)
    ap.add_argument("--output",required=True)
    ap.add_argument("--require-all-samples",action="store_true")
    ap.add_argument("--default-tumour-content",type=float,default=1.0)
    a=ap.parse_args()
    build_pyclone_input(a.vcf_dir,a.cnv_info,a.output,a.require_all_samples,a.default_tumour_content)

if __name__=="__main__": main()
