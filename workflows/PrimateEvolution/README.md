# Primate Evolution

The goal of this workflow is to perform comparative evolutionary analyses on several priamte genomes. Specifically, we will be exploring two loci that are of interest to primate immunity (SAMD9 and HLA). A self-contained Snakemake file and configfile are also included for reference. These use Singularity to run the docker image of pggb.   

## Data collection
Four high-quality primate genomes were selected for this analysis: GRC38 (Human), Clint (Chimpanzee), Mhudiblu (Orangutan), and Kamilah (Gorilla). Genome accessions and ftp links are shown below. Reference genomes were all downloaded from ncbi at these links. 


```
	"assemblies":{
		"GRC38":{
			"genbank_accession":"GCA_000001405.28",
			"species":"human",
			"fn_fna":"GCA_000001405.28_GRCh38.p13_genomic.fna.gz",
			"fn_ftp":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz",
		},
		"Clint_PTRv2":{
			"genbank_accession":"GCA_002880755.3",
			"species":"chimp",
			"fn_fna":"GCF_002880755.1_Clint_PTRv2_genomic.fna.gz",
			"fn_ftp":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/880/755/GCF_002880755.1_Clint_PTRv2/GCF_002880755.1_Clint_PTRv2_genomic.fna.gz"
		},
		"Mhudiblu_PPA_v2":{
			"genbank_accession":"GCA_013052645.3",
			"species":"bornean_orangutan",
			"fn_fna":"GCA_013052645.3_Mhudiblu_PPA_v2_genomic.fna.gz",
			"fn_ftp":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/052/645/GCA_013052645.3_Mhudiblu_PPA_v2/GCA_013052645.3_Mhudiblu_PPA_v2_genomic.fna.gz"
		},
		"Kamilah_GGO_v0":{
			"genbank_accession":"GCA_008122165.1",
			"species":"gorilla",
			"fn_fna":"GCA_008122165.1_Kamilah_GGO_v0_genomic.fna.gz",
			"fn_ftp":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/122/165/GCA_008122165.1_Kamilah_GGO_v0/GCA_008122165.1_Kamilah_GGO_v0_genomic.fna.gz"
		},
	}
```
We download and index these fastqs using the following rules
```

	rule download:
		output:
			fa_gz_out="data/input_genomes/{g}/{fn_fna}.fna.gz"
		run:
			ftp_path = config['assemblies'][wildcards.g]["fn_ftp"]
			shell("wget -O {fn_out} {ftp_path}".format(fn_out=output.fa_gz_out,
													   ftp_path=ftp_path))

	rule idx_fa:
		input:
			get_fa
		output:
			idx_out="data/input_genomes/{g}/{fn_fna}.fna.gz.fai"
		run:
			fa_gz = input[0]
			fa = fa_gz.replace(".gz","")
			shell("gunzip {fa_gz_in}".format(fa_gz_in=fa_gz))
			shell("bgzip {fa_gz_in}".format(fa_gz_in=fa))
			shell("samtools faidx {fa_bz_out}".format(fa_bz_out=fa_gz))
                                    
```


## Pangenome Sequence Naming and Partitioning
The two loci we are interested (SAMD9 and HLA) reside on chr7 and chr6 respectively. We thus extract these homologous chromosomes (termed "communities" below) and create individual fastas for each. 

```
    "contig_communities":{
        "asssembly_order":["GRC38",
                           "Clint_PTRv2",
                           "Mhudiblu_PPA_v2",
                           "Kamilah_GGO_v0"],
        "communities":{
            "chr6": ["CM000668.2",
                     "NC_036885.1",
                     "CM023033.3",
                     "CM017853.1"],
            "chr7": ["CM000669.2",
                     "NC_036886.1",
                     "CM023034.3",
                     "CM017854.1"]
        }
    }
```
Rules to extract the loci of interest and concatenate them into a properly named input.
```

	"""
	https://github.com/pangenome/PanSN-spec
	uses this hack https://github.com/ekg/fastix
	"""
	rule make_input_fa:
		input:
			get_fa_idxs
		output:
			fa_out="data/{contig}.input.fa.gz"
		run:
			fastix="/global/home/users/psudmant/code/fastix/target/debug/fastix"
			shell("> {fa_out}".format(fa_out=output.fa_out.replace(".gz","")))
			assemblies = config['contig_communities']['asssembly_order']
			for i, contig in enumerate(config['contig_communities']['communities'][wildcards.contig]):
				curr_assembly = assemblies[i] 
				fn_fa=config['assemblies'][curr_assembly]['fn_fna']
				shell("samtools faidx "
					  "data/input_genomes/{g}/{fn_fa} "
					  "{contig} | "
					  "{fastix} "
					  "--prefix {g}# "
					  "/dev/stdin "
					  ">>{fa_out}"
					  "".format(fn_fa=fn_fa,
								fastix=fastix,
								contig=contig,
								g=curr_assembly,
								fa_out=output.fa_out.replace(".gz","")))
			shell("bgzip {fa_out}".format(fa_out = output.fa_out.replace(".gz","")))
			shell("samtools faidx {fa_out}".format(fa_out = output.fa_out))


```

## Divergence estimation
In the case of these primates we know the relative divergences with the most recent common ancestor ~10 million years ago and a maximum divergence of ~3.5% genome wide. We thus set ourdivergence conservatively at 5%. 

## Pangenome graph building
We build the pangenome graph using the singularity image of the PGGB. Complete instructions on using singularity with docker images can be found [here](https://github.com/pangenome/pggb#singularity). The rule build

```
	rule run_pggb:
		input:
			fa_input="data/{contig}.input.fa.gz"
		output:
			out="data/out.{contig}.{seg_len}/multiqc_config.yaml"
		run:
			pggb_path=config['pggb_path']
			PWD=config['PWD'] 
			outdir = "data/out.{contig}.{seg_len}".format(contig=wildcards.contig,
														  seg_len=wildcards.seg_len) 
			cmd = ("singularity "
				   "run -B {PWD}/data:/data "
				   "-H {PWD} "
				   "{pggb_path} "
				   "\"pggb -i /{fa_input} "
				   "-p 95 "
				   "-s {seg_len} "
				   "-n 4 "
				   "-t 24 "
				   "-o /{outdir} "
				   #"-M -C cons,100,1000,10000 -m\""
				   "-M -m\""
				   "".format(pggb_path=pggb_path,
							 fa_input = input.fa_input,
							 seg_len = wildcards.seg_len,
							 outdir=outdir,
							 PWD=PWD))
			shell(cmd)

```


## Analyses
We are now ready to analyze our genomes! Let's begin by vizualizing our loci of interest. We do this using several ODGI commands. First we need to extract subgraphs at our locus of interest (`odgi extact`). We can use the genome coordinates from any of our references (see below). We then need to resort this subgraph (`odgi sort`). Finally, we can use `odgi viz` to vizualize these subgraphs. 

```
    "extract_loci":{

        "SAMD9_Clint":{
            "graph":"data/out.chr7.10000/chr7.input.fa.gz.e6f73d2.e34d4cd.20398a9.smooth.final.og",
            "graph_path":"out.chr7.10000",
            "locus":"Clint_PTRv2#NC_036886.1:89120280-89210088"
        },

        "SAMD9_GRC38":{
            "graph":"data/out.chr7.10000/chr7.input.fa.gz.e6f73d2.e34d4cd.20398a9.smooth.final.og",
            "graph_path":"out.chr7.10000",
            "locus":"GRC38#CM000669.2:93060633-93150835"
        },
        "HLA_GRC38":{
            "graph":"data/out.chr6.10000/chr6.input.fa.gz.e6f73d2.e34d4cd.20398a9.smooth.final.og",
            "graph_path":"out.chr6.10000",
            "locus":"GRC38#CM000668.2:29657092-33192467"
        }

    }
```
Rules for using ODGI to extract, sort, and vizualize specific loci of our pangenomes
```
	rule extract_locus:
		input:
			"data/out.chr6.10000/multiqc_config.yaml",
			"data/out.chr7.10000/multiqc_config.yaml"
		output:
			png="data/{contig}/{target}/{target}.sorted.png",
			subgraph="data/{contig}/{target}/{target}.og",
			sorted_subgraph="data/{contig}/{target}/{target}.sorted.og"
		run:
			pggb_path=config['pggb_path']
			PWD=config['PWD'] 
			locus = config["extract_loci"][wildcards.target]["locus"] 
			graph = config["extract_loci"][wildcards.target]["graph"] 
			
			cmd = ("singularity "
				   "run -B {PWD}/data:/data "
				   "-H {PWD} "
				   "{pggb_path} "
				   "\"odgi extract -E -i /{graph} -o /{subgraph} -r {locus}\""
				   "".format(pggb_path=pggb_path,
							 PWD=PWD,
							 locus=locus,
							 graph=graph,
							 subgraph=output.subgraph))
			shell(cmd)
			cmd = ("singularity "
				   "run -B {PWD}/data:/data "
				   "-H {PWD} "
				   "{pggb_path} "
				   "\"odgi sort -O -Y -i /{subgraph} -o /{sorted_subgraph}\""
				   "".format(pggb_path=pggb_path,
							 PWD=PWD,
							 subgraph=output.subgraph,
							 sorted_subgraph=output.sorted_subgraph))
			shell(cmd)

			cmd = ("singularity "
				   "run -B {PWD}/data:/data "
				   "-H {PWD} "
				   "{pggb_path} "
				   "\"odgi viz -i /{sorted_subgraph} -o /{png}\""
				   "".format(pggb_path=pggb_path,
							 PWD=PWD,
							 sorted_subgraph=output.sorted_subgraph,
							 png=output.png))
			shell(cmd)

```
Results of our vizualization. Note that we have extracted the subgraph of the SAMD9 locus using the human linear genome coordinates and the chimpanzee linear genome coordinates with identical results.  

![An ODGI viz visualization of the SAMD9 locus extracted from the Clint](data/out.chr7.10000/SAMD9_Clint/SAMD9_Clint.sorted.png)
![An ODGI viz visualization of the SAMD9 locus extracted against GRC38](data/out.chr7.10000/SAMD9_GRC38/SAMD9_GRC38.sorted.png)
![An ODGI viz visualization of the HLA locus](data/out.chr6.10000/HLA_GRC38/HLA_GRC38.sorted.png)



