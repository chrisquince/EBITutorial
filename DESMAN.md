## DESMAN TARA Tutorial

We are going to locate ecotypes in one of the TARA MAGs, a SAR11 relative TARA_PSW_MAG_00074.
We begin by downloading the Ebame3 repo onto our VMs which has some useful scripts.
Log into the computers:

```
ssh ubuntu@myVM
```

```
cd ~/repos
git clone https://github.com/chrisquince/Ebame3.git
```

We will also need to install the maps package into R:

```
sudo R
>install.packages('maps')
>q()
```

Create a directory to work in:

```
cd ~/Projects
mkdir DESMANTutorial
cd DESMANTutorial
```

Download the pre-computed frequency file on the single copy core genes. 

```
wget https://desmantutorial.s3.climb.ac.uk/TARA_PSW_MAG_00074_scg.freq
```

How many core COGs did this MAG contain?

### Variant selection

Then we run the variant filter to identify variant positions. We also filter for non-core cogs '-c' and set a sample coverage threshold of '-m 1.0':

```
python $DESMAN/desman/Variant_Filter.py TARA_PSW_MAG_00074_scg.freq -p -o TARA_PSW_MAG_00074_scg -f 25.0 -c -sf 0.80 -t 2.5 -m 1.0
```

The '-p' option uses the slower but more precision minor variant frequency optimisation.

How many COGs survived the filtering? How many variant positions are there?

### Haplotype inference

Then we run the haplotype inference on the variant positions that are selected.

```
stub=TARA_PSW_MAG_00074_scg
varFile=${stub}sel_var.csv

eFile=${stub}tran_df.csv

for g in 1 2 3 4  
do
for r in 0 1 2 3 4
    do
        echo $g:
        desman $varFile -e $eFile -o ${stub}_${g}_${r} -g $g -s $r -m 1.0 -i 100 > ${stub}_${g}_${r}.out&
    done
done
```

Create a deviance plot for this MAG:

```
cat */fit.txt | cut -d"," -f2- > Dev.csv
sed -i '1iH,G,LP,Dev' Dev.csv 
```

Which we can then visualise:

```
$DESMAN/scripts/PlotDev.R -l Dev.csv -o Dev.pdf
```

![Posterior deviance](../Figures/Dev.pdf)

This is typical of real data. The number of true strains is not clear. So we try our heuristic estimate:

```
python /home/ubuntu/repos/DESMAN/scripts/resolvenhap.py TARA_PSW_MAG_00074_scg
```

Should see something like:

3,3,0,0.0251628182356,TARA_PSW_MAG_00074_scg_3_0/Filtered_Tau_star.csv

Suggesting three strains are present.


Now we need the TARA sample meta data:

```
wget https://desmantutorial.s3.climb.ac.uk/TARA-samples.csv
```

The haplotypes themselves are encoded in the file:

```
TARA_PSW_MAG_00074_scg_3_0/Filtered_Tau_star.csv
```

The haplotypes differ on the core genes:

```
python $DESMAN/scripts/validateSNP2.py TARA_PSW_MAG_00074_scg_3_0/Filtered_Tau_star.csv TARA_PSW_MAG_00074_scg_3_0/Filtered_Tau_star.csv
```

We can look at abundances and reproducibility of the strains:

```
python /home/ubuntu/repos/DESMAN/scripts/taucomp.py TARA_PSW_MAG_00074_scg_3_0/Gamma_starR.csv TARA_PSW_MAG_00074_scg_3_*/Filtered_Tau_star.csv
```

The position encoding is ACGT so what are the base predictions at each variant position? 
We can turn these into actual sequences with the following commands.
First, we obtain a list of core genes based on the 36 COGs. These we can get from the core COG filtered file:

```
cut -d"," -f1 TARA_PSW_MAG_00074_scgcogf.csv | sed '1d' | sort | uniq > core_genes.txt
```


No we need the original cluster sequences and COGs:

```
wget https://desmantutorial.s3.climb.ac.uk/TARA_PSW_MAG_00074.fa
wget https://desmantutorial.s3.climb.ac.uk/TARA_PSW_MAG_00074F.cogs
```

Select the core cogs:

```
$DESMAN/scripts/SelectContigsPos.pl ~/repos/MAGAnalysis/scgs.txt < TARA_PSW_MAG_00074F.cogs > TARA_PSW_MAG_00074F_core.cogs
```

Create an output directory:

```
mkdir SCG_Fasta_3_0
```

And get the haplotype sequences:

```
python $DESMAN/scripts/GetVariantsCore.py TARA_PSW_MAG_00074.fa TARA_PSW_MAG_00074F_core.cogs ./TARA_PSW_MAG_00074_scg_3_0/Filtered_Tau_star.csv core_genes.txt -o SCG_Fasta_3_0/
```

And now we can see if these strains are significantly associated with location:

```
Rscript ~/repos/Ebame3/scripts/StrainANOVA.R -c TARA_PSW_MAG_00074_scg_3_0/Gamma_star.csv -s TARA-samples.csv -t TARA_PSW_MAG_00074
```

TARA_PSW_MAG_00074,1,0.254687,0.007304,22.544615
TARA_PSW_MAG_00074,2,0.572026,0.003279,24.734444
TARA_PSW_MAG_00074,3,0.173287,0.049916,16.924231

and generate the world maps:

```
~/repos/Ebame3/scripts/TARA_HaploMap.R -g TARA_PSW_MAG_00074_scg_3_0/Gamma_star.csv -c ~/repos/Ebame3/data/TARA-clusters.txt -s ~/repos/Ebame3/data/TARA-samplesR.txt 
```

![H0](../Figures/WorldMap_H0.pdf)

![H1](../Figures/WorldMap_H1.pdf)

![H2](../Figures/WorldMap_H2.pdf)

### Accessory gene inference

We start from the base frequencies on all genes.

```
mkdir ~/Projects/DESMANTutorial/AllFreq
cd ~/Projects/DESMANTutorial/AllFreq
wget https://desmantutorial.s3.climb.ac.uk/TARA_PSW_MAG_00074.freq
```

Then we find SNVs notice the lack of '-p' option for efficiency. With more time I would probably apply that.

```
python $DESMAN/desman/Variant_Filter.py TARA_PSW_MAG_00074.freq -o TARA_PSW_MAG_00074 -m 0. -f 25.0
```

Then we need gene coverages:

```
python $DESMAN/scripts/CalcGeneCov.py TARA_PSW_MAG_00074.freq > TARA_PSW_MAG_00074_gene_cov.csv
```

Then we calculate core gene coverages:

```
python $DESMAN/scripts/CalcDelta.py TARA_PSW_MAG_00074_gene_cov.csv ../core_genes.txt TARA_PSW_MAG_00074_core
```

and finally:

```
stub=TARA_PSW_MAG_00074
sel_run=../TARA_PSW_MAG_00074_scg_3_0/

python $DESMAN/desman/GeneAssign.py ${stub}_coremean_sd_df.csv ${sel_run}/Gamma_star.csv ${stub}_gene_cov.csv ${sel_run}/Eta_star.csv -m 20 -v ${stub}sel_var.csv -o $stub --assign_tau

```

But we will not run this ourselves it is too slow:

```
cd ~/Projects/DESMANTutorial
rm -r AllFreq
cp ~/AllFreq.tar.gz .
tar -xvzf AllFreq.tar.gz
```

We can look at the overlap in terms of accessory genes:
```
cd AllFreq
python EtaGamma.py TARA_PSW_MAG_00074etaS_df.csv ../TARA_PSW_MAG_00074_scg_3_0/Gamma_star.csv ../TARA_PSW_MAG_00074_scg_3_0/Filtered_Tau_star.csv
```