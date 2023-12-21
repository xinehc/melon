# Melon
Melon: **me**tagenomic **lon**g-read-based taxonomic identification and quantification

## Quick start
### Installation
Create a new conda environment and install:
```bash
conda create -n melon -c conda-forge -c bioconda melon
conda activate melon
```

### Database setup
Download either the NCBI or the GTDB database:
```bash
## GTDB
# wget -q --show-progress https://figshare.com/ndownloader/files/42694702/database.tar.gz

## NCBI
wget -q --show-progress https://figshare.com/ndownloader/files/42694969/database.tar.gz
tar -zxvf database.tar.gz
```

Index the files: 
```bash
## If you encounter memory issue please consider manually lowering cpu_count or simply set cpu_count=1
cpu_count=$(python -c 'import os; print(os.cpu_count())')

diamond makedb --in database/prot.fa --db database/prot --quiet
ls database/nucl.*.fa | sort | xargs -P $cpu_count -I {} bash -c '
    filename=${1%.fa*}; \
    filename=${filename##*/}; \
    minimap2 -x map-ont -d database/$filename.mmi ${1} 2> /dev/null' - {}

## remove unnecessary files to save space
rm -rf database/*.fa
```

### Run Melon
> [!NOTE]  
> Melon takes **quality-controlled** and **decontaminated** long reads as input. We suggest to remove low-quality raw reads before running Melon with e.g., `nanoq -q 10 -l 1000` (minimal quality score 10; minimal read length 1,000 bp). If your sample is known to have a large proportion of human DNAs or known eukaryotes/viruses, please consider removing them via proper mapping. If the origin of contamination is unknown, or if you want to estimate the mean genome size of prokaryotes, you may consider enabling the simple pre-filtering module. See [Run Melon with pre-filtering of non-prokaryotic reads](#run-melon-with-pre-filtering-of-non-prokaryotic-reads) for more details.

We provide an example file comprising 10,000 quality-controlled (processed with `Porechop` and `nanoq`), prokaryotic reads (fungal and other reads removed with `minimap2`) randomly selected from the R10.3 mock sample of [Loman Lab Mock Community Experiments](https://lomanlab.github.io/mockcommunity/r10.html).

```bash
wget -q --show-progress https://figshare.com/ndownloader/files/42847672/example.fa.gz
melon example.fa.gz -d database -o .
```

You should see:
```
INFO: Estimating genome copies ...
INFO: ... found 27.5 copies of genomes (bacteria: 27.5; archaea: 0).
INFO: Assigning taxonomy ...
INFO: ... found 8 unique species (bacteria: 8; archaea: 0).
INFO: Done.
```

The output file `*.tsv` contains the estimated genome copies for individual species, their corresponding relative abundances and gap-compressed ANI (average nucleotide identity between marker-gene-containing reads and reference genome clusters) values:
```
...    species                               copy     abundance              identity
...    287|Pseudomonas aeruginosa            2.125    0.07727272727272727    0.9570235294117647
...    96241|Bacillus spizizenii             2.875    0.10454545454545454    0.9617130434782607
...    1351|Enterococcus faecalis            3.0      0.10909090909090909    0.9616333333333333
...    28901|Salmonella enterica             3.125    0.11363636363636363    0.9525159999999999
...    562|Escherichia coli                  3.5      0.12727272727272726    0.9589107142857143
...    1639|Listeria monocytogenes           3.75     0.13636363636363635    0.9627300000000001
...    1280|Staphylococcus aureus            3.875    0.1409090909090909     0.9599290322580646
...    1613|Limosilactobacillus fermentum    5.25     0.19090909090909092    0.9644999999999998
```

The output file `*.json` contains the lineage and remark of each processed reads.
```
{
    "002617ff-697a-4cd5-8a97-1e136a792228": {
        "remark": "marker-gene-containing",
        "lineage": "2|Bacteria;1239|Bacillota;91061|Bacilli;186826|Lactobacillales;81852|Enterococcaceae;1350|Enterococcus;1351|Enterococcus faecalis"
    },
    ...
    "ffe73d61-55eb-4ad8-9519-e38c364fc11d": {
        "remark": "marker-gene-containing",
        "lineage": "2|Bacteria;1224|Pseudomonadota;1236|Gammaproteobacteria;91347|Enterobacterales;543|Enterobacteriaceae;590|Salmonella;28901|Salmonella enterica"
    }
}
```

### Run Melon with pre-filtering of non-prokaryotic reads
To enable the pre-filtering module, you need to download a database of Kraken that includes at least human and fungi (PlusPF, PlusPFP, or their capped versions). Using the PlusPF-8 (ver. 2023-06-05, capped at 8 GB) as an example:

```bash
## https://benlangmead.github.io/aws-indexes/k2
mkdir database_kraken
wget -q --show-progress https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08gb_20230605.tar.gz -O database_kraken/db.tar.gz
tar -zxvf database_kraken/db.tar.gz -C database_kraken

## remove temporary files to save space
rm -rf database_kraken/db.tar.gz
```

Run with argument `-k/-db-kraken`:
```bash
melon *.fa -d database -o . -k database_kraken
```

## Miscellaneous
### Statistics of influent/effluent samples with/without pre-filter
Both samples were collected from the Shatin wastewater treatment plant, sequenced with ONT SQK-NBD114 & R10.4.1 on PromethION, basecalled with Guppy v6.5.7.

<table>
   <thead>
      <tr>
         <th></th>
         <th></th>
         <th>no pre-filter</th>
         <th>PlusPF-8</th>
         <th>PlusPF-16</th>
         <th>PlusPF</th>
      </tr>
   </thead>
   <tbody>
      <tr>
         <td rowspan="7">influent
         <div>
  <ul>
    <li>read: 1,091,772</li>
    <li>bp: 7,815,976,467</li>
    <li>n50: 9,500</li>
    <li>mean quality: 18.3</li>
  </ul>
</div></td>
         <td>genome copy</td>
         <td>1,809</td>
         <td>1,801</td>
         <td>1,800</td>
         <td>1,797</td>
      </tr>
      <tr>
         <td>species richness</td>
         <td>2,101</td>
         <td>2,098</td>
         <td>2,099</td>
         <td>2,097</td>
      </tr>
      <tr>
         <td>number of filtered reads</td>
         <td>-</td>
         <td>8,312</td>
         <td>10,212</td>
         <td>14,988</td>
      </tr>
      <tr>
         <td>mean genome size</td>
         <td>4.322</td>
         <td>4.304</td>
         <td>4.300</td>
         <td>4.292</td>
      </tr>
      <tr>
         <td>ARG abudnance</td>
         <td>0.551</td>
         <td>0.553</td>
         <td>0.553</td>
         <td>0.554</td>
      </tr>
      <tr>
         <td>real time (sec)</td>
         <td>2,056</td>
         <td>2,271</td>
         <td>2,279</td>
         <td>-</td>
      </tr>
      <tr>
         <td>peak resident set size (GB)</td>
         <td>10.629</td>
         <td>10.665</td>
         <td>17.294</td>
         <td>-</td>
      </tr>
      <tr>
         <td rowspan="7">effluent
           <ul>
    <li>read: 1,118,502</li>
    <li>bp: 5,157,996,045</li>
    <li>n50: 6,577</li>
    <li>mean quality: 18.2</li>
  </ul>
</td>
         <td>genome copy</td>
         <td>1,348</td>
         <td>1,336</td>
         <td>1,331</td>
         <td>1,315</td>
      </tr>
      <tr>
         <td>species richness</td>
         <td>1,704</td>
         <td>1,697</td>
         <td>1,700</td>
         <td>1,696</td>
      </tr>
      <tr>
         <td>number of filtered reads</td>
         <td>-</td>
         <td>16,602</td>
         <td>29,300</td>
         <td>54,774</td>
      </tr>
      <tr>
         <td>mean genome size</td>
         <td>3.826</td>
         <td>3.789</td>
         <td>3.757</td>
         <td>3.715</td>
      </tr>
      <tr>
         <td>ARG abundance</td>
         <td>0.507</td>
         <td>0.511</td>
         <td>0.512</td>
         <td>0.519</td>
      </tr>
      <tr>
         <td>real time (sec)</td>
         <td>1,341</td>
         <td>1,496</td>
         <td>1,495</td>
         <td>-</td>
      </tr>
      <tr>
         <td>peak resident set size (GB)</td>
         <td>10.893</td>
         <td>10.902</td>
         <td>17.067</td>
         <td>-</td>
      </tr>
   </tbody>
</table>

Tested on MacBook Pro 2021, Apple M1 Max, 64 GB memory, macOS Sonoma. Melon v0.1.0, NCBI database (ver. 2023-07-31) and Kraken database (ver. 2023-06-05). Mean genome size is in unit of Mb. ARG abundance is in unit of copies per cell (exlcuding multidrug ARGs). Real time and peak resident set size are measured with `time`.

## Citation
Chen, Xi, Xiaole Yin, Xianghui Shi, Weifu Yan, Yu Yang, Lei Liu, and Tong Zhang. "Melon: metagenomic long-read-based taxonomic identification and quantification using marker genes." bioRxiv (2023).
