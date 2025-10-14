# PCA da file VCF

Questo progetto illustra come eseguire un'**Analisi delle Componenti Principali (PCA)** partendo da uno o più file VCF (aggregati o non aggregati).

## Introduzione

Questo documento descrive il processo per effettuare una PCA a partire da file VCF.
È **necessario** che i file VCF contengano il campo `GT` (genotipo), poiché gli altri campi non vengono considerati in questa analisi.

La PCA è uno strumento di clustering semplice ma potente, basato sulla decomposizione della matrice dei dati.
Nel nostro caso, le **righe** rappresentano le varianti e le **colonne** rappresentano i genotipi dei campioni.
L’obiettivo è trovare autovalori (eigenvalues) e autovettori (eigenvectors).

Gli autovalori vengono ordinati in modo decrescente: il primo autovalore identifica la prima componente principale (PC1) e la sua direzione nello spazio delle componenti principali è determinata dall’autovettore corrispondente.
In sintesi, la PCA consente di individuare uno spazio a dimensionalità ridotta che descrive al meglio la distribuzione dei dati.

---

## Prerequisiti e Strumenti necessari

Per eseguire questa pipeline sono richiesti:

* **PLINK v1.9** (installato direttamente o tramite [Docker container](https://github.com/asherkhb/plink-docker))

* **BCFtools** (installato direttamente o tramite [Docker container](https://github.com/samtools/bcftools))

* **File di conversione dei nomi dei cromosomi:**

  1. File `no_chr_name_convention.txt`, con due colonne:

     ```
     chr1     1
     chr2     2
     ...
     chr22    22
     chrX     X
     ```
  2. File `chr_name_conv.txt`, con le colonne invertite:

     ```
     1        chr1
     2        chr2
     ...
     22       chr22
     23       chrX
     ```

* **Genoma di riferimento** in formato `.fa` e relativo indice `.fai` (ad es. hg19 o hg38).
  Disponibile da [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/).

---

## Nota su PLINK nel nostro caso

`Plink` accetta in input file `VCF.gz` accompagnati dal relativo file di indice `.tbi`.
È fondamentale che i file contengano il campo `GT` (genotipo).

Nel nostro caso studio, la colonna `CHROM` deve essere in formato `chrN` e la colonna `ID` deve essere costruita come:

```
CHROM:POS:REF:ALT
```

Questo formato è essenziale per eseguire correttamente l’intersezione tra diverse coorti.

---

## Descrizione delle Coorti

* **SNP_array**

  * Build: **hg19**
  * Numero di campioni: **555**
  * Numero di varianti selezionate: **14.861**

* **your cohort**

  * Seconda coorte da confrontare con SNP_array.

---

## 1) Normalizzazione dei file VCF

In questa fase vogliamo assicurarci che:

* la colonna `CHROM` inizi con il prefisso `chr`;
* la colonna `ID` sia nel formato `CHROM:POS:REF:ALT`.

Esempio di comando per creare correttamente la colonna `ID`:

```bash
bcftools annotate --rename-chrs /path/to/chr_name_conv.txt \
  /path/to/fixref_splt_merged_file_your_cohort.vcf.gz | \
  bcftools norm -Ou -f /path/to/hg19.fa | \
  bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' \
  -o /path/to/newID_fixref_splt_merged_file_your_cohort.vcf.gz
```

E per garantire che la colonna `CHROM` abbia il prefisso `chr`:

```bash
bcftools annotate --rename-chrs /path/to/no_chr_name_convention.txt \
  /path/to/newID_fixref_splt_merged_file_your_cohort.vcf.gz \
  -Oz -o /path/to/noCHR_newID_fixref_splt_merged_file_your_cohort.vcf.gz
```

---

## 2) Unione delle Coorti

Una volta che entrambi i file VCF sono normalizzati (`SNP_array` e `your cohort`), possiamo cercare le varianti comuni.
Le coorti verranno indicate come:

```
/path/to/noCHR_newID_fixref_splt_merged_file_SNP_array.vcf.gz \
/path/to/noCHR_newID_fixref_splt_merged_file_your_cohort.vcf.gz
```

### Intersezione delle varianti comuni

```bash
bcftools isec -n=2 -w1 \
  /path/to/noCHR_newID_fixref_splt_merged_file_SNP_array.vcf.gz \
  /path/to/noCHR_newID_fixref_splt_merged_file_your_cohort.vcf.gz \
  -Oz | bcftools query -f '%ID\n' \
  > /path/to/common_var/common_ids.txt
```

Nel nostro esempio, il numero di varianti comuni è **14.661**.

Ora estraiamo solo le varianti comuni da ciascun file:

```bash
bcftools view --include 'ID=@/path/to/common_var/common_ids.txt' \
  /path/to/noCHR_newID_fixref_splt_merged_file_SNP_array.vcf.gz \
  -Oz -o /path/to/common_var/common_var_SNP_array.vcf.gz
```

```bash
bcftools view --include 'ID=@/path/to/common_var/common_ids.txt' \
  /path/to/noCHR_newID_fixref_splt_merged_file_your_cohort.vcf.gz \
  -Oz -o /path/to/common_var/common_var_your_cohort.vcf.gz
```

Infine, uniamo le due coorti:

```bash
bcftools merge --threads 64 \
  /path/to/common_var/common_var_SNP_array.vcf.gz \
  /path/to/common_var/common_var_your_cohort.vcf.gz \
  -Oz -o /path/to/common_var/commonVar_SNParray_yourCohort_hg19.vcf.gz
```

---

## 3) Esecuzione della PCA con PLINK 1.9

Una volta ottenuto il file VCF finale unificato (`commonVar_SNParray_yourCohort_hg19.vcf.gz`), possiamo eseguire la PCA:

```bash
plink --vcf /path/to/commonVar_SNParray_yourCohort_hg19.vcf.gz \
      --pca --double-id --out commonVar_SNParray_yourCohort_hg19
```

### Opzioni principali:

* `--vcf`: specifica il file di input;
* `--pca`: esegue la PCA;
* `--double-id`: mantiene Family ID e Sample ID identici;
* `--out`: prefisso dei file di output (`.eigenvec`, `.eigenval`, `.log`).

---

## 4) Esecuzione della PCA con PLINK 2

```bash
plink2 --vcf /path/to/commonVar_SNParray_yourCohort_hg19.vcf.gz \
       --make-pgen --out /path/to/pca/commonVar_SNParray_yourCohort_hg19
```

Poi:

```bash
plink2 --pfile /path/to/pca/commonVar_SNParray_yourCohort_hg19 \
       --freq counts \
       --pca allele-wts vcols=chrom,ref,alt \
       --out /path/to/pca/pca_commonVar_SNParray_yourCohort_hg19
```

