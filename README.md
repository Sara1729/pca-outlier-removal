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

* **PLINK v1.9/v2** (installato direttamente o tramite Docker container [1](https://github.com/asherkhb/plink-docker) [2](https://hub.docker.com/r/miguelpmachado/plink_2.0))

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
  > /path/to/common_ids.txt
```

Nel nostro esempio, il numero di varianti comuni è **14.661**.

Ora estraiamo solo le varianti comuni da ciascun file:

```bash
bcftools view --include 'ID=@/path/to/common_ids.txt' \
  /path/to/noCHR_newID_fixref_splt_merged_file_SNP_array.vcf.gz \
  -Oz -o /path/to/common_var_SNP_array.vcf.gz
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

## 5) Identificazione e Rimozione degli Outlier

Dopo aver eseguito la PCA, un passaggio fondamentale è l’identificazione degli **outlier**.
Questi campioni possono rappresentare individui con diversa origine ancestrale, contaminazioni del campione o artefatti tecnici.
Rimuoverli è importante perché possono distorcere le componenti principali e compromettere le analisi successive.

Per identificare outlier multidimensionali è possibile utilizzare la **[distanza di Mahalanobis](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/mahalanobis)**.

Questa metrica misura la distanza di ciascuna osservazione dal centro della distribuzione dei dati (cioè il vettore delle medie), tenendo conto della struttura di covarianza tra le variabili.
In questo modo è possibile rilevare outlier considerando contemporaneamente più componenti principali (ad esempio le prime 10 PC), fornendo un criterio più robusto e affidabile rispetto agli approcci univariati.

---

### Implementazione in R

Il seguente script R può essere utilizzato per identificare gli outlier a partire dal file `.eigenvec` generato da PLINK.
Lo script esegue i seguenti passaggi:

1. Carica il file degli autovettori (`.eigenvec`).
2. Calcola la distanza di Mahalanobis per ciascun campione utilizzando le prime 10 componenti principali (PC).
   *(Si consiglia di utilizzare le prime 10 PC poiché sono quelle salvate di default da PLINK).*
3. Imposta la soglia basata sulla distribuzione del chi-quadro, specificando quantile e gradi di libertà.
   In questo caso viene utilizzato il quantile **0.999** e il numero di gradi di libertà pari al numero di PC utilizzate.
4. Segnala come outlier tutti i campioni che superano la soglia.
5. Genera un grafico PC1 vs PC2 con gli outlier evidenziati.
6. Crea un file di testo (`outliers_to_remove.txt`) contenente gli ID dei campioni da rimuovere.

---

### Prerequisiti in R

Assicurarsi di avere installate le librerie **MASS** e **ggplot2**:

```r
install.packages("MASS")
install.packages("ggplot2")
```

---

### Codice R per l’Identificazione degli Outlier

```r
# Carica le librerie necessarie
library(MASS)
library(ggplot2)

# Specifica il percorso del file .eigenvec
eigenvectors <- "merged_cohorts.eigenvec" # <-- MODIFICARE CON IL PROPRIO PERCORSO
# Nota: Plink non scrive un'intestazione, quindi deve essere creata manualmente

# --- Analisi degli Outlier ---

# Seleziona le prime 10 Componenti Principali per l'analisi
pcs_to_use <- 10
pcs <- as.matrix(eigenvectors[, paste0("PC", 1:pcs_to_use)])

# Calcola il centro (media) e la matrice di covarianza delle PC
center <- colMeans(pcs)
covmat <- cov(pcs)

# Calcola la distanza di Mahalanobis per ciascun campione
mahal <- mahalanobis(pcs, center, covmat)

# Imposta la soglia per definire un outlier
# Usiamo il quantile 0.999 della distribuzione Chi-quadro
# I gradi di libertà (df) corrispondono al numero di PC utilizzate
threshold <- qchisq(0.999, df = pcs_to_use)

# Aggiunge una colonna booleana (TRUE/FALSE) per identificare gli outlier
eigenvectors$is_outlier <- mahal > threshold

# Stampa il numero di outlier identificati nella console
cat("Numero di outlier rilevati:", sum(eigenvectors$is_outlier), "\n")

# Crea un grafico a dispersione PC1 vs PC2
# I campioni outlier verranno colorati in rosso
ggplot(eigenvectors, aes(x = PC1, y = PC2, color = is_outlier)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), name = "Outlier?") +
  theme_minimal() +
  labs(title = "Rilevamento degli Outlier con Distanza di Mahalanobis",
       subtitle = paste(sum(eigenvectors$is_outlier), "outlier rilevati"),
       x = "Componente Principale 1 (PC1)",
       y = "Componente Principale 2 (PC2)") +
  coord_fixed()

# Estrae gli ID dei campioni outlier
# Il flag --remove di Plink richiede un file a due colonne (FID e IID) senza intestazione
outliers_to_remove <- eigenvectors[eigenvectors$is_outlier, c("ID")]

# Scrive il file che Plink utilizzerà per la rimozione
output_file <- "outliers_to_remove.txt"
write.table(outliers_to_remove,
            file = output_file,
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

cat("File '", output_file, "' creato con", nrow(outliers_to_remove), "campioni da rimuovere.\n")
```


