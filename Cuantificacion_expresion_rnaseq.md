# 🧬Cargar datos
# 1.	Instalar awscli:
sudo apt-get install awscli
# 2.	Configurar awscli con las credenciales necesarias:
aws configure
# 3.	Descargar los datos del bucket de S3:
aws s3 cp s3://tgsaraparra/[directorio]/ /Datos/Sara/rawseq/[directorio]/ --recursive
# 4.	Convertir los archivos SRA a formato FASTQ si es necesario:
	  fastq-dump --split-files [archivo]

# 🧬FastQC

fastqc /Datos/Sara/rawseq/[directorio]/*.fastq.gz --outdir /Datos/Sara/fastqc/

# 🧬MultiQC

multiqc /Datos/Sara/rawfastqc -o /Datos/Sara/multiqc_results


# 🧬Limpieza

cutadapt -q 30 --max-n 0 -m 40 \
    -o "/Datos/Sara/trimmedReads/[archivo]_1P.fastq.gz" \
    -p "/Datos/Sara/trimmedReads/[archivo]_2P.fastq.gz" \
    "/Datos/Sara/rawseq/[archivo]_1.fastq.gz" \
    "/Datos/Sara/rawseq/[archivo]_2.fastq.gz"

# 🧬Alineación con STAR

STAR --runThreadN 32 \
    --genomeLoad LoadAndKeep \
    --genomeDir /Data/Data/starIndex/ \
    --chimSegmentMin 12 \
    --chimOutType Junctions \
    --readFilesCommand zcat \
    --readFilesIn "/Datos/Sara/trimmedReads/[archivo]_1P.fastq.gz" "/Datos/Sara/trimmedReads/[archivo]_2P.fastq.gz" \
    --outFileNamePrefix "/Datos/Sara/starAligned/[archivo]" \
    --limitBAMsortRAM 80000000000 \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx

# 🧬Conteo de lecturas con FeatureCounts

featureCounts -T 8 -p --countReadPairs -t exon -g gene_id \
    -a /Data/Data/ncbi_dataset/data/GCF_000001405.40/genomic.gtf \
    -o /Datos/Sara/readCounts/readCounts.txt /Datos/Sara/starAligned/*.sortedByCoord.out.bam > counts.log


# 🧬Subir los resultados a S3

aws s3 cp /Datos/Sara/readCounts/readCounts.txt s3

