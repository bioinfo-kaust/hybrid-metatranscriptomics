/*
 * Stage 0 — reference acquisition and indexing.
 * Nothing here is species-specific: every input is a parameter, and anything that can be
 * fetched from a public source is fetched rather than assumed to be staged already.
 */

process DOWNLOAD_GENOME {
    label 'process_download'
    tag "$accession"

    input:
    val accession
    val name

    output:
    path "${name}.fna", emit: fasta
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    datasets download genome accession ${accession} --include genome --filename ${accession}.zip
    unzip -o -q ${accession}.zip -d ${accession}_dir
    cat ${accession}_dir/ncbi_dataset/data/${accession}/*_genomic.fna > ${name}.fna
    rm -rf ${accession}.zip ${accession}_dir
    grep -c '^>' ${name}.fna | sed 's/^/[${name}] scaffolds: /'
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    touch ${name}.fna
    emit_versions.sh "${task.process}"
    """
}

/*
 * Contaminant genomes to subtract. Accepts a comma-separated accession list; the genomes are
 * downloaded, concatenated and their scaffolds prefixed so a contaminant hit is a prefix test.
 */
process DOWNLOAD_CONTAMINANT_GENOMES {
    label 'process_download'
    tag "${accessions}"

    input:
    val accessions

    output:
    path 'contaminants_combined.fna', emit: fasta
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    : > contaminants_combined.fna
    for acc in \$(echo "${accessions}" | tr ',' ' '); do
        echo "[*] datasets download \$acc"
        datasets download genome accession "\$acc" --include genome --filename "\$acc.zip"
        unzip -o -q "\$acc.zip" -d "\$acc"_dir
        cat "\$acc"_dir/ncbi_dataset/data/"\$acc"/*_genomic.fna >> contaminants_combined.fna
        rm -rf "\$acc.zip" "\$acc"_dir
    done
    grep -c '^>' contaminants_combined.fna | sed 's/^/[contaminants] scaffolds: /'
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    touch contaminants_combined.fna
    emit_versions.sh "${task.process}"
    """
}

/* A directory of FASTA files, combined as-is. */
process COMBINE_CONTAMINANTS {
    label 'process_low'
    tag "${fasta_dir}"

    input:
    path fasta_dir

    output:
    path 'contaminants_combined.fna', emit: fasta
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    shopt -s nullglob
    files=( ${fasta_dir}/*.fa ${fasta_dir}/*.fna ${fasta_dir}/*.fasta ${fasta_dir}/*.fa.gz ${fasta_dir}/*.fna.gz ${fasta_dir}/*.fasta.gz )
    [ \${#files[@]} -gt 0 ] || { echo "[ERROR] no FASTA files in ${fasta_dir}" >&2; exit 1; }
    : > contaminants_combined.fna
    for f in "\${files[@]}"; do
        case "\$f" in *.gz) zcat "\$f" ;; *) cat "\$f" ;; esac >> contaminants_combined.fna
    done
    grep -c '^>' contaminants_combined.fna | sed 's/^/[contaminants] scaffolds: /'
    emit_versions.sh "${task.process}" seqkit
    """

    stub:
    """
    touch contaminants_combined.fna
    emit_versions.sh "${task.process}" seqkit
    """
}

/*
 * rRNA reference for k-mer depletion. Default set = SortMeRNA v4.3 (SILVA 16S/18S/23S/28S +
 * Rfam 5S/5.8S); any other set of URLs works via --rrna_db_sets / --rrna_db_base.
 */
process DOWNLOAD_RRNA_DB {
    label 'process_download'

    output:
    path 'rrna_db.fasta', emit: fasta
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    : > rrna_db.fasta
    for f in \$(echo "${params.rrna_db_sets}" | tr ',' ' '); do
        echo "[*] \$f"
        curl -sL --max-time 600 "${params.rrna_db_base}/\$f.fasta" >> rrna_db.fasta
    done
    [ -s rrna_db.fasta ] || { echo "[ERROR] rRNA database download produced an empty file" >&2; exit 1; }
    grep -c '^>' rrna_db.fasta | sed 's/^/[rRNA db] sequences: /'
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    touch rrna_db.fasta
    emit_versions.sh "${task.process}"
    """
}

/* A protein FASTA for DIAMOND (Swiss-Prot, a specialist set, or anything else). */
process DOWNLOAD_PROTEIN_DB {
    label 'process_download'
    tag "$name"

    input:
    tuple val(name), val(url)

    output:
    tuple val(name), path("${name}.fasta.gz"), emit: fasta
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    curl -sL --max-time 1800 "${url}" -o ${name}.fasta.gz
    # UniProt streams can arrive uncompressed despite compressed=true
    if ! gzip -t ${name}.fasta.gz 2>/dev/null; then
        mv ${name}.fasta.gz ${name}.fasta.tmp && gzip -c ${name}.fasta.tmp > ${name}.fasta.gz && rm ${name}.fasta.tmp
    fi
    [ -s ${name}.fasta.gz ] || { echo "[ERROR] empty download for ${name}" >&2; exit 1; }
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    touch ${name}.fasta.gz
    emit_versions.sh "${task.process}"
    """
}

process MAKE_DIAMOND_DB {
    label 'process_medium'
    tag "$name"

    input:
    tuple val(name), path(fasta)

    output:
    tuple val(name), path("${name}.dmnd"), emit: db
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    diamond makedb --in ${fasta} -d ${name} --threads ${task.cpus} --quiet
    emit_versions.sh "${task.process}" diamond
    """

    stub:
    """
    touch ${name}.dmnd
    emit_versions.sh "${task.process}" diamond
    """
}

process STAR_INDEX {
    label 'process_high'
    tag "${genome.name}"

    input:
    path genome

    output:
    path 'star_index', emit: index
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p star_index
    # No GTF: junctions are discovered de novo by the 2-pass alignment. Supplying an unreliable
    # annotation for a draft genome would bias the discovery this pipeline exists to perform.
    STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ${genome} \\
         --runThreadN ${task.cpus} ${args}
    emit_versions.sh "${task.process}" STAR
    """

    stub:
    """
    mkdir -p star_index && touch star_index/SA
    emit_versions.sh "${task.process}" STAR
    """
}

process SEQ_STATS {
    label 'process_single'
    tag "$label"

    input:
    tuple val(label), path(fasta)

    output:
    tuple val(label), path("${label}.stats.tsv"), emit: stats
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    seqkit stats -T ${fasta} > ${label}.stats.tsv
    cat ${label}.stats.tsv
    emit_versions.sh "${task.process}" seqkit
    """

    stub:
    """
    printf "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tN50\n${label}\tFASTA\tDNA\t100\t100000\t300\t1000\t5000\t1200\n" > ${label}.stats.tsv
    emit_versions.sh "${task.process}" seqkit
    """
}
