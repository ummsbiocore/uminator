$HOSTNAME = ""
params.outdir = 'results'  


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 

if (params.reads){
Channel
	.fromFilePairs( params.reads,checkExists:true , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 ) 
	.set{g_2_0_g_1}
  } else {  
	g_2_0_g_1 = Channel.empty()
 }

Channel.value(params.mate).set{g_5_1_g_1}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* autofill

process readsFiltering {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${sample}$/) "readsFiltering/$filename"}
input:
 tuple val(sample), file(reads)
 val mate

output:
 tuple val(sample), file("${sample}")  ,emit:g_1_outputDir00_g_4 

script:
if(params.readsFiltering == "yes")
"""
mkdir -p ${sample}
#filter reads
zless ${reads} | NanoFilt -q ${params.minQ} --length ${params.minLen} --maxlength ${params.maxLen} > ${sample}/${sample}_filtered.fastq
"""
else
"""
mkdir -p ${sample}
zless ${reads} > ${sample}/${sample}_filtered.fastq
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* autofill

process candidateUMIsExtraction {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${sample}$/) "candidateUMIsExtraction/$filename"}
input:
 tuple val(sample), file(dir)

output:
 tuple val(sample), file("${sample}")  ,emit:g_4_outputDir00_g_6 

script:
"""
mkdir readsFiltering && mv $dir readsFiltering/.
    mkdir -p candidateUMIsExtraction/${sample}

    #obtain reverse complement of primers and adapters sequences
    FW_primer_R=\$(echo -e \">tmp\\n\" ${params.FW_primer} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')
    RV_primer_R=\$(echo -e \">tmp\\n\" ${params.RV_primer} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')
    FW_adapter_R=\$(echo -e \">tmp\\n\" ${params.FW_adapter} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')
    RV_adapter_R=\$(echo -e \">tmp\\n\" ${params.RV_adapter} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')

    #double UMI design
    #obtain first and last bases of reads
    READS_START=candidateUMIsExtraction/${sample}/first_${params.searchLen}bp.fastq
    READS_END=candidateUMIsExtraction/${sample}/last_${params.searchLen}bp.fastq
    seqtk trimfq -L ${params.searchLen} readsFiltering/${sample}/${sample}_filtered.fastq > \$READS_START
    seqtk seq -r readsFiltering/${sample}/${sample}_filtered.fastq | seqtk trimfq -L ${params.searchLen} - | seqtk seq -r - > \$READS_END
    
    #search candidate UMIs with exact length between adapters and primers
    cutadapt -j ${task.cpus} -e ${params.tolCutadaptErr} -O ${params.minLenOvlp} -m ${params.UMILen} -M ${params.UMILen} \
    --discard-untrimmed --match-read-wildcards \
    -g ${params.FW_adapter}...${params.FW_primer} -g ${params.RV_adapter}...${params.RV_primer} \
    -G \$RV_primer_R...\$RV_adapter_R -G \$FW_primer_R...\$FW_adapter_R \
    -o candidateUMIsExtraction/${sample}/UMI_part1_db_tmp1.fastq \
    -p candidateUMIsExtraction/${sample}/UMI_part2_db_tmp1.fastq \
    \$READS_START \$READS_END

    #collapse the 5' and 3' end UMIs of each read
    paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' candidateUMIsExtraction/${sample}/UMI_part1_db_tmp1.fastq ) \
    <( sed -n '1~4s/^@/>/p;2~4p' candidateUMIsExtraction/${sample}/UMI_part2_db_tmp1.fastq ) | cut -d " " -f1 \
    > candidateUMIsExtraction/${sample}/UMI_db_tmp1.fasta

    #search candidate UMIs with approximate length between adapters and primers
    cutadapt -j ${task.cpus} -e ${params.tolCutadaptErr} -O ${params.minLenOvlp} -m ${params.UMILen} -l ${params.UMILen} \
    --discard-untrimmed --match-read-wildcards \
    -g ${params.FW_adapter} -g ${params.RV_adapter} \
    -G \$RV_primer_R -G \$FW_primer_R \
    -o candidateUMIsExtraction/${sample}/UMI_part1_candidates.fastq \
    -p candidateUMIsExtraction/${sample}/UMI_part2_candidates.fastq \
    \$READS_START \$READS_END

    #collapse the 5' and 3' end UMIs of each read
    paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' candidateUMIsExtraction/${sample}/UMI_part1_candidates.fastq ) \
    <( sed -n '1~4s/^@/>/p;2~4p' candidateUMIsExtraction/${sample}/UMI_part2_candidates.fastq ) |  \
    cut -d " " -f1 > candidateUMIsExtraction/${sample}/UMI_candidates.fastq
      
    #convert UMI candidates to fasta
    seqtk seq -A candidateUMIsExtraction/${sample}/UMI_candidates.fastq > candidateUMIsExtraction/${sample}/UMI_candidates.fasta
    mv candidateUMIsExtraction/${sample} .
    """
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* autofill

process candidateUMIsFiltering {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /$sample$/) "candidateUMIsFiltering/$filename"}
input:
 tuple val(sample), file(dir)

output:
 tuple val(sample), file("$sample")  ,emit:g_6_outputDir00_g_9 

script:
"""
mkdir candidateUMIsExtraction && mv $dir candidateUMIsExtraction/.
mkdir -p candidateUMIsFiltering
mkdir -p candidateUMIsFiltering/${sample}
    
    #filter candidate UMIs for pattern
    cat candidateUMIsExtraction/${sample}/UMI_db_tmp1.fasta | grep -B1 -E ${params.UMIPattern} | sed \'/^--\$/d\' \
    > candidateUMIsFiltering/${sample}/UMI_db_tmp2.fasta

    #evaluate the total UMI length (UMI1 + UMI2)
    totUMILen=\$(echo 2*${params.UMILen} | bc)
    
    #dereplicate high-quality UMIs
    vsearch --derep_fulllength candidateUMIsFiltering/${sample}/UMI_db_tmp2.fasta \
    --minseqlength \$totUMILen --sizeout --relabel umi --strand both \
    --output candidateUMIsFiltering/${sample}/UMI_db_tmp3.fasta

    #cluster dereplicated high-quality UMIs to obtain a database of high-quality UMIs
    vsearch --cluster_smallmem candidateUMIsFiltering/${sample}/UMI_db_tmp3.fasta \
    --usersort --id ${params.UMIClustID} --iddef 1 --strand both --clusterout_sort --minseqlength \$totUMILen \
    --sizein --sizeout --consout candidateUMIsFiltering/${sample}/UMI_db_tmp4.fasta \
    --minwordmatches 0

    #index the database of tmp high-quality UMIs
    bwa index candidateUMIsFiltering/${sample}/UMI_db_tmp4.fasta

    #align candidate UMIs to tmp high-quality UMIs
    bwa aln \
    candidateUMIsFiltering/${sample}/UMI_db_tmp4.fasta \
    candidateUMIsExtraction/${sample}/UMI_candidates.fasta \
    -n ${params.maxDiff} \
    -t ${task.cpus} \
    -N > candidateUMIsFiltering/${sample}/UMI_candidates_map_tmp1.sai

    bwa samse \
    -n 10000000 \
    candidateUMIsFiltering/${sample}/UMI_db_tmp4.fasta \
    candidateUMIsFiltering/${sample}/UMI_candidates_map_tmp1.sai \
    candidateUMIsExtraction/${sample}/UMI_candidates.fasta | \
    samtools view -F 4 - \
    > candidateUMIsFiltering/${sample}/UMI_candidates_map_tmp1.sam

    #retain only UMIs from the database that are supported by at least 3 candidate UMIS
    cat candidateUMIsFiltering/${sample}/UMI_candidates_map_tmp1.sam \
    | cut -f3 | sort | uniq -c | sort -nr | awk \'BEGIN {FS=\" \"} { if ( \$1 > 2 ) print "\t"\$2"\t"}\' \
    > candidateUMIsFiltering/${sample}/UMI_candidates_map_supp3.txt

    #extract from fasta only UMIs supported by at least 3 candidate UMIs
    cat candidateUMIsFiltering/${sample}/UMI_candidates_map_supp3.txt \
    | cut -f2 > candidateUMIsFiltering/${sample}/UMI_candidates_map_supp3_noTab.txt

    seqtk subseq candidateUMIsFiltering/${sample}/UMI_db_tmp4.fasta \
    candidateUMIsFiltering/${sample}/UMI_candidates_map_supp3_noTab.txt \
    > candidateUMIsFiltering/${sample}/UMI_db.fasta

    #split the database of high-quality UMIs in part 1 and part 2
    #trim UMILen from right
    seqtk trimfq -e ${params.UMILen} candidateUMIsFiltering/${sample}/UMI_db.fasta \
    > candidateUMIsFiltering/${sample}/UMI_db_p1_fw.fasta

    seqtk seq -r candidateUMIsFiltering/${sample}/UMI_db_p1_fw.fasta \
    | sed \'/^>/ s/\$/_rc/\' > candidateUMIsFiltering/${sample}/UMI_db_p1_rv.fasta

    #trim UMILen from left
    seqtk trimfq -b ${params.UMILen} candidateUMIsFiltering/${sample}/UMI_db.fasta \
    > candidateUMIsFiltering/${sample}/UMI_db_p2_fw.fasta

    seqtk seq -r candidateUMIsFiltering/${sample}/UMI_db_p2_fw.fasta \
    | sed \'/^>/ s/\$/_rc/\' > candidateUMIsFiltering/${sample}/UMI_db_p2_rv.fasta
    
    #concatenate UMI1 and reverse complement of UMI2
    cat candidateUMIsFiltering/${sample}/UMI_db_p1_fw.fasta \
    candidateUMIsFiltering/${sample}/UMI_db_p2_rv.fasta \
    > candidateUMIsFiltering/${sample}/UMI_db_p1.fasta

    #concatenate UMI2 and reverse complement of UMI1
    cat candidateUMIsFiltering/${sample}/UMI_db_p2_fw.fasta \
    candidateUMIsFiltering/${sample}/UMI_db_p1_rv.fasta \
    > candidateUMIsFiltering/${sample}/UMI_db_p2.fasta

    mv candidateUMIsFiltering/${sample} .
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 30
}
//* autofill

process readsUMIsAssignment {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${sample}_readsUMIsAssignment$/) "readsUMIsAssignment/$filename"}
input:
 tuple val(sample), file(dir)
 path "candidateUMIsExtraction/*"
 path "readsFiltering/*"

output:
 tuple val(sample), env(UMI_all)  ,emit:g_9_umi00_g_11 
 path "${sample}_readsUMIsAssignment"  ,emit:g_9_outputDir11_g_11 

stageInMode 'copy'

script:
"""
mkdir candidateUMIsFiltering && mv $dir candidateUMIsFiltering/.
mkdir -p readsUMIsAssignment
mkdir -p readsUMIsAssignment/${sample}
path=\$(which Bin_reads.R) && cp \$path .
path=\$(which Filter_UMIs.R) && cp \$path .


    #align high-quality UMIs (part 1 and part 2) to the terminal portion of reads
    READS_START_FQ=candidateUMIsExtraction/${sample}/first_${params.searchLen}bp.fastq
    READS_START_FA=candidateUMIsExtraction/${sample}/first_${params.searchLen}bp.fasta
    READS_END_FQ=candidateUMIsExtraction/${sample}/last_${params.searchLen}bp.fastq
    READS_END_FA=candidateUMIsExtraction/${sample}/last_${params.searchLen}bp.fasta
    
    #convert fastq to fasta
    seqtk seq -A \$READS_START_FQ > \$READS_START_FA
    seqtk seq -A \$READS_END_FQ > \$READS_END_FA

    #index db
    bwa index \$READS_START_FA
    bwa index \$READS_END_FA

    maxDiffSingle=\$(echo ${params.maxDiff}/2 | bc)
    
    #map UMI_db_p1 to reads start
    bwa aln \
    \$READS_START_FA \
    candidateUMIsFiltering/${sample}/UMI_db_p1.fasta \
    -n \$maxDiffSingle \
    -t ${task.cpus} \
    -N > readsUMIsAssignment/${sample}/UMI_db_p1.sai

    bwa samse \
    -n 10000000 \
    \$READS_START_FA \
    readsUMIsAssignment/${sample}/UMI_db_p1.sai \
    candidateUMIsFiltering/${sample}/UMI_db_p1.fasta | \
    samtools view -F 20 - \
    > readsUMIsAssignment/${sample}/UMI_db_p1.sam

    #map UMI_db_p2 to reads end
    bwa aln \
    \$READS_END_FA \
    candidateUMIsFiltering/${sample}/UMI_db_p2.fasta \
    -n \$maxDiffSingle \
    -t ${task.cpus} \
    -N > readsUMIsAssignment/${sample}/UMI_db_p2.sai

    bwa samse \
    -n 10000000 \
    \$READS_END_FA \
    readsUMIsAssignment/${sample}/UMI_db_p2.sai \
    candidateUMIsFiltering/${sample}/UMI_db_p2.fasta | \
    samtools view -F 20 - \
    > readsUMIsAssignment/${sample}/UMI_db_p2.sam

    #filter alignments
    /opt/conda/envs/UMInator_env/bin/Rscript Filter_UMIs.R alignment_file_1=readsUMIsAssignment/${sample}/UMI_db_p1.sam alignment_file_2=readsUMIsAssignment/${sample}/UMI_db_p2.sam map_file=readsUMIsAssignment/${sample}/UMI_read_map.txt max_NM_mean=${params.max_NM_mean} max_NM_sd=${params.max_NM_sd}

    #read the current fastq file and, if a read matches a single UMI, assign it to it
    /opt/conda/envs/UMInator_env/bin/Rscript Bin_reads.R fastq_file=readsFiltering/${sample}/${sample}_filtered.fastq map_file=readsUMIsAssignment/${sample}/UMI_read_map.txt outdir=readsUMIsAssignment/${sample}

    #extract all UMIs
    UMI_all_tmp=\$(cat readsUMIsAssignment/${sample}/UMI_read_map.txt | cut -f2 | sed \'s/^centroid=//\' | sed \'s/;seqs=.*//\' | sort | uniq)

    #add sample name before UMI ID
    sn=${sample}
    UMI_all=\$(echo \$UMI_all_tmp | sed \"s/umi/\$sn\"\\|umi\"/g\")
    mv readsUMIsAssignment/${sample} ${sample}_readsUMIsAssignment
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 30
}
//* autofill

process draftConsensusCalling {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${sample}_${UMI}$/) "draftConsensusCalling/$filename"}
input:
 tuple val(sample), val(UMI)
 path "readsUMIsAssignment/*"

output:
 tuple val(sample), val(UMI)  ,emit:g_11_umi00_g_12 
 path "${sample}_${UMI}"  ,emit:g_11_outputDir11 
 tuple val(sample), val(UMI), file("${sample}_${UMI}")  ,emit:g_11_umi20_g_14 

stageInMode 'copy'

script:
"""
mkdir -p draftConsensusCalling
mkdir -p draftConsensusCalling/${sample}
path=\$(which Obtain_draft_consensus.R) && cp \$path .
mv "readsUMIsAssignment/${sample}_readsUMIsAssignment" readsUMIsAssignment/${sample}
find readsUMIsAssignment/${sample}/

#concatenate files with the same UMI obtained from different reads chunks

chunks_binned_files_curr_UMI=\$(find readsUMIsAssignment/${sample}/ | grep ${UMI}_chunk);
echo "chunks_binned_files_curr_UMI: \$chunks_binned_files_curr_UMI"
cat \$chunks_binned_files_curr_UMI > readsUMIsAssignment/${sample}/${UMI}.fastq
rm \$chunks_binned_files_curr_UMI

#convert fastq to fasta
echo "seqtk seq -A readsUMIsAssignment/${sample}/${UMI}.fastq > readsUMIsAssignment/${sample}/${UMI}.fasta"
seqtk seq -A readsUMIsAssignment/${sample}/${UMI}.fastq > readsUMIsAssignment/${sample}/${UMI}.fasta

#obtain draft consensus sequence
echo "ready to Obtain_draft_consensus.R"
ls readsUMIsAssignment/${sample}/${UMI}.fastq
find .
/opt/conda/envs/UMInator_env/bin/Rscript Obtain_draft_consensus.R fastq_file=readsUMIsAssignment/${sample}/${UMI}.fastq TRC=${params.target_reads_consensus} PLUR=${params.plurality} num_threads=${task.cpus} fast_alignment_flag=${params.fast_alignment_flag} fast_consensus_flag=${params.fast_consensus_flag} min_UMI_freq=${params.min_UMI_freq} 
mv draftConsensusCalling/${sample} ${sample}_${UMI}
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 30
}
//* autofill

process QC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /$sample$/) "QC/$filename"}
input:
 val sample
 path "readsUMIsAssignment/*"

output:
 path "$sample"  ,emit:g_12_outputDir00 

stageInMode 'copy'

when:
params.QC == "yes"

script:
"""
    mkdir -p QC
    mkdir -p QC/${sample}
    mv "readsUMIsAssignment/${sample}_readsUMIsAssignment" readsUMIsAssignment/${sample}
    find readsUMIsAssignment/
    
    #concatenate files with the same UMI obtained from different reads chunks
    chunks_unbinned_files=\$(find readsUMIsAssignment/${sample}/ -name \"*unbinned_chunk*\");
    unbinned_reads_files=readsUMIsAssignment/${sample}/unbinned.fastq
    binned_reads_files=readsUMIsAssignment/${sample}/binned.fastq
    if [[ -f "\$chunks_unbinned_files" ]]; then cat \$chunks_unbinned_files > \$unbinned_reads_files; rm \$chunks_unbinned_files; fi
    fastq_files_binned=\$(find readsUMIsAssignment/${sample} -name \"umi*\\.fastq\")
    fastq_files=\$(find readsUMIsAssignment/${sample} -name \"*\\.fastq\")
    if [[ ! -z "\$fastq_files_binned" ]]; then
      #cat \$fastq_files_binned > \$binned_reads_files
      for f in \$fastq_files_binned; do
        cat \$f >> \$binned_reads_files;
      done
    fi
    
    #do QC plot for unbinned reads
    if [[ -f "\$unbinned_reads_files" ]]; then
      NanoPlot -t ${task.cpus} --fastq \$unbinned_reads_files -o QC/${sample}/QC_unbinned_reads
    fi

    #do QC plot for binned reads
    if [[ -f "\$binned_reads_files" ]]; then
      NanoPlot -t ${task.cpus} --fastq \$binned_reads_files -o QC/${sample}/QC_binned_reads 
      rm \$binned_reads_files
    fi
    
    #produce tsv files with read-UMI assignment stats
    for f in \$fastq_files; do
      reads_names=\$(seqtk seq -A \$f | grep \"^>\" | sed \'s/>//\' | paste -sd ",")
      num_reads=\$(seqtk seq -A \$f | grep \"^>\" | sed \'s/>//\' | wc -l)
      echo -e \"\$(basename \$f)\t\$num_reads\t\$reads_names\" >> QC/${sample}/${sample}_UMI_stats_tmp.tsv
    done
    
    cat QC/${sample}/${sample}_UMI_stats_tmp.tsv | sort -k2,2nr > QC/${sample}/${sample}_UMI_stats.tsv

    rm QC/${sample}/${sample}_UMI_stats_tmp.tsv
    mv QC/${sample} .
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 20
}
//* autofill

process consensusPolishing {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${sample}$/) "consensusPolishing/$filename"}
input:
 tuple val(sample), val(UMI), file(draftConsensusCalling)
 path "readsUMIsAssignment/*"

output:
 tuple val(sample), file("${sample}/${UMI}")  ,emit:g_14_umi00_g_15 
 path "${sample}"  ,emit:g_14_outputDir11 

stageInMode 'copy'

script:
  if(params.consensusPolishing == "yes")
  """
  find .
  mkdir draftConsensusCalling
  cp -R ${draftConsensusCalling} draftConsensusCalling/${sample}
  echo "consensusPolishing started "
  mkdir -p consensusPolishing/${sample}
  path=\$(which Polish_consensus.R) && cp \$path .

  mv "readsUMIsAssignment/${sample}_readsUMIsAssignment" readsUMIsAssignment/${sample}
    
    #polish consensus sequence with racon and medaka
    /opt/conda/envs/UMInator_env/bin/Rscript Polish_consensus.R draft_consensus=draftConsensusCalling/${sample}/${UMI}/${UMI}_draft_consensus.fasta fastq_file=readsUMIsAssignment/${sample}/${UMI}.fastq TRP=${params.target_reads_polishing}  num_threads=${task.cpus} fast_polishing_flag=${params.fast_polishing_flag} medaka_model=${params.medaka_model}
    mv consensusPolishing/${sample} ${sample}
    find .
  """
  else
  """
  find .
  echo "INFO: files loaded"
  mkdir draftConsensusCalling
  mkdir -p consensusPolishing/${sample}/${UMI}
  cp -R ${draftConsensusCalling} draftConsensusCalling/${sample}
  echo "INFO: files ready"
  find .
  echo "consensusPolishing Skipped "
  cp draftConsensusCalling/${sample}/${UMI}/${UMI}_draft_consensus.fasta consensusPolishing/${sample}/${UMI}/${UMI}_polished_consensus.fasta
  mv consensusPolishing/${sample} ${sample}
  find .
  """
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 20
}
//* autofill

process primersTrimming {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /$sample$/) "primersTrimming/$filename"}
input:
 tuple val(sample), path(consDir)

output:
 path "$sample"  ,emit:g_15_outputDir00 

stageInMode 'copy'

script:


  if(params.primersTrimming == "yes")
  """
   find .
    mkdir -p primersTrimming/${sample} consensusPolishing/$sample
    mv $consDir consensusPolishing/$sample

    #find consensus sequences for all UMIs of one sample and concatenate them
    echo "files ready"
    find .
    polished_consensus=\$(find consensusPolishing/${sample} | grep "_polished_consensus.fasta")
    echo "polished_consensus files: \$polished_consensus"
    cat \$polished_consensus > primersTrimming/${sample}/${sample}_consensus_polished.fasta

    #obtain reverse complement of primers and adapters sequences
    FW_primer_R=\$(echo -e \">tmp\\n\" ${params.FW_primer} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')
    RV_primer_R=\$(echo -e \">tmp\\n\" ${params.RV_primer} | seqtk seq -r - | grep -v \"^>\" | tr -d ' ')
  
    #trim PCR primers and external sequence
    if [[ -f "primersTrimming/${sample}/${sample}_consensus_polished.fasta" ]]; then
      cutadapt -j ${task.cpus} -e 0 \
      --discard-untrimmed --match-read-wildcards \
      -g ${params.FW_primer} -g ${params.RV_primer}  \
      -a \$FW_primer_R -a \$RV_primer_R \
      -m ${params.minLen} -M ${params.maxLen} \
      -o primersTrimming/${sample}/${sample}_consensus_polished_primersTrimmed.fasta \
      primersTrimming/${sample}/${sample}_consensus_polished.fasta
    fi
    mv primersTrimming/${sample} .
  """
  else
  """
    mkdir -p primersTrimming/${sample} consensusPolishing/$sample
    mv $consDir consensusPolishing/$sample
    
    #find consensus sequences for all UMIs of one sample and concatenate them
    polished_consensus=\$(find consensusPolishing/${sample} | grep "_polished_consensus.fasta")
    cat \$polished_consensus > primersTrimming/${sample}/${sample}_consensus_polished.fasta
    mv primersTrimming/${sample} .
  """
}


workflow {


readsFiltering(g_2_0_g_1,g_5_1_g_1)
g_1_outputDir00_g_4 = readsFiltering.out.g_1_outputDir00_g_4
(g_1_outputDir02_g_9) = [g_1_outputDir00_g_4]


candidateUMIsExtraction(g_1_outputDir00_g_4)
g_4_outputDir00_g_6 = candidateUMIsExtraction.out.g_4_outputDir00_g_6
(g_4_outputDir01_g_9) = [g_4_outputDir00_g_6]


candidateUMIsFiltering(g_4_outputDir00_g_6)
g_6_outputDir00_g_9 = candidateUMIsFiltering.out.g_6_outputDir00_g_9


readsUMIsAssignment(g_6_outputDir00_g_9,g_4_outputDir01_g_9.map{ file -> return file[1] }.flatten().toList(),g_1_outputDir02_g_9.map{ file -> return file[1] }.flatten().toList())
g_9_umi00_g_11 = readsUMIsAssignment.out.g_9_umi00_g_11
g_9_outputDir11_g_11 = readsUMIsAssignment.out.g_9_outputDir11_g_11
(g_9_outputDir11_g_12,g_9_outputDir11_g_14) = [g_9_outputDir11_g_11,g_9_outputDir11_g_11]


draftConsensusCalling(g_9_umi00_g_11.groupTuple(by:0).map { it -> it[1]}.flatten().distinct().splitCsv( sep: ' ').distinct().flatten().distinct().splitCsv( sep: '|'),g_9_outputDir11_g_11.flatten().toList())
g_11_umi00_g_12 = draftConsensusCalling.out.g_11_umi00_g_12
g_11_outputDir11 = draftConsensusCalling.out.g_11_outputDir11
g_11_umi20_g_14 = draftConsensusCalling.out.g_11_umi20_g_14


if (!(params.QC == "yes")){
g_9_outputDir11_g_12.set{g_12_outputDir00}
} else {

QC(g_11_umi00_g_12.groupTuple(by:0).map{ it -> it[0]},g_9_outputDir11_g_12.flatten().toList())
g_12_outputDir00 = QC.out.g_12_outputDir00
}


consensusPolishing(g_11_umi20_g_14,g_9_outputDir11_g_14.flatten().toList())
g_14_umi00_g_15 = consensusPolishing.out.g_14_umi00_g_15
g_14_outputDir11 = consensusPolishing.out.g_14_outputDir11


primersTrimming(g_14_umi00_g_15.groupTuple(by:0))
g_15_outputDir00 = primersTrimming.out.g_15_outputDir00


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
