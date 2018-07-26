#!/bin/bash
# pipe usage:
# user@domain: path_to_pipe/pipe.sh <file1> <file2>
date
pipe_version='v1.3'

# executives 
GENOME_PATH="/home/zhangbo/Genome/mm10/Bismark_bowtie2_index_mm10"
LAMBDA_PATH='/home/zhangbo/Genome/Lambda'
G_size='/home/Resource/Genome/mm10/mm10.chrom.sizes'
BISMARK_PATH='/home/zhangbo/Tools/bismark_v0.14.5'
PRESEQ_PATH='/usr/bin'
Rscript_PATH='/home/zhangbo/Tools/IHEC_WGBS'
PhiX_PATH='/home/Resource/Genome/PhiX174/bwa_index/phiX174.fa'
find_id='/home/shaopengliu/pipe_script/wgbs/find_image_ID_digest.sh'
cutadapt='/home/shaopengliu/tools/anaconda3/bin/cutadapt'

# read parameters
while getopts o:p:a:b:h  opts
do 
    case "$opts" in
    o) file1="$OPTARG";;    # PE read 1, or the SE file, or the sra file
    p) file2="$OPTARG";;    # PE read 2. 
    a) adapter_1="$OPTARG";;   # add adapter1
    b) adapter_2="$OPTARG";;    # add adapter2
    h) echo "
usage:  path-to-pipe/pipe.sh -o <read_file1>  -p <ead_file2> -a <adapter_1> -b <adapter_2> 

Options:    -a Sequence of an adapter ligated to the 3' end of the first read. Default TruSeq adapter "AGATCGGAAGAGC"

            -b 3' adapter to be removed from second read in a pair. Default TruSeq adapter "AGATCGGAAGAGC"        
"
exit;;
    [?]) echo "Please use -h to see help information";;
    esac
done

host="zhanglab/wgbs base"
md5=`md5sum $0 | awk '{print $1}'`

if [ -z "$adapter_1" ]
    then
    adapter_1="AGATCGGAAGAGC"
    adapter_2="AGATCGGAAGAGC"
fi

# preparation
if [[ $file1 == *.fastq* ]]; then
    name=`echo ${file1%.fastq*}`
    name=`echo ${name%_1}`
    elif [[ $file1 == *fq.gz ]]; then
    name=`echo ${file1%.fq.gz}`
    name=`echo ${name%_1}`
    else
    echo "please check if the file name has suffix .fastq.gz, .fastq, or .fq.gz"
    exit
fi

s0_wgbs_pre(){
    mkdir 'Processed_'$name
    ln -rs $file1 ./'Processed_'$name
    ln -rs $file2 ./'Processed_'$name
    cd ./'Processed_'$name/
    mkdir 'QC_WGBS_data_collection_'$name
    touch QC_pipe_processing.log

    # start record
    date >> QC_pipe_processing.log
    echo "Target file is $file1 $file2" >> QC_pipe_processing.log
    echo "mm10 WGBS PE data" >> QC_pipe_processing.log
    echo " " >> QC_pipe_processing.log
}


s1.1_cutadapt() {
    # TruSeq adapter
    $cutadapt -j 6 -a $adapter_1 -A $adapter_2 -q 15,10 --minimum-length 36  -o 'Trimed_'$name'_1.fq' -p 'Trimed_'$name'_2.fq'  $file1 $file2  > step1.1_$name'_cutadapt_PE.trimlog' \
        && echo "step1.1, cutadapt with adapter $adapter_1 $adapter_2 successful" >> QC_pipe_processing.log \
        || echo "step1.1, cutadapt with adapter $adapter_1 $adapter_2 fail......" >> QC_pipe_processing.log

    # Nextera kit
    #cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -q 15,10 --minimum-length 36  -o 'Trimed_'$name'_1.fq' -p 'Trimed_'$name'_2.fq'  *fq.gz  > step1.1_$name'_cutadapt_PE.trimlog'

    find . -type l | xargs rm
    echo 'finish removing adapter, start mapping with bismark...'
}


s2_bismark(){
    head -40000000 'Trimed_'$name'_1.fq' > test1.fq
    head -40000000 'Trimed_'$name'_2.fq' > test2.fq

    # 2.1
    echo 'Test phiX contamination...'
    bwa mem -t 6 $PhiX_PATH 'test'*'.fq' | samtools view -bS -F 4 - -o step2.1_phiX_$name'.mapped.bam' \
        && echo "step2.1, bwa phiX mapping successful" >> QC_pipe_processing.log \
        || echo "step2.1, bwa phiX mapping fail......" >> QC_pipe_processing.log
    phi=`samtools view step2.1_phiX_$name'.mapped.bam' |less |wc -l`
    phi_rate=`echo "scale=2; $phi/10000000" |bc -l `
    echo -e "Reads map to phiX genome:\t$phi" > step2.1_phiX_$name'.mapping.txt'
    echo -e "Reads map to phiX genome ratio:\t$phi_rate" >> step2.1_phiX_$name'.mapping.txt' \
        && mv step2.1_phiX_$name'.mapping.txt' 'QC_WGBS_data_collection_'$name

    # 2.2
    echo 'start mapping reads to Lambda genome with bismark...'
    $BISMARK_PATH/bismark --bowtie2  -p 4 --bam -B 'step2.2_Lambda_genome_alignment'  $LAMBDA_PATH  -1 test1.fq  -2 test2.fq  \
        && rm test*.fq \
        && mv step2.2_Lambda_genome_alignment_PE_report.txt 'QC_WGBS_data_collection_'$name \
        && echo "step2.2, bowtie2 Lambda mapping successful" >> QC_pipe_processing.log \
        || echo "step2.2, bowtie2 Lambda mapping fail......" >> QC_pipe_processing.log

    # 2.3
    echo 'start mapping reads to mm10 genome with bismark...'
    $BISMARK_PATH/bismark --bowtie2 --multicore 2 -p 2 --bam $GENOME_PATH -1 'Trimed_'$name'_1.fq'  -2 'Trimed_'$name'_2.fq' \
        && mv "Trimed_"$name"_1.fq_bismark_bt2_PE_report.txt" "step2.3_trimed_"$name"_bismark_bt2_PE_report.txt" \
        && mv 'Trimed_'$name'_1.fq_bismark_bt2_pe.bam' 'step2.3_trimed_'$name'_bismark_bt2_pe.bam' \
        && mv "step2.3_trimed_"$name"_bismark_bt2_PE_report.txt" 'QC_WGBS_data_collection_'$name \
        && echo "step2.3, bismark bowtie2 mapping successful" >> QC_pipe_processing.log \
        || echo "step2.3, bismark bowtie2 mapping fail......" >> QC_pipe_processing.log

    # 2.4
    echo 'start to estimate PCR duplicates in alinged bam file...'
    $BISMARK_PATH/deduplicate_bismark -p --bam 'step2.3_trimed_'$name'_bismark_bt2_pe.bam' \
        && mv step2.3_trimed_$name'_bismark_bt2_pe.deduplicated.bam'  step2.4_trimed_$name'_bismark_bt2_pe.deduplicated.bam' \
        && mv step2.3_trimed_$name'_bismark_bt2_pe.deduplication_report.txt'  step2.4_trimed_$name'_bismark_bt2_pe.deduplication_report.txt' \
        && samtools view -h -b -q 10 -U step2.4_trimed_$name'_Deduplcated_Qless10.bam' 'step2.4_trimed_'$name'_bismark_bt2_pe.deduplicated.bam' -o 'step2.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bam' \
        && mv step2.4_trimed_$name'_bismark_bt2_pe.deduplication_report.txt' 'QC_WGBS_data_collection_'$name \
        && rm step2.4_trimed_$name'_bismark_bt2_pe.deduplicated.bam' \
        && echo "step2.4, bismark estimate PCR dup successful" >> QC_pipe_processing.log \
        || echo "step2.4, bismark estimate PCR dup fail......" >> QC_pipe_processing.log

    # 2.5
    echo 'start to estimate library complexity...'
    awk '{ if (NR%4==2) print substr($0,1,20); }' 'Trimed_'$name'_1.fq' | sort | uniq -c| awk '{ print $1 }' > counts_1.txt
    awk '{ if (NR%4==2) print substr($0,1,20); }' 'Trimed_'$name'_2.fq' | sort | uniq -c| awk '{ print $1 }' > counts_2.txt

    $PRESEQ_PATH/preseq lc_extrap -x 100 -v -Q -o Complexity_future_yield_1.txt -V counts_1.txt \
        && $PRESEQ_PATH/preseq lc_extrap -x 100 -v -Q -o Complexity_future_yield_2.txt -V counts_2.txt \
        && rm counts_1.txt counts_2.txt \
        && rm Trimed_$name*'.fq' \
        && echo "step2.5, preseq estimate successful" >> QC_pipe_processing.log \
        || echo "step2.5, preseq estimate fail......" >> QC_pipe_processing.log

    # 2.6 insertion size
    echo 'start to estimate insertion size...'
    samtools stats 'step2.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bam' | grep ^IS | cut -f 2-3 > 'QC_WGBS_data_collection_'$name/step2.6_insertion_size_$name'.txt' \
        && echo "step2.6, insertion size estimate successful" >> QC_pipe_processing.log \
        || echo "step2.6, insertion size estimate fail......"
}


s3_qc(){
    # 3.1
    echo 'start to estimate GC-bias ...'
    samtools view -h 'step2.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bam' |less |head -100000000> temp.sam \
        && samtools view -Sb temp.sam >temp.bam \
        && bedtools bamtobed -i temp.bam -bedpe > temp.10M.deduplicated.bedpe \
        && awk '{$7=$6-$2} 1' temp.10M.deduplicated.bedpe |cut -d ' ' -f 1,2,6,7 >temp.insert \
        && bedtools bamtobed -i temp.bam  > temp.10M.deduplicated.bed \
        && awk '$1=="chr1"' temp.10M.deduplicated.bed > temp.chr1.bed \
        && bedtools coverage -a $Rscript_PATH/mm10.GCpercentage.chr1.100bp.bin.bg -b temp.chr1.bed -counts > step3.1_cov.chr1.txt \
        && Rscript $Rscript_PATH/GCbias_all.R step3.1_cov.chr1.txt temp.insert \
        && rm temp.* \
        && mv Insertion_size_distribution.pdf step3.1_insertion_size_distribution_$name'.pdf' \
        && rename "s/GC_Percentage_vs_coverage_rho/step3.1_GC_Percentage_vs_coverage_rho_$name/" GC_Percentage_vs_coverage_rho*pdf \
        && rm step3.1_cov.chr1.txt \
        && mv step3.1*pdf 'QC_WGBS_data_collection_'$name \
        && echo "step3.1, GC-bias estimate successful" >> QC_pipe_processing.log \
        || echo "step3.1, GC-bias estimate fail......" >> QC_pipe_processing.log

    # 3.2
    echo 'start to extract methylation level...'
    $BISMARK_PATH/bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --multicore 6 --no_overlap  --gzip --comprehensive --merge_non_CpG --report --bedGraph --buffer_size 10G --genome_folder $GENOME_PATH 'step2.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bam' \
        && mv CpG_context_step2.4_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.txt.gz' step3.2_CpG_context_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.txt.gz' \
        && mv Non_CpG_context_step2.4_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.txt.gz' step3.2_Non_CpG_context_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.txt.gz' \
        && mv step2.4_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.bam_splitting_report.txt' 'QC_WGBS_data_collection_'$name/step3.2_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.bam_splitting_report.txt' \
        && mv step2.4_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.M-bias.txt' 'QC_WGBS_data_collection_'$name/step3.2_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.M-bias.txt' \
        && rename 's/step2.4_/step3.2_/' step2.4*gz \
        && echo "step3.2, GC content estimate successful" >> QC_pipe_processing.log \
        || echo "step3.2, GC content estimate fail......" >> QC_pipe_processing.log

    # 3.2.2 plot CpG results
    if [ -f 'QC_WGBS_data_collection_'$name/step3.2_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.M-bias.txt' ]; then
        cd 'QC_WGBS_data_collection_'$name 
        input_file=step3.2_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.M-bias.txt' 
        total_length=`wc -l $input_file | awk '{print $1}'` 
        input_length=$((total_length/6 - 4 - 5))  # 4 for symbols and space, 5 for trim 
        grep_length=$((input_length+2))

        grep "CpG context (R1)" -A $grep_length $input_file | sed '1,3d' > cpg_r1.txt
        grep "CpG context (R2)" -A $grep_length $input_file | sed '1,3d' > cpg_r2.txt

        grep "CHG context (R1)" -A $grep_length $input_file | sed '1,3d' > chg_r1.txt
        grep "CHG context (R2)" -A $grep_length $input_file | sed '1,3d' > chg_r2.txt

        grep "CHH context (R1)" -A $grep_length $input_file | sed '1,3d' > chh_r1.txt
        grep "CHH context (R2)" -A $grep_length $input_file | sed '1,3d' > chh_r2.txt

        cat <<EOF > temp_cpg_plot.R
        # plot wgbs bias 
        args = commandArgs()
        input_length=as.integer(args[6])

        # read table
        read.table("cpg_r1.txt", header=F, col.names = c("pos","methyl","unmethyl","percent","coverage")) -> cpg1
        read.table("cpg_r2.txt", header=F, col.names = c("pos","methyl","unmethyl","percent","coverage")) -> cpg2
        cpg2\$pos=cpg2\$pos+input_length

        read.table("chg_r1.txt", header=F, col.names = c("pos","methyl","unmethyl","percent","coverage")) -> chg1
        read.table("chg_r2.txt", header=F, col.names = c("pos","methyl","unmethyl","percent","coverage")) -> chg2
        chg2\$pos=chg2\$pos+input_length

        read.table("chh_r1.txt", header=F, col.names = c("pos","methyl","unmethyl","percent","coverage")) -> chh1
        read.table("chh_r2.txt", header=F, col.names = c("pos","methyl","unmethyl","percent","coverage")) -> chh2
        chh2\$pos=chh2\$pos+input_length

        # plot together
        jpeg("step3.2_methylation_percentage_plot_$name.jpg", width=10, height=6, units = 'in', res = 300)
        plot(cpg1\$pos, cpg1\$percent, col="darkseagreen3", lty=1, pch=1, type='b', lwd=2, ylim=c(-10, 100), xlim=c(-10, 2*input_length), xlab=paste("position",input_length, sep="_"), ylab="percentage", main="methylation_percentage_plot")
        lines(cpg2\$pos, cpg2\$percent,col="darkseagreen3", lty=1, pch=1, type='b', lwd=2)

        lines(chg1\$pos, chg1\$percent, col="goldenrod1", lty=3,pch=2, type='b', lwd=2)
        lines(chg2\$pos, chg2\$percent, col="goldenrod1", lty=3,pch=2, type='b', lwd=2)

        lines(chh1\$pos, chh1\$percent, col="cornflowerblue", lty=4,pch=3, type='b', lwd=2)
        lines(chh2\$pos, chh2\$percent, col="cornflowerblue", lty=4,pch=3, type='b', lwd=2)

        legend(1, 70, legend=c("CpG", "CHG","CHH"), col=c("darkseagreen3","goldenrod1","cornflowerblue"), lty=c(1,3,4),lwd=c(2,2,2),cex=0.8)
        dev.off()
EOF

        Rscript temp_cpg_plot.R $input_length && rm temp_cpg_plot.R \
            && rm cpg_r1.txt cpg_r2.txt chg_r1.txt chg_r2.txt chh_r1.txt chh_r2.txt \
            && echo "step3.2.2, plot CpG sucessful" >> ../QC_pipe_processing.log

        cd ..
    else
        echo "step3.2.2, plot CpG fail......" >> QC_pipe_processing.log
    fi

    # 3.3 complexity
    total_read=`grep "Pairs written" *trimlog | awk '{print $5}' | sed 's/,//g'`
    Rscript $Rscript_PATH/WGBS.complex_plot.PE.r  $name $total_read \
        && mv $name'_Complexity_lc_extrap.pdf' step3.3_trimed_$name'_complexity_lc_extrap.pdf' \
        && mv step3.3_trimed_$name'_complexity_lc_extrap.pdf' 'QC_WGBS_data_collection_'$name \
        && echo "step3.3, complexity plot from preseq successful" >> QC_pipe_processing.log \
        || echo "step3.3, complexity plot from preseq fail......" >> QC_pipe_processing.log

    # 3.4 methylation level
    less step3.2_trimed_$name'_bismark_bt2_pe.deduplicated.Q10.bismark.cov.gz' | awk -v OFS='\t' '$7=$5+$6' |awk -v OFS='\t' '$2=$3-1'> 'Temp.cov' \
        && cut -f 1-4 'Temp.cov' > 'step3.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bismark.methyl_calling.bg' \
        && cut -f 1-3,7 'Temp.cov' > 'step3.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bismark.Read_density.bg' \
        && bedSort 'step3.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bismark.methyl_calling.bg' 'step3.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bismark.methyl_calling.bg' \
        && bedSort 'step3.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bismark.Read_density.bg' 'step3.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bismark.Read_density.bg' \
        && detect_c=`wc -l 'step3.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bismark.Read_density.bg' | cut -d ' ' -f 1` \
        && C_cover=`awk '{sum += $4 } END {print sum}' 'step3.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bismark.Read_density.bg'` \
        && C_over20=`awk '$4>20' 'step3.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bismark.Read_density.bg' |wc -l ` \
        && echo -e "detect_c\tC_cover\tC_over20" > 'QC_WGBS_data_collection_'$name/step3.4_temp_store.txt \
        && echo -e "$detect_c\t$C_cover\t$C_over20" >> 'QC_WGBS_data_collection_'$name/step3.4_temp_store.txt \
        && Rscript $Rscript_PATH/CG.coverage.r Temp.cov $name \
        && rm Temp.cov \
        && mv $name'_C_in_CpG_coverage.pdf' step3.4_trimed_$name'_C_in_CpG_coverage.pdf' \
        && mv $name'_methylation_level_Coverage10.pdf' step3.4_trimed_$name'_methylation_level_Coverage10.pdf' \
        && mv step3.4*pdf 'QC_WGBS_data_collection_'$name \
        && echo "step3.4, methylation level estimate successful" >> QC_pipe_processing.log \
        || echo "step3.4, methylation level estimate fail......" >> QC_pipe_processing.log

    # 3.4, gat C_over bin
    for i in `seq 1 1 10`; do
        echo "C_over$((i*2))" >> cover_item.txt
        thres=$((i*2))
        echo `awk -v thres=$thres '$4>thres' 'step3.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bismark.Read_density.bg' |wc -l` >> cover_value.txt
    done

    cat <(cat cover_item.txt | tr "\n" "\t" | sed -e '$a\')  <(cat cover_value.txt | tr "\n" "\t" | sed -e '$a\' ) > 'QC_WGBS_data_collection_'$name/step3.4_CpG_coverage_$name'.txt' \
        && rm cover_item.txt cover_value.txt


    for ea in *bg ; do bgzip $ea; done
    for ea in *bg.gz ; do tabix -p bed $ea; done
}


s4_summarize() {
    # 4.1
    echo 'clean files'
    [ -f Complexity_future_yield_1.txt ] && rename 's/Complexity_future_yield/step2.5_preseq_yield/' Complexity_future_yield_*.txt

    # 4.2
    echo 'summarize results'
    cd 'QC_WGBS_data_collection_'$name

    Lambda_rate=`grep 'Mapping efficiency' step2.2_Lambda_genome_alignment_PE_report.txt |cut -f 2`
    L_total=`grep "Total number of C's analysed:" step2.2_Lambda_genome_alignment_PE_report.txt |cut -f 2`
    mc_cg=`grep "Total unmethylated C's in CpG context" step2.2_Lambda_genome_alignment_PE_report.txt |cut -f 2`
    mc_chg=`grep "Total unmethylated C's in CHG context" step2.2_Lambda_genome_alignment_PE_report.txt |cut -f 2`
    mc_chh=`grep "Total unmethylated C's in CHH context" step2.2_Lambda_genome_alignment_PE_report.txt |cut -f 2`

    BS_c_rate=`echo "scale=4; ($mc_cg+$mc_chg+$mc_chh)/$L_total" |bc -l`
    head -13 "step2.3_trimed_"$name"_bismark_bt2_PE_report.txt" > sm2

    echo "Percentage of Lambda reads in library:    $Lambda_rate" > sm3
    echo "Bisulfite convertion rate:    $BS_c_rate" >> sm3

    q10_read=`samtools view ../'step2.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bam' | wc -l `
    r_length=`samtools view ../'step2.4_trimed_'$name'_bismark_bt2_pe.deduplicated.Q10.bam' | head -1 |cut -f 6 |cut -d 'M' -f 1`
    genome_coverage=`echo "scale=2; $q10_read*2*$r_length/2600000000" | bc -l`
    c_coverage=`echo "scale=2; $C_cover/$detect_c" | bc -l`
    echo "Non-redundant uniquely mapped read (useful reads with Q >10):   $q10_read" >sm4
    echo "Genome coverage:   $genome_coverage" >> sm4
    echo "Detected C in CpGs:   $detect_c" >> sm4
    echo "C in CpG averaged coverage:    $c_coverage" >> sm4
    echo "C in CpG coverage >20X:   $C_over20" >>sm4 
    cat step2.1_phiX_$name'.mapping.txt' sm2 sm3 step2.4_trimed_$name'_bismark_bt2_pe.deduplication_report.txt' sm4> QC_WGBS_$name'.summary.txt'
    rm sm*

    # 4.3
    echo 'get one line summarize table'
    ## cutadapt data
    temp=`grep "Total read pairs processed:" ../step1.1_*trimlog | awk '{print $5}'` \
        && raw_reads=`echo ${temp//,}` \
        && temp2=`grep "Pairs written" ../step1.1_*trimlog | awk '{print $5}'` \
        && written_reads=`echo ${temp2//,}`

    ## phiX
    phiX_read=`head -1 step2.1_phiX_$name'.mapping.txt' | cut -f 2`
    phiX_ratio=`tail -1 step2.1_phiX_$name'.mapping.txt' | cut -f 2`

    ## mapping status
    unique_mapped=`head -13 "step2.3_trimed_"$name"_bismark_bt2_PE_report.txt" | grep "Number of paired-end alignments with a unique best hit" | cut -f 2`
    unique_percent=`head -13 "step2.3_trimed_"$name"_bismark_bt2_PE_report.txt" | grep "Mapping efficiency" | cut -f 2`
    multiple_mapped=`head -13 "step2.3_trimed_"$name"_bismark_bt2_PE_report.txt" | grep "Sequence pairs did not map uniquely" | cut -f 2`
    non_mapped=`head -13 "step2.3_trimed_"$name"_bismark_bt2_PE_report.txt" | grep "Sequence pairs with no alignments under any condition" | cut -f 2`

    ## after-align dup
    after_dup=`cat step2.4_trimed_$name'_bismark_bt2_pe.deduplication_report.txt' | grep "Total number duplicated alignments removed:" | cut -f 2 | awk '{print $2}' | sed 's/[)(]//g'`

    ## methylation level
    detect_c=`sed -n '2p' step3.4_temp_store.txt | cut -f 1`
    C_cover=`sed -n '2p' step3.4_temp_store.txt | cut -f 2`
    C_over20=`sed -n '2p' step3.4_temp_store.txt | cut -f 3`

    ## adding:
    echo -e "name\traw_reads\twritten_reads\tmultiple_mapped\tunique\tuniq_percent\tnon_mapped\tnon_redun_uniquely_mapped_single_ends\tafter_align_dup_percent\tgenome_coverage\tratio_map2phiX\tLambda_reads_percent\tBS_convertion_rate\tdetected_C_in_CpG\tC_in_CpG_avg_coverage\tC_in_CpG_cov_gt20X" > QC_summary_$name'.txt'
    echo -e "$name\t$raw_reads\t$written_reads\t$multiple_mapped\t$unique_mapped\t$unique_percent\t$non_mapped\t$q10_read\t$after_dup\t$genome_coverage\t$phiX_ratio\t$Lambda_rate\t$BS_c_rate\t$detect_c\t$c_coverage\t$C_over20" >> QC_summary_$name'.txt'

    if [ -z "$name" ] || [ -z "$raw_reads" ] || [ -z "$written_reads" ] || [ -z "$multiple_mapped" ] || [ -z "$unique_mapped" ] || [ -z "$unique_percent" ] || [ -z "$non_mapped" ] || [ -z "$q10_read" ] || [ -z "$genome_coverage" ] || [ -z "$phiX_ratio" ] || [ -z "$Lambda_rate" ] || [ -z "$BS_c_rate" ] || [ -z "$detect_c" ] || [ -z "$c_coverage" ] || [ -z "$C_over20" ]; then
        echo "step4.2, checking variables fail......" >> ../QC_pipe_processing.log 
    else
        echo "step4.2, checking variables successful" >> ../QC_pipe_processing.log
    fi

    # 4.4 generate json
    cat <<EOF > wgbs_summarize.R
#!/usr/bin/env Rscript
args=commandArgs()
library(jsonlite)

name=args[6]
pipev=args[7]
run_time=args[8]
image_id=args[9]
md5=args[10]
genome="mm10"
read_type="PE"

### delete
print(paste("name","is",name))
print(paste("version","is",pipev))
print(paste("time","is",run_time))
print(paste("id","is",image_id))
print(paste("md5","is", md5))
###

# read QC table
read.table(paste("QC_summary_", name, ".txt", sep=""), header=T, stringsAsFactors = FALSE) -> qc_table

# part1, data information
part1=data.frame(name, genome, read_type, pipev, image_id, md5, run_time)
colnames(part1)=c("file_name","genome","read_type","pipe_version","Docker_image_id","bash_script_MD5","running_time")
file=list('data_information'=part1)

# part2, pre-alignment
part2=data.frame("cutadapt","1.16",qc_table\$raw_reads,qc_table\$written_reads, round(qc_table\$written_reads/qc_table\$raw_reads,digit=2))
colnames(part2)=c("program1", "program1_version","raw_reads","written_reads_by_cutadapt", "written_reads_ratio")
file=append(file, list('pre_alignment_stats'=part2))

# part3, mapping status
part3=data.frame("bismark_bw2", "v0.14.4", "-q --score-min L,0,-0.2 -p 2 --reorder --ignore-quals --no-mixed --no-discordant --maxins 500" , "samtools", "1.6", qc_table\$raw_reads, qc_table\$non_mapped, qc_table\$unique, as.numeric(sub("%","", qc_table\$uniq_percent))/100 , qc_table\$multiple_mapped, qc_table\$non_redun_uniquely_mapped_single_ends)
colnames(part3)=c("alignment_program", "alignment_program_version", "alignment_program_parameters", "post_alignment_program", "post_alignment_program_version","total_reads", "reads_aligned_0_time", "reads_aligned_exactly_1_time","unique_mapped_ratio","reads_aligned_greater_than_1_time", "non_redun_uniquely_mapped_single_ends")
file=append(file, list("mapping_stats"=part3))

# part4, library complexity
part4=data.frame(qc_table\$ratio_map2phiX, as.numeric(sub("%","", qc_table\$Lambda_reads_percent))/100, as.numeric(sub("%","", qc_table\$after_align_dup_percent))/100)
colnames(part4)=c("reads_ratio_mapped_to_phiX", "Lambda_reads_ratio", "after_alignment_duplicates_ratio_by_bismark")
file=append(file, list("library_complexity"=part4))

# part5, bisulfite
part5=data.frame(qc_table\$BS_convertion_rate, qc_table\$genome_coverage, qc_table\$detected_C_in_CpG, qc_table\$C_in_CpG_avg_coverage, qc_table\$C_in_CpG_cov_gt20X )
colnames(part5)=c("bisulfite_convertion_date", "genome_coverage", "detected_C_in_CpGs", "C_in_CpG_averaged_coverage", "C_in_CpG_coverage>20X")
file=append(file, list("bisulfite"=part5))

# part6, CpG coverage saturation
read.table(paste("step3.4_CpG_coverage_", name, ".txt", sep=""), header=T, stringsAsFactors = FALSE) -> cpg_table
part6=data.frame(cpg_table)
colnames(part6)<- paste(gsub("C_over","C_in_CpG_coverage>", colnames(part6)), "X", sep="")
file=append(file, list("CpG_coverage"=part6))


capture.output(toJSON(file,pretty=T),file=paste(name,"report.json",sep='_'))
EOF

    time=`head -1 ../QC_pipe_processing.log | sed 's/ /_/g'`
    image_id=`bash $find_id $host  2> /dev/null | awk '{print $2}'`
    if [ -z "$image_id" ]
        then
        image_id="failed_to_get_id"
    fi
    Rscript wgbs_summarize.R $name $pipe_version $time $image_id $md5 \
        && echo "step4.2, summarizing results successful" >> ../QC_pipe_processing.log \
        || echo "step4.2, summarizing results fail......" >> ../QC_pipe_processing.log
    rm wgbs_summarize.R

    # edit json(delete []; add "," after the ending }; replace "! !" with {} which I don't find on WGBS json; replace "? ?" with {} which I don't find on WGBS json; )
    sed 's/\[/{/g' $name'_report.json' | sed '/    {/d' | sed '/\]/d' | sed 's/\\r//g' |  sed 's/    }/    },/g' | tac | sed '2s/},/}/g' | tac > QC_$name'.json' && rm $name'_report.json' && mv QC_$name'.json' ../ && mv *summary*txt ../
    #didn't use:  
    #sed 's/"!/{/g' | sed 's/!"/}/g' | sed 's/"?/[/g' | sed 's/?"/]/g' | sed 's/@/"/g' | tac | sed '3s/},/}/g' | sed '1,2d' | tac | cat - <(echo "  },") <(sed '1d' $pipe_path'/../atac_ref/mm10_encode_pe/encode_pe.json') | sed 's/\\r//g' " > QC_$name'.json'

    echo "processing finished"
    echo "finished" >> ../QC_pipe_processing.log
    date >> ../QC_pipe_processing.log
    cd ..
}




#run pipe
########################################
s0_wgbs_pre
s1.1_cutadapt
s2_bismark
s3_qc
s4_summarize









