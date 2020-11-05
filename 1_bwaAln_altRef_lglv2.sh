#BWA ALN

REF_DIR=/DATA/share/pipelines/lgleon/References/Remco_ref

for file in *.f*q.gz
do

   inputfqgz=$file
   #outputbam=${file/%.fastq.gz/.bam}  #  % means that suffixes are matched, # would mean prefixes"
   filenamestem=${file/%.f*q.gz/}
   
   echo "!!!!!!!!!!!!!! Expected final output {$filenamestem.markdup.bam}"

   bwa aln -t 32 -n 2 $REF_DIR/hg19/hg19.fa  $inputfqgz | bwa samse -r '@RG\tID:inputfqgz\tSM:inputfqgz' $REF_DIR/hg19/hg19.fa - $inputfqgz | samtools view -Su - | samtools sort - -o $filenamestem.sorted.bam > $filenamestem.sorted.bam

   picard MarkDuplicates -I $filenamestem.sorted.bam -O $filenamestem.markdup.bam -M $filenamestem.markdup_metrics.txt -ASSUME_SORTED true -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT

done

#mkdir fastq
#mv *fastq.gz fastq

#mkdir bam.sorted
#mv *sorted.bam bam.sorted
#rm *sorted.bam

#mkdir bam.markdup
#mv *markdup.bam bam.markdup

