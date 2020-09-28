# WGBS_analysis
whole genome bisulfite sequencing analysis:
WGBS pipeline is designed for standardized data processing of mouse WGBS data. It incorporates automatic quality control, generates user friendly files for computational analysis and outputs genome browser tracks for data visualization. To ensure consistent and reproducible data processing, the entire workflow, associated software and libraries are built into a singularity image, which can be run on computational clusters with job submission as well as on stand-alone machines. Pipeline installation requires minimal user input. All the software and genome references used for WGBS data processing are included in the pipeline image. The pipeline supports paired-end data, it accepts FASTQ files, performs alignments, features summary and data visualization.

Current version: `V20191126`   

Advisor: Bo Zhang  
Contributor: Shaopeng Liu  

For any question, please contact Wustl.Zhanglab@gmail.com  

## Usage: 
2-step guideline for pipeline usage:

Please note that this is for **mm10 PE WGBS data only** for now.  

Step1. download singularity container (you only need download the containcer for once, then you can use them directly):  
####  
```bash
# download image from local server:  
wget http://regmedsrv1.wustl.edu/Public_SPACE/bmiao/Public_html/TaRGET_II_pipeline/WGBS/mm10_TaRGET_WGBS_20191126.simg
```

Step2. process data by the singularity image: 
#### Please run at same directory with your data OR the soft link of your data    
```bash
singularity run -B ./:/process <path-to-image> -r <PE> -o <read_file1> -p <read_file2> 
```
**soft link introduction:** For soft link of data, need to add one bind option for singularity, which is ```-B <full-path-of-original-position>:<full-path-of-original-position>```. 
For example:
soft link the data from /scratch to run on folder /home/example. **Please make sure you use the absolute path.**
```
ln -s /scrach/mydata.fastq.gz /home/example;
cd /home/example
singularity run -B ./:/process -B /scratch:/scratch /home/image/mm10_TaRGET_WGBS_20191126.simg -r PE -o read1.fastq.gz -p read2.fastq.gz
```
That's it!

#parameters:  
`-h`: help information  
`-o`: reads file 1, must be ended by .fastq or .fastq.gz  
`-p`: reads file 2, must be ended by .fastq or .fastq.gz  
`-a`: ADAPT1 for cutadapt  
`-b`: ADAPT2 for cutadapt  

# Output files

After running the pipeline, there will be a folder called Processed_${name}, all intermediate files and final output files are stored there. The major outputs produced by the pipeline are:
   1)	JSON file with quality control measurement of user supplied dataset. 
   2)	Log file with processing status of each step .
   3)	Processed files after alignment (.bam) and bismark output files. 
   4)	BedGraph file for visualization purpose

