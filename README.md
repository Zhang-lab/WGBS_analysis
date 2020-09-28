# WGBS_analysis
whole genome bisulfite sequencing analysis

This is for the QC matrix construction, data analysis and visualization for WGBS data.  
Current version: `V1.3`   

Advisor: Bo Zhang  
Contributor: Shaopeng Liu  

For any question, please contact Wustl.Zhanglab@gmail.com  

## Usage: 
Singularity 2-step solution (easiest way)  
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
singularity run -H ./:/scratch wgbs.simg -o <read_file1>  -p <read_file2>  -a <adapter_1> -b <adapter_2>  
```

That's it!

#parameters:  
`-h`: help information  
`-o`: reads file 1, must be ended by .fastq or .fastq.gz  
`-p`: reads file 2, must be ended by .fastq or .fastq.gz  
`-a`: ADAPT1 for cutadapt  
`-b`: ADAPT2 for cutadapt  

