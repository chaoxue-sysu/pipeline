### pipline for phasing and imputation
# phased by SHAPEIT (https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
# imputed by IMPUTE2 (https://mathgen.stats.ox.ac.uk/impute/impute_v2.htm)

import subprocess
import time
import os
import sys

SHAPEIT='shapeit'
IMPUTE2='impute2'
PLINK='plink'
LOG=os.environ['HOME']+'/tmp'

GENETIC_MAP_HAPLOTYPE='/public1/data/resources/db_tools/phasing_imputation/1000GP_Phase3'


def split_chr_bed(bed_prefix,out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    log_path='%s/%s.%s.log'%(LOG,sys._getframe().f_code.co_name,str(time.time()).replace('.',''))
    log(f'shell log at {log_path}')
    cmds=[]
    for i in range(27):
        cmd=[PLINK]
        cmd.append(f'--make-bed --bfile {bed_prefix} --chr {i} --out {out_dir}/chr{i}')
        command=' '.join(cmd)
        cmds.append(command)
        # print(command)
    batchShellTask(cmds,30,log_path)
    pass

def split_chr_vcf(vcf,out_dir):
    ## ref GRCh37
    ChrX_PAR1=(60001,2699520,'chrX_PAR1')
    ChrX_PAR2=(154931044,155260560,'chrX_PAR2')
    chrids={23:'chrX_nonPAR',24:'chrY_nonPAR',26:'chrMT',0:'chrUnknown'}
    for i in range(1,23):
        chrids[i]=f'chr{i}'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    log_path='%s/%s.%s.log'%(LOG,sys._getframe().f_code.co_name,str(time.time()).replace('.',''))
    log(f'shell log at {log_path}')
    cmds=[]
    for i in range(27):
        if i==25:
            for start_end in (ChrX_PAR1,ChrX_PAR2):
                cmd=[PLINK]
                cmd.append(f'--recode vcf --vcf {vcf} --chr {i} --from-bp {start_end[0]} --to-bp {start_end[1]} --out {out_dir}/{start_end[2]}')
                command=' '.join(cmd)
                cmds.append(command)
        else:
            cmd=[PLINK]
            cmd.append(f'--recode vcf --vcf {vcf} --chr {i} --out {out_dir}/{chrids[i]}')
            command=' '.join(cmd)
            cmds.append(command)
        # print(command)
    batchShellTask(cmds,30,log_path)
    pass


def split_chr_ped(ped_prefix,out_dir):
    ## ref GRCh37
    ChrX_PAR1=(60001,2699520,'chrX_PAR1')
    ChrX_PAR2=(154931044,155260560,'chrX_PAR2')
    chrids={23:'chrX_nonPAR',24:'chrY_nonPAR',26:'chrMT',0:'chrUnknown'}
    for i in range(1,23):
        chrids[i]=f'chr{i}'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    log_path='%s/%s.%s.log'%(LOG,sys._getframe().f_code.co_name,str(time.time()).replace('.',''))
    log(f'shell log at {log_path}')
    cmds=[]
    for i in range(27):
        if i==25:
            for start_end in (ChrX_PAR1,ChrX_PAR2):
                cmd=[PLINK]
                cmd.append(f'--make-bed --file {ped_prefix} --chr {i} --from-bp {start_end[0]} --to-bp {start_end[1]} --out {out_dir}/{start_end[2]}')
                command=' '.join(cmd)
                cmds.append(command)
        else:
            cmd=[PLINK]
            cmd.append(f'--make-bed --file {ped_prefix} --chr {i} --out {out_dir}/{chrids[i]}')
            command=' '.join(cmd)
            cmds.append(command)
        # print(command)
    batchShellTask(cmds,30,log_path)
    pass




def phasing_chrs(type,chr_dir,genetic_map,out_dir,log_dir,thread,batch_nt):
    thread=int(thread)
    batch_nt=int(batch_nt)
    sex_chr=('chrX_nonPAR','chrY_nonPAR')
    for dir in [out_dir,log_dir]:
        if not os.path.exists(dir):
            os.makedirs(dir)
    log_path='%s/%s.%s.log'%(LOG,sys._getframe().f_code.co_name,str(time.time()).replace('.',''))
    cmds=[]
    for chrid in [f'chr{x}' for x in range(1,23)]+['chrX_nonPAR','chrX_PAR1','chrX_PAR2']:
        sex=False
        if sex_chr.__contains__(chrid):
            sex=True
        cmd=phasing_cmd(type,f'{chr_dir}/{chrid}',f'{genetic_map}/genetic_map_{chrid}_combined_b37.txt',f'{out_dir}/{chrid}',thread,f'{log_dir}/{chrid}',sex)
        cmds.append(cmd)
    batchShellTask(cmds,batch_nt,log_path)


def phasing_cmd(type,bed_prefix,genetic_map,out_prefix,thread,log_file,sex_chr=False):
    cmd=[SHAPEIT]
    cmd.append('--force')
    if sex_chr:
        cmd.append('--chrX')
    if type=='vcf':
        cmd.append('-V %s.vcf'%bed_prefix)
    elif type=='bed':
        cmd.append('-B %s'%bed_prefix)
    else:
        raise Exception(f'Not support type: {type}')
    cmd.append('-M %s'%genetic_map)
    cmd.append('-O %s'%out_prefix)
    cmd.append('-T %s'%thread)
    cmd.append('-L %s'%log_file)
    command=' '.join(cmd)
    # os.system(command)
    return command


def impute_chr_haps(haps_dir,ref_haps_dir,genetic_map,out_dir,chr_length,batch_thread):
    ## get the length of all chr
    batch_thread=int(batch_thread)
    chr_len={}
    for line in open(chr_length,'r'):
        arr=line.strip().split('\t')
        try:
            chr_len[f'chr{arr[0]}']=(0,int(arr[1]))
        except:
            continue
    chr_len['chrX_nonPAR']=chr_len['chrX']
    chr_len['chrX_PAR1']=(60001,2699520)
    chr_len['chrX_PAR2']=(154931044,155260560)

    sex_chr=('chrX_nonPAR','chrY_nonPAR')
    mkdirs(out_dir)
    cmds=[]
    block_size=int(5e6)
    for chrid in [f'chr{x}' for x in range(1,23)]+['chrX_nonPAR','chrX_PAR1','chrX_PAR2']:
        sex=False
        if sex_chr.__contains__(chrid):
            sex=True
        chrid_ref=chrid
        if chrid_ref=='chrX_nonPAR':
            chrid_ref='chrX_NONPAR'
        ## split the whole genome into 5M blocks for imputation
        pos=chr_len[chrid]
        blocks=(pos[1]-pos[0])//block_size
        if (pos[1]-pos[0])%block_size!=0:
            blocks+=1
        mkdirs(f'{out_dir}/{chrid}')
        for bidx in range(blocks):
            start_pos=int(pos[0]+block_size*bidx)
            end_pos=int(pos[0]+block_size*(bidx+1))
            out_prefix=f'{out_dir}/{chrid}/IMPUTE2_{chrid}.{start_pos//int(1e6)}-{end_pos//int(1e6)}Mb.gen'
            cmd= impute_haps(out_prefix,
                            f'{haps_dir}/{chrid}.haps',
                            f'{haps_dir}/{chrid}.sample',
                            f'{ref_haps_dir}/1000GP_Phase3_{chrid_ref}.hap.gz',
                            f'{ref_haps_dir}/1000GP_Phase3_{chrid_ref}.legend.gz',
                            f'{genetic_map}/genetic_map_{chrid}_combined_b37.txt',
                            start_pos, end_pos, sex)
            cmds.append(cmd)
    ## submit to run
    log_path='%s/%s.%s.log'%(LOG,sys._getframe().f_code.co_name,str(time.time()).replace('.',''))
    log(f'shell log at {log_path}')
    log(f'all cmds: {len(cmds)}; batch thread: {batch_thread}')
    batchShellTask(cmds,batch_thread,log_path)



def impute_haps(out_prefix,phased_haps,ref_sample,ref_haps,ref_legend,genetic_map,start_pos,end_pos,sex_chr=False):
    cmd=[IMPUTE2]
    if sex_chr:
        cmd.append('-chrX ')
    cmd.append('-phase')
    cmd.append(f'-sample_g {ref_sample}')
    cmd.append('-use_prephased_g')
    cmd.append('-known_haps_g %s'%phased_haps)
    cmd.append('-h %s'%ref_haps)
    cmd.append('-l %s'%ref_legend)
    cmd.append('-m %s'%genetic_map)
    cmd.append('-int %s %s'%(start_pos,end_pos))
    cmd.append('-Ne %s'%20000)
    cmd.append(f'-o {out_prefix}')
    command=' '.join(cmd)
    return command



def mkdirs(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def batchShellTask(all_task,limit_task,Log_file):
    '''
    batch shell task
    :param all_task: all shell task list
    :param limit_task: limit task number simultaneously running
    :type all_task: list
    :type limit_task: int
    :return: null
    '''
    log = open(Log_file, "w")
    task_pool=[]
    task_remain=len(all_task)
    for task in all_task:
        task_remain+=-1
        break_out = True
        p = subprocess.Popen(task, shell=True, stdin=log, stdout=log, stderr=log, close_fds=True)
        task_pool.append(p)
        print(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' '+str(p.pid)+': '+task+' start ...',flush=True)
        if len(task_pool)==limit_task or task_remain==0:
            while break_out:
                for intask_Popen in task_pool:
                    if intask_Popen.poll()!=None:
                        print(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' '+str(intask_Popen.pid)+': '+' finish...',flush=True)
                        task_pool.remove(intask_Popen)
                        break_out = False
                        if task_remain==0:
                            break_out=True
                        if len(task_pool)==0:
                            break_out=False
                        break
    log.close()

def qc_plink(ped_prefix,out_ped_prefix):
    log_path='%s/%s.%s.log'%(LOG,sys._getframe().f_code.co_name,str(time.time()).replace('.',''))
    log(f'shell log at {log_path}')
    cmd=f'plink --recode ped --file {ped_prefix} --out {out_ped_prefix} --geno 0.05 --hwe 1e-5 --maf 0.01 --mind 0.05'
    batchShellTask([cmd],1,log_path)

## static log method
def log(content):
    content=time.strftime("%Y-%m-%d %H:%M:%S [INFO] ", time.localtime(time.time()))+ "%s"%(content)
    # wt = FF.getWriter(self.logPath, True)
    print(content)

def main():
    # import vcf_handler as VH
    # ## remove duplicate site
    # VH.remove_ped_duplicate_snp(sys.argv[1],sys.argv[2])
    # ## QC
    # qc_plink(sys.argv[2],sys.argv[3])
    # ## split by chr
    # split_chr_ped(sys.argv[3],sys.argv[4])
    ## phasing
    # phasing_chrs(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])
    impute_chr_haps(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])




if __name__=='__main__':
    # split_chr_bed(sys.argv[1],sys.argv[2])
    # phasing_chrs_vcf(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
    # split_chr_vcf(sys.argv[1],sys.argv[2])
    # split_chr_ped(sys.argv[1],sys.argv[2])
    main()