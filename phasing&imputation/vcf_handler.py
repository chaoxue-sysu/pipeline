## VCF operation

import sys
import re
def remove_duplicate_snp(vcf,out):
    snps=set()
    bw=open(out,'w')
    count=0
    with open(vcf,'r') as br:
        for line in br:
            if line.startswith('#'):
                bw.write(line)
                continue
            arr=line.strip().split('\t')
            snp='_'.join([arr[x] for x in [0,1,3,4]])
            if snps.__contains__(snp):
                print(f'remove duplicate: {snp}')
                count+=1
                continue
            snps.add(snp)
            bw.write(line)
    bw.close()
    print(f'summary: remove {count} SNPs')

def add_sex_info(ped_prefix,sex,out_prefix):
    id_sex={}
    with open(sex,'r') as br:
        for line in br:
            arr=line.strip().split('\t')
            id_sex[arr[0]]=arr[1]
    bw=open(out_prefix+'.ped','w')
    ic=0
    ac=0
    with open(ped_prefix+'.ped','r') as br:
        for line in br:
            ac+=1
            arr=line.strip().split('\t')
            iid=arr[1]
            if not id_sex.__contains__(iid):
                continue
            bw.write('\t'.join(arr[0:4]+[id_sex[iid]]+arr[5:])+'\n')
            ic+=1
    bw.close()

    ## copy .map file
    bw=open(out_prefix+'.map','w')
    with open(ped_prefix+'.map','r') as br:
        for line in br:
            bw.write(line)
    bw.close()

    print(f'remain {ic} samples, original {ac}')


def remove_ped_duplicate_snp(pef_prefix,out_prefix):
    snps=set()
    bw=open(f'{out_prefix}.map','w')
    count=0
    idx=-1
    idxs=[]
    with open(f'{pef_prefix}.map','r') as br:
        for line in br:
            idx+=1
            arr=line.strip().split('\t')
            snp='_'.join([arr[x] for x in [0,3]])
            if snps.__contains__(snp):
                print(f'remove duplicate: {snp}')
                count+=1
                continue
            snps.add(snp)
            idxs.append(idx)
            bw.write(line)
    bw.close()
    print(f'summary: remove {count} SNPs')
    bw=open(f'{out_prefix}.ped','w')
    with open(f'{pef_prefix}.ped','r') as br:
        for line in br:
            idx+=1
            arr=re.split('\s+',line.strip())
            sep1=' '
            sep2=' '
            bw.write(sep1.join(arr[0:6]+[sep2.join(arr[x*2+6:x*2+8]) for x in idxs])+'\n')
    bw.close()


def get_chr_length(ref_fa,out):
    bw=open(out,'w')
    for line in open(ref_fa,'r'):
        if not line.startswith('>'):
            continue
        arr=re.split('\s+',line.strip())
        try:
            info=arr[0].replace('>','')+'\t'+arr[2].split(':')[-2]
        except:
            info=arr[0].replace('>','')+'\t'+'NA'
        print(info)
        bw.write(info+'\n')
    bw.close()




if __name__=='__main__':
    # remove_duplicate_snp(sys.argv[1],sys.argv[2])
    # add_sex_info(sys.argv[1],sys.argv[2],sys.argv[3])
    # 'plink --geno 0.1 --hwe 1e-5 --maf 0.01 --mind 0.1'
    # remove_ped_duplicate_snp(sys.argv[1],sys.argv[2])
    get_chr_length(sys.argv[1],sys.argv[2])
