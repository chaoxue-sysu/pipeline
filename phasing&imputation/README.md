## GWAS基因型数据填充流程
基因型填充是指利用参考的单体型数据推测样本中未被测到的位点基因型。
对于非全基因组测序来讲，基因型数据填充是进行GWAS等下游分析有着重要的意义。

本文基于一组精神分裂症数据为例进行填充流程说明。

#### 1. 预处理
获得数据为SNP芯片数据，已经预处理为plink 格式的数据（.ped/.map）。但是存在以下几个问题，
需要处理：.ped不包含性别信息（在性染色体的phasing和imuptation中需要用到），
存在重复的位点（SNP位置,等位基因完全相同但SNP ID不同的位点）。
##### 1.1 添加性别信息
用脚本处理：将性别信息添加到.ped文件的第5列
（.ped格式：http://zzz.bwh.harvard.edu/plink/data.shtml#ped）
##### 1.2 移除冗余位点
用脚本处理.ped/.map文件，去除位置信息冗余位点（保留唯一）。
（.map格式：http://zzz.bwh.harvard.edu/plink/data.shtml#map）

#### 2. 质控
基因型数据的质量对下游分析有关键的影响，质控步骤必不可少。质控主要从样本和maker(SNP)
质量两个角度进行。
根据 [Gamazon](https://www.nature.com/articles/ng.3367#methods) 的质控方法:
> Approximately 650,000 SNPs (minor allele frequency (MAF) > 0.05, 
>in Hardy-Weinberg equilibrium (P > 0.05) and with non-ambiguous strand 
>mapping (no A/T or C/G SNPs)) comprised the input set of SNPs for 
>imputation, which was performed on the University of Michigan Imputation 
>Server39,40 with the following parameters: 1000G Phase 1 v3 ShapeIt2 
>(no singletons) reference panel, SHAPEIT phasing and the EUR (European) 
>population. Non-ambiguous-stranded SNPs with MAF > 0.05 and imputation 
>R2 > 0.8 were retained for subsequent analysis. 

这里利用plink进行质控实现：
```bash 
plink --recode ped --file $ped_prefix --out $out_ped_prefix --geno 0.05 --hwe 5e-2 --maf 0.05 --mind 0.05
```
特别提出，当SNP芯片探针设计的无从考证时（测得的SNP可能是对应正链或负链的基因型），
在严格要求正确的情况下，我们只能删除链歧义的SNP位点，即A/T，C/G等位基因的情况。
而非歧义的SNP位点可以通过后续与参考基因组相同位置碱基的比较来判断改位点探针的正负链。
如A/C，参考基因组正(+)链为T，那么可以判定此位点的基因型基于负(-)链。

顺便提一下，Gamazon使用的[University of Michigan Imputation Server](https://imputationserver.sph.umich.edu/)
，这是一个用于基因型填充的界面友好的网页服务器，用户只需要上传基因型文件，并设置相关参数即可。
等到任务在服务器后台执行完成后，用户就可以下载填充结果。其使用的pipeline采用
[Eagle](https://www.nature.com/articles/ng.3571) (Phasing) +
[Minimac4](https://github.com/statgen/Minimac4) (Imputation)，
详见
[imputationserver](https://github.com/genepi/imputationserver/blob/master/docs/pipeline.md)
。参考文献见 [Next-generation genotype imputation service and methods, Nat. Genet.](https://www.nature.com/articles/ng.3656#ref8)
。

| Feature | As summary statistic |A s inclusion criteria |
|--|--|--|
|Missingness per individual|--missing|	--mind N|
|Missingness per marker|	--missing|	--geno N|
|Allele frequency|	--freq	|--maf N|
|Hardy-Weinberg equilibrium|	--hardy|	--hwe Nv|
|Mendel error rates|	--mendel|	--me N M|

#### 3. 分割染色体
按染色体分割基因型文件，可以提升后续分析的效率。
在按染色体分割时，将ped/map里面的染色体编号[0,26]改为容易理解的标识
(与参考单体型染色体一致，方便后续分析):
* __0__，chrUnknown; 
* __1-22__, chr1-22; 
* __23__, chrX_nonPAR; 
* __24__, chrY_nonPAR;
* __25__, 分两段：[60001,2699520],chrX_PAR1; [154931044,155260560],chrX_PAR2;
坐标基于GRCh37;
* __26__, chrMT.

#### 4. 定相（Phasing）
根据 [REF](#ref) ,做基因型填充的最佳流程是：先用比较准确的定相方法对样本基因型进行定相。
本文使用 [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
进行定相。
```python
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
    return command
```

#### 5. 基因型填充（Imputation）
利用 [IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)
进行基因型填充。

```python
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
```



#### <span id = "ref">6. 参考<span>
1. https://zhuanlan.zhihu.com/p/53564008
2. https://www.plob.org/article/13447.html
3. https://www.nature.com/articles/ng.2354


