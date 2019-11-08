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

plink进行质控：
```bash 
plink --recode ped --file $ped_prefix --out $out_ped_prefix --geno 0.05 --hwe 1e-5 --maf 0.01 --mind 0.05
```

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


