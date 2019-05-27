# pdb2vall
# Rosetta pdb2vall: 构建自己的fragments_vall数据库

> 参考: 一篇关于Fragment Picker的教程(日文): https://qiita.com/Ag_smith/items/58a7e28c75d2a51367f8
>
> 使用笔记本计算的同学，可以忽略本文，因为构建最新的PBD非冗余库时间过长。

## 前言

Rosetta自带的`vall.jul19.2011.gz`于2011年7月19日发布。截止2018年9月30日期间，蛋白质数据库的数据量已经增长了1倍。在Rosetta3.9版本中，用户可以通过pdb2vall模块自行构建vall数据库，但是操作并不友好，需要大量的修改和配置。因此本文将提供修改后的代码以及提供安装和使用方法。



## 更新方法

**务必注意**|**务必注意**|**务必注意**: 

使用本文的配置方法前，请正确配置好Rosetta Fragment tools的环境，如nr数据库，make_fragments.pl脚本环境等等，详细请参考本人另外一篇博文:

<https://awakenwu.github.io/2019/05/18/Rosetta-Fragment-Picker-短肽片段分选器/>



本文直接采用前人的修改成果并进行CentOS7适配性修改, 首先先下载前人的成果:(如果想自己配置的话)

```shell
cd $ROSETTA/tools/fragment_tools/
git clone https://github.com/BILAB/fragment_tools.git
# 将fragment_tools中的所有文件替换至$ROSETTA/tools/fragment_tools/并覆盖。
```



此处提供经过本人修改pdb2vall:

```shell
cd $ROSETTA/tools/fragment_tools/
git clone https://github.com/awakenwu/pdb2vall.git
# 直接替换掉原文件的pdb2vall即可。
```



以下为必要配置操作步骤:

安装32位环境依赖:

```shell
yum -y install libg2c.so.0
```



#### 1 获取非冗余PDB数据列表

利用PISCES获取非冗余PDB的数据列表。(http://dunbrack.fccc.edu/Guoli/PISCES_ChooseInputPage.php)，其中Maximum percentage identity设置为60()，留下邮箱后，即可获取结果。

打开邮件: 点击 Download sequence ID list file下面的url即可获取满足条件的非冗余PDB ID列表。可用于后面的脚本编写。



#### 2 获取pdb_seqres.txt，ss_dis.txt

下载全网pdb的序列信息和其二级结构信息。

```shell
cd $ROSETTA/tools/fragment_tools/pdb2vall/database/rcsb_data/
mkdir derived_data && cd derived_data
mkdir $ROSETTA/tools/fragment_tools/pdb2vall/database/PDB_uncompressed
wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt
gunzip pdb_seqres.txt.gz
# 运行 get_ss_dis.sh
cd ..
chmod +X get_ss_dis.sh
./get_ss_dis.sh
```



#### 3 安装BioPerl

安装时间还挺久的，半小时左右，取决于网速。

```shell
yum -y install perl-CPAN perl-CPANPLUS
perl -MCPAN -e shell  # 全部选yes即可

# 安装完毕后自动进入CPAN: 并输入以下命令
install CPAN
reload cpan
install Module::Build
o conf commit MB
o conf commit
install Algorithm::Munkres Array::Compare Convert::Binary::C Data::Stag Graph GraphViz Math::Random Bio::AlignIO
force install Bio::Perl
```



#### 4 安装DEPTH

自带的DEPTH貌似有些问题，此处作者推荐重新编译安装:

```shell
cd cd $ROSETTA/tools/fragment_tools
wget http://cospi.iiserpune.ac.in/depth/htdocs/program/depth-2.0.tar.gz
tar zxvf depth-2.0.tar.gz
cd depth-2.0/bin
make
```



#### 5 配置pdb2vall.cfg

只要是配置需要调用的第三方应用的路径。

```shell
cd $ROSETTA/tools/fragment_tools/pdb2vall
vi pdb2vall.cfg
# 根据实际情况进行路径的修改:

[pdb2vall]

rosetta_path = /mnt/sdd/software/rosetta_src_2019.14.60699_bundle
fragment_picker_num_cpus = 40

# 使用old blast
blastpgp_num_cpus = 8
blastpgp = /mnt/sdd/software/rosetta_src_2019.14.60699_bundle/tools/fragment_tools/blast/bin/blastpgp
bl2seq =  /mnt/sdd/software/rosetta_src_2019.14.60699_bundle/tools/fragment_tools/blast/bin/bl2seq
formatdb =  /mnt/sdd/software/rosetta_src_2019.14.60699_bundle、tools/fragment_tools/blast/bin/formatdb

nr =  /mnt/sdd/software/rosetta_src_2019.14.60699_bundle/tools/fragment_tools/databases/nr

depth =  /mnt/sdd/software/rosetta_src_2019.14.60699_bundle/tools/fragment_tools/depth-2.0/bin/DEPTH
depth_num_cpus = 1
```



设置pdb2vall的必要环境:

```shell
vi ~/.zshrc
export PDB_DIR=$ROSETTA/tools/fragment_tools/pdb2vall/database/PDB_uncompressed
#export INET_HOST="localhost"
export BLASTDB=$ROSETTA/tools/fragment_tools/nr
export BLASTMAT=$ROSETTA/tools/blast/data
```



重写/pdb2vall/structure_profile_scripts/make_sequence_fragments.pl

```shell
# 约29行
$ROSETTA_PATH/main/source/bin/fragment_picker.boost_thread.linuxgccrelease
的修改Fragment Picker的实际app名称为$ROSETTA_PATH/main/source/bin/fragment_picker.boost_thread.linuxclangrelease
```



修改pdb2vall.py

```python
#66 67行，将relax和idealize_jd2的app名称修改为实际名称。
idealization_app = ROSETTA_PATH + "main/source/bin/idealize_jd2.mpi.linuxclangrelease"
relax_app = ROSETTA_PATH + "main/source/bin/relax.mpi.linuxclangrelease"
```



#### 6 其他细节修改(如果从我提供的github地址中进行下载的pdb2vall，可以忽视这个步骤。)

1 重写pdb2vall.py

注意脚本应该在Python2环境运行。

```python
# 31行 def download_pdb(pdb_id,dest_dir)函数修改为以下:
def download_pdb(pdb_id,dest_dir):
	print "downloading %s" % ( pdb_id )
	url      = 'https://files.rcsb.org/download/%s.pdb' % (pdb_id.lower())
	dest     = '%s/%s.pdb' % ( os.path.abspath(dest_dir), pdb_id.lower() )
	wget_cmd = 'wget -c -t 0 %s' % ( url )
	os.system(wget_cmd)

	#if remote_host:
	#	subdirectory = pdb_id[1:3]
    #	wget_cmd = 'scp %s:%s/%s/%s.pdb %s' % ( remote_host, labdatabases, subdirectory, pdb_id, dest_dir )

	#lines = popen( wget_cmd ).readlines()
	if ( exists(dest) ):
		return dest
	else:
		print "Error: didn't download file!"

    
#62 63行 将 idealization_app/relax_app修改实际环境中app名称
idealization.mpi.linuxclangrelease
relax.mpi.linuxclangrelease

# 修改 184行，使用relaxed&idealized pdb来进行dssp预测，可能会导致dssp命名错误。不修改的话又会导致CA来源的模型无法正确重构。
# dssp_Dict = parse_dssp_results( options.pdb_fn[:-4]+'_0001_0001.pdb', fasta_dict )

# 在         stderr.write("ERROR:pdb2vall.py: relax app failed at making %s_0001_0001.pdb\n" % pdb ); exit() 这一行下插入:
# IDEALIZE PDB AND GET ITS TORSION ANGLE again. 因为有时候会relax会导致非理想化的结构出现，从而抹掉REMARK的信息，bug。
if exists( pdb + "_0001_0001_0001.pdb"):
    print "pdb2vall.py: detected the relaxed pdb - %s_0001_0001_0001.pdb" % pdb
else:
    print "pdb2vall.py: detected the relaxed pdb - %s_0001_0001.pdb" % pdb
    cmd = idealization_app + " -in:file:s " + pdb + "_0001_0001.pdb -out:file:output_torsions -chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm -chemical:exclude_patches ActeylatedProteinNterm -chainbreaks -overwrite -database %s" % rosetta_database
    system( cmd )
    os.system('mv ' + pdb + '_0001_0001_0001.pdb ' + pdb + '_0001_0001.pdb')


```



2 重写fetch_raw_pdb.py

重写整个脚本。

```python
#!/usr/bin/env python
from sys import argv
from os import system,path
import os
import ConfigParser

labdatabases = os.environ["PDB_DIR"]

# remote host for downloading pdbs
#remote_host = os.environ["INET_HOST"] 处理掉，这个会导致无法正常从pdb.org下载数据

pdbname = argv[1]

netpdbname = labdatabases
netpdbhost = ""
if remote_host:
    netpdbhost = remote_host

if not netpdbname.endswith(path.sep):
    netpdbname += path.sep

netpdbname += pdbname[1:3] + '/' + pdbname

url = 'https://files.rcsb.org/download/%s.pdb' % (pdbname)
wget_cmd = 'wget -c -t 0 %s' % (url)
os.system(wget_cmd)

#if netpdbhost:
#    system("scp %s:%s.pdb ." % ( netpdbhost, netpdbname) )
#else:
#    system("cp %s.pdb ." % netpdbname )
```



3 修改/pdb2vall/jump_over_missing_density.py

这个可能有问题。修改后会导致一部分的pdb无法正确输出vall。

```python
# 修改97行往下的判断条件，可能是因为新版的rosetta输出的pdb header格式改变了，不进行修改的话，python索引切片报错。
                if len(line.split()) > 8:     # modify1
                    idealized_rsd = line.split()[4]
                    print(idealized_rsd)
                    if idealized_rsd != "X":
                        secstr = line.split()[5]
                        phi = line.split()[6]
                        psi = line.split()[7]
                        omega = line.split()[8]

                        torsion_line_edit = "%s %s %s %s %s %s\n" % (
                        newnum_torsion, idealized_rsd, secstr, phi, psi, omega)
                        pdb_torsion_file.append(torsion_line_edit)
                        newnum_torsion = newnum_torsion + 1
```



4 修改make_sequence_fragments.pl

不修改，会卡死在Fragment picker的2%进度。。

```shell
#修改130行为
my $shell = "$PICKER \@picker_sequence_cmd_$id.txt -ignore_unrecognized_res true"; # -j $PICKER_NUM_CPUS >& picker_sequence_cmd_$id.log
```



5 修改/pdb2vall/pdb_scripts/rebuild_pdb_from_CA.sh为一下内容: 因为找不到`runmaxsprout.pl`的绝对路径。

```shell
#!/bin/bash
Bin=$(dirname $(readlink -f $0))
cwd=$(pwd)
pdb=$(basename $1 .pdb )

cp $1 $Bin/maxsprout/$pdb.alpha
cd $Bin/maxsprout/test/
$Bin/maxsprout/test/runmaxsprout.pl ../$pdb.alpha $pdb.rebuilt
cd $cwd
mv  $Bin/maxsprout/test/$pdb.rebuilt ./$pdb.pdb.rebuilt
rm $Bin/maxsprout/$pdb.alpha
rm $Bin/maxsprout/test/$pdb*.log
```



6 修改/structure_profile/make_depthfile.py

目的是避免depth计算时经常出现错误。

```python
# 47行 # GET single chain PDB copy from ../;往后修改为:
# GET single chain PDB copy from ../;
os.system('cp ../../'+pdbname+pdbchain+'.pdb ./')
os.system('mv '+pdbname+pdbchain+'.pdb '+pdbname+'.pdb')

## GET WHOLE CHAIN PDB AND CLEAN THAT UP.
#cmd = PDB2VALL_PATH + 'pdb2vall/pdb_scripts/fetch_raw_pdb.py %s' % pdbname
#system( cmd )

print 'cleaning up pdb'
system( PDB2VALL_PATH + 'pdb2vall/pdb_scripts/clean_pdb.py %s.pdb > /dev/null' % pdbname+pdbchain )
system('rm -f %s.pdb' % argv[1][:4] )
system('ln -s %s.pdb.clean.pdb %s.pdb' % ( pdbname, pdbname ))
system('mv ' + pdbname + '.pdb.clean.pdb ' + pdbname+pdbchain+'.pdb')

wholechain_results  = pdbname + '.depth-residue.depth'
singlechain_results = pdbname+pdbchain + '.depth-residue.depth'

if not exists( wholechain_results ):
    print 'fail to get the depth data from running the whole chains pdb - now try running with pdb with only a chain'

    print 'running DEPTH'
    cmd_depth = DEPTH + ' -thread ' + DEPTH_THREADS + ' -i %s.pdb -o %s.depth' % ( pdbname+pdbchain, pdbname+pdbchain )
    system( cmd_depth )

if exists( singlechain_results ):
    chdir( basedir )
    system("ln -s ./%s/%s ." %( pdbname+pdbchain, singlechain_results ))
else:
    print "job failed - something bad happened here"
```



#### 7 运行

运行一下命令，这样就创建了一个vall数据点了。

注意: 此处pdbid 一定要前4个字母小写，最后一个字幕大写。5code格式，因此大规模计算前一定要先对pdb name进行预处理。

-  --no_structure_profile用于加速计算

```shell
# 以添加2e7dA生成vall为例:
pdbid=2e7dA
mkdir -p ${pdbid} && cd ${pdbid}
python $ROSETTA/tools/fragment_tools/pdb2vall/pdb2vall.py -p ${pdbid} --no_structure_profile
```



最后成功构建后，输出2e7dA.vall，运行时间大致是1小时左右，如果调用44个核心数，处理28000个数据点，大概是二周的时间左右。

```shell
---------------------------------------------------------------------------
pdb2vall.py is trying to combine 2e7dA.vall from:
     1. fasta_dict               (from pdb_seqres.txt)
     2. sectr_dict               (from dssp in rosetta idealizer)
     3. idealized_pdb_xyz_dict
     4. idealized_torsion_dict
     5. dssp_Dict
     6. seq_pro_checkpoint_dict
     7. str_pro_checkpoint_dict / seq_pro_checkpoint_dict
---------------------------------------------------------------------------
pdb2vall.py is done making 2e7dA.vall
```



最后我们通过对非冗余PDB数据列表中的Sequence ID进行循环操作，就可以得到每个pdb的vall数据，将它们合并起来，就得到了自己构建的vall数据库了。(注意，使用pdb2vall处理的vall数据点，是含有CB数据，对于甘氨酸来说，CB坐标默认为9999.999, 并非出错。)

注意: 使用更大的vall数据库，可能存在结构过采样的问题。参考:https://www.rosettacommons.org/node/3335



已知问题:

Bug1: CA_only model虽然修复了，但是依然可能dssp报错导致无法输出vall数据。

Bug2: relax有时候会自动把REMARK的信息抹掉，原因不明。
