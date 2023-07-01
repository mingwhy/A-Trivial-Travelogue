https://bioinf.shenwei.me/taxonkit/chinese/

*install*
conda install taxonkit -c bioconda -y
# 表格数据处理，推荐使用 csvtk 更高效
conda install csvtk -c bioconda -y

*download NCBI Taxonomy data*
# 有时下载失败，可多试几次；或尝试浏览器下载此链接
wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -zxvf taxdump.tar.gz

# 解压文件存于家目录中.taxonkit/，程序默认数据库默认目录
mkdir -p $HOME/.taxonkit
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit

*list 列出指定TaxId所在子树的所有TaxID*
# 以人属(9605)和肠道中著名的Akk菌属(239934)为例
$ taxonkit list --show-rank --show-name --indent "    " --ids 9605,239934

*lineage 根据TaxID获取完整谱系*
$ head taxids.txt #create taxids.txt first
9606
9913
376619
# 查找指定taxids列表的物种信息，tee可输出屏幕并写入文件
$ taxonkit lineage taxids.txt | tee lineage.txt 

*name2taxid 将分类单元名称转化为TaxID*
将分类单元名称转化为TaxID非常容易理解，唯一要注意的是某些TaxId对应相同的名称，比如
# -i指定列，-r显示级别，-L不显示世系
$ echo Drosophila | taxonkit name2taxid | taxonkit lineage -i 2 -r -L

*lca 计算最低公共祖先(LCA)*
TaxID的分隔符可用-s/--separater指定，默认为" "。
# 计算两个物种的最近共同祖先，以上面尼安德特人亚种和海德堡人种
$ echo 63221 2665953 | taxonkit lca


*实战，LCA of human, mouse, rat* 
$ echo Homo sapiens | taxonkit name2taxid | taxonkit lineage -i 2 -r -L
Homo sapiens	9606	species
$ echo Mus Musculus | taxonkit name2taxid | taxonkit lineage -i 2 -r -L
Mus Musculus	10090	species
$ echo Rattus norvegicus | taxonkit name2taxid | taxonkit lineage -i 2 -r -L
Rattus norvegicus	10116	species

$ echo 9606 10090 10116 | taxonkit lca
9606 10090 10116	314146

$echo 314146 | taxonkit lineage | taxonkit reformat | cut -f 1,3
314146	Eukaryota;Chordata;Mammalia;;;;

