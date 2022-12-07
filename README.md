# SVGMM
本程序为无监督PolSAR图像分类方法
在测试你自己的数据时，你需要做以下步骤：
1）下载VLfeat9_20库，并将其放在"VLfeat9_20"文件夹下
2）下载数据，如欧空局网站：https://earth.esa.int/web/polsarpro/airborne-data-sources
              或AFS网站：https://uavsar.jpl.nasa.gov/cgi-bin/data.pl
3）使用PolSAR pro软件对数据进行Lee滤波预处理，将数据保存为C3形式的mat文件，然后存入data文件夹
   注意：.mat数据中数据X的格式为：[Sz(1)*Sz(2),9]，其中，Sz(1)和Sz(2)分别为图像尺寸
4)程序输出为分类标签及分类结果图，这2个结果将保存在Result文件夹下的文件夹中，文件夹与数据集名称请保持相同
