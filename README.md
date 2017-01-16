After cloning the repo, start a virtual environment with something like:

```
python3 -m venv ./.env
source ./.env/bin/activate
```

and then:

```
pip install -r requirements.txt
```

to get the required modules and then:

```
jupyter-notebook
```

to open up this directoy in a jupyter notebook page in your browser.

## Work in progress notebook reproducing the work from: https://www.toptal.com/python/comprehensive-introduction-your-genome-scipy
## and then extending the suggested questions at the end of the post


```python
import pandas as pd
```


```python
pd.__version__
```




    '0.19.1'



Download the human genome annotation file in GFF3 format.


```python
!wget ftp://ftp.ensembl.org/pub/release-85/gff3/homo_sapiens/Homo_sapiens.GRCh38.85.gff3.gz
```

    --2017-01-16 14:28:59--  ftp://ftp.ensembl.org/pub/release-85/gff3/homo_sapiens/Homo_sapiens.GRCh38.85.gff3.gz
               => 'Homo_sapiens.GRCh38.85.gff3.gz'
    Resolving ftp.ensembl.org... 193.62.203.85
    Connecting to ftp.ensembl.org|193.62.203.85|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /pub/release-85/gff3/homo_sapiens ... done.
    ==> SIZE Homo_sapiens.GRCh38.85.gff3.gz ... 38469475
    ==> PASV ... done.    ==> RETR Homo_sapiens.GRCh38.85.gff3.gz ... done.
    Length: 38469475 (37M) (unauthoritative)
    
    Homo_sapiens.GRCh38 100%[===================>]  36.69M  2.42MB/s    in 13s     
    
    2017-01-16 14:29:13 (2.78 MB/s) - 'Homo_sapiens.GRCh38.85.gff3.gz' saved [38469475]
    


It is about 37 MB, a very small file compared to the information content of a human genome, which is about 3 GB in plain text. That’s because the GFF3 file only contains the annotation of the sequences, while the sequence data is usually stored in another file format called FASTA.


```python
columnNames = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
```


```python
df = pd.read_csv(
    'Homo_sapiens.GRCh38.85.gff3.gz',
    compression='gzip',
    sep='\t', 
    comment='#', 
    low_memory=False,
    header=None, 
    names=columnNames
)
```


```python
df.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>seqid</th>
      <th>source</th>
      <th>type</th>
      <th>start</th>
      <th>end</th>
      <th>score</th>
      <th>strand</th>
      <th>phase</th>
      <th>attributes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>GRCh38</td>
      <td>chromosome</td>
      <td>1</td>
      <td>248956422</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>ID=chromosome:1;Alias=CM000663.2,chr1,NC_00000...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>.</td>
      <td>biological_region</td>
      <td>10469</td>
      <td>11240</td>
      <td>1.3e+03</td>
      <td>.</td>
      <td>.</td>
      <td>external_name=oe %3D 0.79;logic_name=cpg</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1</td>
      <td>.</td>
      <td>biological_region</td>
      <td>10650</td>
      <td>10657</td>
      <td>0.999</td>
      <td>+</td>
      <td>.</td>
      <td>logic_name=eponine</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1</td>
      <td>.</td>
      <td>biological_region</td>
      <td>10655</td>
      <td>10657</td>
      <td>0.999</td>
      <td>-</td>
      <td>.</td>
      <td>logic_name=eponine</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1</td>
      <td>.</td>
      <td>biological_region</td>
      <td>10678</td>
      <td>10687</td>
      <td>0.999</td>
      <td>+</td>
      <td>.</td>
      <td>logic_name=eponine</td>
    </tr>
  </tbody>
</table>
</div>




```python
df.info()
```

    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 2601849 entries, 0 to 2601848
    Data columns (total 9 columns):
    seqid         object
    source        object
    type          object
    start         int64
    end           int64
    score         object
    strand        object
    phase         object
    attributes    object
    dtypes: int64(2), object(7)
    memory usage: 178.7+ MB



```python
df.seqid.unique()
```




    array(['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
           '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9',
           'GL000008.2', 'GL000009.2', 'GL000194.1', 'GL000195.1',
           'GL000205.2', 'GL000208.1', 'GL000213.1', 'GL000214.1',
           'GL000216.2', 'GL000218.1', 'GL000219.1', 'GL000220.1',
           'GL000221.1', 'GL000224.1', 'GL000225.1', 'GL000226.1',
           'KI270302.1', 'KI270303.1', 'KI270304.1', 'KI270305.1',
           'KI270310.1', 'KI270311.1', 'KI270312.1', 'KI270315.1',
           'KI270316.1', 'KI270317.1', 'KI270320.1', 'KI270322.1',
           'KI270329.1', 'KI270330.1', 'KI270333.1', 'KI270334.1',
           'KI270335.1', 'KI270336.1', 'KI270337.1', 'KI270338.1',
           'KI270340.1', 'KI270362.1', 'KI270363.1', 'KI270364.1',
           'KI270366.1', 'KI270371.1', 'KI270372.1', 'KI270373.1',
           'KI270374.1', 'KI270375.1', 'KI270376.1', 'KI270378.1',
           'KI270379.1', 'KI270381.1', 'KI270382.1', 'KI270383.1',
           'KI270384.1', 'KI270385.1', 'KI270386.1', 'KI270387.1',
           'KI270388.1', 'KI270389.1', 'KI270390.1', 'KI270391.1',
           'KI270392.1', 'KI270393.1', 'KI270394.1', 'KI270395.1',
           'KI270396.1', 'KI270411.1', 'KI270412.1', 'KI270414.1',
           'KI270417.1', 'KI270418.1', 'KI270419.1', 'KI270420.1',
           'KI270422.1', 'KI270423.1', 'KI270424.1', 'KI270425.1',
           'KI270429.1', 'KI270435.1', 'KI270438.1', 'KI270442.1',
           'KI270448.1', 'KI270465.1', 'KI270466.1', 'KI270467.1',
           'KI270468.1', 'KI270507.1', 'KI270508.1', 'KI270509.1',
           'KI270510.1', 'KI270511.1', 'KI270512.1', 'KI270515.1',
           'KI270516.1', 'KI270517.1', 'KI270518.1', 'KI270519.1',
           'KI270521.1', 'KI270522.1', 'KI270528.1', 'KI270529.1',
           'KI270530.1', 'KI270538.1', 'KI270539.1', 'KI270544.1',
           'KI270548.1', 'KI270579.1', 'KI270580.1', 'KI270581.1',
           'KI270582.1', 'KI270583.1', 'KI270584.1', 'KI270587.1',
           'KI270588.1', 'KI270589.1', 'KI270590.1', 'KI270591.1',
           'KI270593.1', 'KI270706.1', 'KI270707.1', 'KI270708.1',
           'KI270709.1', 'KI270710.1', 'KI270711.1', 'KI270712.1',
           'KI270713.1', 'KI270714.1', 'KI270715.1', 'KI270716.1',
           'KI270717.1', 'KI270718.1', 'KI270719.1', 'KI270720.1',
           'KI270721.1', 'KI270722.1', 'KI270723.1', 'KI270724.1',
           'KI270725.1', 'KI270726.1', 'KI270727.1', 'KI270728.1',
           'KI270729.1', 'KI270730.1', 'KI270731.1', 'KI270732.1',
           'KI270733.1', 'KI270734.1', 'KI270735.1', 'KI270736.1',
           'KI270737.1', 'KI270738.1', 'KI270739.1', 'KI270740.1',
           'KI270741.1', 'KI270742.1', 'KI270743.1', 'KI270744.1',
           'KI270745.1', 'KI270746.1', 'KI270747.1', 'KI270748.1',
           'KI270749.1', 'KI270750.1', 'KI270751.1', 'KI270752.1',
           'KI270753.1', 'KI270754.1', 'KI270755.1', 'KI270756.1',
           'KI270757.1', 'MT', 'X', 'Y'], dtype=object)



The seqids starting with KI and GL are DNA sequences – or scaffolds – in the genome that have not been successfully assembled into the chromosomes.

Although the first human genome draft came out more than 15 years ago, the current human genome is still incomplete. The difficulty in assembling these sequences is largely due to complex repetitive regions in the genome.


```python
df.source.value_counts()
```




    havana            1441093
    ensembl_havana     745065
    ensembl            228212
    .                  182510
    mirbase              4701
    GRCh38                194
    insdc                  74
    Name: source, dtype: int64



### How Much of the Genome Is Incomplete?


```python
gdf = df[df.source == 'GRCh38']
```


```python
gdf.seqid.unique().shape
```




    (194,)



### In summary

1. About 0.37% of the human genome is still incomplete even though the first draft came out over 15 years ago.
2. There are about 42,000 genes in the human genome based on this particular GFF3 file we used.
3. The length of a gene can range from a few dozen to over two million bases.
4. Genes are not evenly distributed among the chromosomes. Overall, the larger the chromosome, the more genes it hosts, but for a subset of the chromosomes, the correlation can be negative.
5. The GFF3 file is very rich in annotation information, and we have just scratched the surface. If you are interested in further exploration, here are a few questions you can play with:

### Next Questions

1. How many transcripts does a gene typically have ? 
2. What percentage of genes have more than 1 transcript ?
3. How many isoforms does a transcript typically have ?
4. How many exons, CDS, and UTRs does a transcript typically have? What sizes are they ?
5. Is it possible to categorize the genes based on their function as described in the description column ?

How many transcripts does a gene typically have ?


```python

```

What percentage of genes have more than 1 transcript ?


```python

```

How many isoforms does a transcript typically have ?



How many exons, CDS, and UTRs does a transcript typically have? What sizes are they ?


```python

```

Is it possible to categorize the genes based on their function as described in the description column ?


```python

```
