
#######################################################################################
<h1>PLEASE NOTE</h1>
<p>ExInAtor is ONLY designed to be used with mutations from Whole Genome Sequences.
 It is NOT designed to be used with Exome data.</p>
#######################################################################################
<h1>Table Of Contents
<ul>
<li><a class="reference internal" href="#authors">Authors</a></li>
<li><a class="reference internal" href="#description">Description</a></li>
<li><a class="reference internal" href="#requirements">Requirements</a></li>
<li><a class="reference internal" href="#installation">Installation</a></li>
<li><a class="reference internal" href="#inputs">Inputs</a></li>
<li><a class="reference internal" href="#outputs">Outputs</a></li>
<li><a class="reference internal" href="#qqplots">QQplots</a></li>
</ul>

<h1>Authors<a class="headerlink" href="#authors" title="Permalink to this headline">¶</a></h1>
<p><a class="reference external" href="mailto:andres.lanzos@dbmr.unibe.ch">Andrés Lanzós</a> <sup>1,2,3</sup>
<a class="reference external" href="mailto:joana.carlevaro@dbmr.unibe.ch">Joana Carlevaro-Fita</a> <sup>1,2,3</sup>
<a class="reference external" href="mailto:loris.mularoni@irbbarcelona.org">Loris Mularoni</a> <sup>3,4</sup>
<a class="reference external" href="mailto:freverter@ub.edu">Ferran Reverter</a> <sup>1,2,3</sup>
<a class="reference external" href="mailto:emilio.palumbo@crg.eu">Emilio Palumbo</a> <sup>1,2,3</sup>
<a class="reference external" href="mailto:roderic.guigo@crg.eu">Roderic Guigó</a> <sup>1,2,3</sup>
<a class="reference external" href="mailto:rory.johnson@dbmr.unibe.ch">Rory Johnson</a> <sup>1,2,3,5 *</sup></p>
<ol class="arabic simple">
<li>Centre for Genomic Regulation (CRG), The Barcelona Institute of Science and Technology, Dr. Aiguader 88, Barcelona, 08003, Spain.</li>
<li>Universitat Pompeu Fabra (UPF), Barcelona, Spain.</li>
<li>Institut Hospital del Mar d’Investigacions Mèdiques (IMIM), 08003 Barcelona, Spain.</li>
<li>Research Unit on Biomedical Informatics, Department of Experimental and Health Sciences, Universitat Pompeu Fabra, Dr. Aiguader 88, Barcelona, Spain.</li>
<li>Present address: Department of Clinical Research, University of Bern, Murtenstrasse 35, 3010 Bern, Switzerland.</li>
</ol>
<ul class="simple">
<p>* Corresponding author.</p>
</ul>
</div>
<div class="section" id="description">
<h1>Description<a class="headerlink" href="#description" title="Permalink to this headline">¶</a></h1>
<p>ExInAtor is designed to detect cancer driver long non-coding RNA genes.
The general approach is to identify genes that have an excess of mutations in their exons compared to local background regions (including intronic sequences)..
The statistic method is based on
the Hypergeometric distribution and generating a Q (False Discovery Rate) value for each gene.
The main output file is called &#8220;ExInAtor_Gene_List.txt&#8221;, which contains
comprehensive results for all analyzed genes.</p>
<p>The general approach is represented in the following figure:</p>
<p><strong>A.</strong> Example gene to analyse and two flanking genes.</p>
<p><strong>B.</strong> Merging of exons for each gene.</p>
<p><strong>C.</strong> Extension of the background region while removing any exon of any flanking gene.</p>
<p><strong>D.</strong> Mutation mapping.</p>
<p><strong>E.</strong> Definition of the new background region to obtain the same trinucleotide content (i.e. percentage) than the exonic region.</p>
<p><strong>F.</strong> Mutation mapping on the new background region.</p>
<a class="reference internal image-reference" href="workflow"><img alt="workflow" class="align-center" src="https://github.com/alanzos/ExInAtor/blob/master/exinator_detailed_workflow.jpg" style="width: 1123.0px; height: 794.0px;" /></a>
</div>
<div class="section" id="requirements">
<h1>Requirements<a class="headerlink" href="#requirements" title="Permalink to this headline">¶</a></h1>
<p><strong>1. Linux</strong>:</p>
<p>Tested in Ubuntu 14.04</p>
<p><strong>2. Bedtools</strong>:</p>
<p>Please use version 2.19.1. More recent versions will not work because of changes in some bedtools commands.</p>
<p><strong>3. Python</strong>:</p>
<p>Tested on version 2.7.6</p>
<p><strong>4. R</strong>:</p>
<p>Tested on versions 3.3.1</p>
<p><strong>5. Awk</strong>:</p>
<p>Tested on versions mawk 1.3.3</p>
</div>
<div class="section" id="installation">
<h1>Installation<a class="headerlink" href="#installation" title="Permalink to this headline">¶</a></h1>
<p><strong>1. Download the compressed file &#8220;All_files.zip&#8221; from:</strong></p>
<p><a class="reference external" href="https://github.com/alanzos/ExInAtor">https://github.com/alanzos/ExInAtor</a></p>
<p><strong>2. Uncompress it in path</strong></p>
<p><strong>3. Ready to run</strong></p>
<p>For each file you must provide the full path. Example command:</p>
<div class="highlight-Description"><div class="highlight"><pre>$  python2.7 Main.py -i /home/All_files/Inputs/Breast.bed -o /home/All_files/Outputs -f /home/All_files/Inputs/Genome_v19.fasta -g /home/All_files/Inputs/gencode.v19.long_noncoding_RNAs.gtf -s /home/All_files/Inputs/chromosomes.txt -k /home/All_files/Inputs/3mers.txt -w /home/All_files/Inputs/gencode.v19.annotation.gtf -c 6 -n 119 -b 10000
</pre></div>
</div>
<p>In the folder &#8220;Outputs&#8221; you can find the expected results.</p>
</div>
<div class="section" id="inputs">
<h1>Inputs<a class="headerlink" href="#inputs" title="Permalink to this headline">¶</a></h1>
<p>All of the following files are provided in the compressed file &#8220;All_files.zip&#8221; except the Fasta file of the whole genome (&#8220;-f | &#8211;fasta_file&#8221;) which can be obtained <a class="reference external" href="https://www.dropbox.com/s/a6vthezotm6iaih/Genome_v19.fasta.gz?dl=0">here (818 MB)</a></p>
<dl class="docutils">
<dt><strong>1. Mandatory</strong>:</dt>
<dd><ul class="first last simple">
<li>-i | &#8211;input_file -&gt; File containing the localization of the cancer mutations in BED format.</li>
<li>-o | &#8211;output_folder -&gt; Name of the folder where files and plots are saved.</li>
<li>-g | &#8211;gtf_file -&gt; File containing the genes and exons to analyse in GTF format. Other features, including transcripts, will be ignored.</li>
<li>-f | &#8211;fasta_file -&gt; Fasta of the whole genome.</li>
<li>-s | &#8211;chr_sizes -&gt; Two-column tab-separated text file containing assembly sequence names and sizes.</li>
<li>-k | &#8211;kmers_file -&gt; Txt file containing all the possible trinucleotides.</li>
<li>-w | &#8211;whole_genome -&gt; File containing the genes and exons of all the genome in GTF format. Other features, including transcripts, will be ignored.</li>
<li>-n | &#8211;number_of_genomes -&gt; The number of samples or genomes corresponding to the &#8220;input_file&#8221;.</li>
</ul>
</dd>
<dt><strong>2. Optional</strong>:</dt>
<dd><ul class="first last simple">
<li>-e | &#8211;exonic_filter -&gt; Minimun number of exonic mutations a gene must have to be analyzed.</li>
<li>-x | &#8211;background_filter -&gt; Minimun number of background mutations a gene must have to be analyzed.</li>
<li>-b | &#8211;background_size -&gt; the extension length of the background region that includes all introns.</li>
<li>-c | &#8211;cores -&gt; the number of CPU cores to use in the analysis.</li>
</ul>
</dd>
</dl>
</div>
<div class="section" id="outputs">
<h1>Outputs<a class="headerlink" href="#outputs" title="Permalink to this headline">¶</a></h1>
<p><strong>1. ExInAtor_Gene_List.txt</strong>: Final output. List of driver candidates with exon_mutations, exon_length, intron_mutations, intron_length, pval, qval. Sorted by Qval.</p>
<p><strong>2. genes.bed</strong>: BED file containing the localization of each gene.</p>
<p><strong>3. counts.txt</strong>: File containing the number of exonic mutations, exonic length, background mutations and background length for each gene.</p>
<p><strong>4. table_kmer_counts.txt</strong>: File containing the number of exonic mutations, exonic length, background mutations and background length for each trinucleotide in each gene.</p>
<p><strong>5. qqplot.png</strong>: QQplot of the expected and observed pvalues</p>
</div>
<div class="section" id="qqplots">
<h1>QQplots<a class="headerlink" href="#qqplots" title="Permalink to this headline">¶</a></h1>
<p>These plots are used to evaluate whether the Pvalues follow a Uniform distribution, and hence the overall statistical behaviour of the results.
&#8220;QQ&#8221; stands for Quantile-Quantile plot. The point of these figures is to compare the distribution of observed pvalues to the expected distribution under the null hypothesis of no association.
The null hypothesis in this case would generate a uniform distribution, this means, a flat histogram over all statistical tests with a total density of 1.</p>
<p>For more information see:</p>
<p><a class="reference external" href="https://en.wikipedia.org/wiki/Q%E2%80%93Q_plot">https://en.wikipedia.org/wiki/Q%E2%80%93Q_plot</a></p>
<p><a class="reference external" href="http://www.jstor.org/stable/2987987?origin=crossref&amp;seq=1#page_scan_tab_contents">http://www.jstor.org/stable/2987987?origin=crossref&amp;seq=1#page_scan_tab_contents</a></p>
<p>The following is the expected QQplot of the example data provided with ExInAtor in the GitHub link: <a class="reference external" href="https://github.com/alanzos/ExInAtor">https://github.com/alanzos/ExInAtor</a></p>
<a class="reference internal image-reference" href="qqplot"><img alt="qqplot" class="align-center" src="https://github.com/alanzos/ExInAtor/blob/master/QQplot.png" style="width: 800.0px; height: 800.0px;" /></a>
