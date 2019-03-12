Two diff processes:

```
ipsc ------> RFP --> GFP


		or 

	 RFP
	  ^
	  |
ipsc -------> GFP
```

Data we have:

	- iPSC / GFP / RFP

		- RNA-seq (bulk)
		- scRNA-seq for iPSC and GFP/RFP
		- ATAC-seq
	
External:

	- H3K27Ac for RPE

Can we prove a model? I don't think so...could try trajectory analysis with the scRNA-seq data?

	- trajectory analysis *may* show RFP -> GFP trajectory branching...unlikely 

Rough pipeline is:

	- call peaks, intersect RPE peaks with h3k27ac
	- GFP unique peaks (against iPSC) -> Homer TFBS enrichment
	- RFP unique peaks (against iPSC) -> Homer TFBS enrichment
	- GFP unique peaks (against RFP) -> Homer TFBS enrichment

Take TFBS lists

	- two types of enriched TFBS 
		- no differential gene expression
			- presumably some **other** TF / process is altering chromatin accessability, which allows TF to alter activity
			- in other words, the TF is downstream of something else
		- differential gene expression
			- more likely to be a driver, as expression and # of chromatin sites are changed
	
Take top* TF, find set of genes with:

	- delta {comparison} (GFP not iPSC, GFP not RFP, etc) gene expression 
	- **and** 
	- matched TF motif in a {comparison} unique peak

* where top is TFBS with < -100 ln(p-value). Also consider whether they TF has delta gene expression.

Analysis:

	- Take gene lists and iteratively compare to find TF with overlapping sets of genes
		- this will merge (semi?)-redundant TF(BS) 

	- From the overlapping gene sets, do GO term enrichment

	- Use the TF as potential targets for *in vitro* perturbation
		- nature of TF and GO term enrichment of potential targets to order the TF list
		- **nature** being literature search, known biology, etc.
		- also give higher priority to TF with differentiatial expression

Do Analysis in two ways:

	- GFP and RFP shared TF / diff expressed genes compared to iPSC
		- this will find processes from iPSC --> RPE like (the assumption is that RFP --> GFP is simply more time / more of the same factors)
	- GFP compared to RFP
		- this assumes that something has changed (more chromatin site / diff gene expression) to give GFP (well high TYR) state compared to RFP

