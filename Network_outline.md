## ATAC	| RNA-Seq Network Construction Proposal

1. Align 
2. Peak call (homer narrow?)
3. subtract likely TSS/promoter regions
4. ID motifs
    	- YAMDA for novel
    	- BAO (unavailable?) or CENTIPEDE for known?
5. motif <	-> TF 
  	- unsure how...HOMER does this built in
6. Assign motifs to gene
  	- closest 1 (or 2) within 1mb
    	- need to do some reading on enhancer / gene distances
    	- 1 mb is often used, but I'd *like* to get something a bit more rigorous
    	- Quote from someone at BoG18:
      	- "just use closest gene, otherwise you'll never get anything done""
  	- perhaps overlay TAD?
    	- would be useful to parition with that
    	- even if other tissue
      	- often claimed that TAD are 'highly correlated' between tissues
  	- quantify gene expression 
    	- or group low / med / high?
    	- on / off?
    	- will be used below
7. Rank motifs by gene expression
  	- sum?
  	- average?
  	- variance?
  	- rank/weight vertice by strength of gene expression?
    	- or by change from expected?
      	- hmm, this could be directed then....
      	- I do have three states (iPSC, GFP+, and RFP+)
8. Retain top n motifs
  	- density plot to assess
  	- ideally would be bimodal
  	- un-ideal would be huge tail
    	- almost certainly will be this
9. Collapse motif to likely TF?
  	- or assign top n TF to motif
  	- or p val cutoff?
10. Draw (un)directed vertices from gene to gene
  	- one vertice per motif
  	- so some gene <-> gene are multicolored for each motif / TF
  	- if I use delta expression, then this could be directed
    	- motif/TF associated with increase gene expression are activating...
    	- or repressive