# A Notebook for current phases of troubleshooting

## 5/10/2024 
### Current status of project: Names aligned against WFO and subsequent catalogues (COL, ITIS, GBIF in order of priority). <br> 
Number of Original Names in Ulloa Ulloa et al. 2017: 34538 <br>
Number of Names after removal of nonparsable entities + Hybrids: 34375 (-163 Names) <br> 
Number of Names matched to WFO: 30872 Accepted status, 4612 Synonym status <br>
Number of Names Unmatched to WFO: 1583 Unchecked or NA status <br>
*Note* the sum of the above being greater than the number of original names is due some names capable of mapping to different taxonomicStatuses depending on the underlying Authorship <br>
Number of instances of names containing multiple mappings: 3923 instances <br>
Number of Names affected by multiple mapping in taxonomicStatus/nameAlignment: 1679 Names <br>
Number of Names considered 'recoverable' given authorship present in the data: 1451 Names  <br>
Number of Names 'unrecoverable' even with authorship present in the data: -228 Names <br>

### Overall Summary for WFO names 
Original Names: 34538 <br> 
Names Matched: 33222 <br>
Names Unmatched: 925 <br>
Names lost due to unparsable/hybrid + unresolvableMultipleMappings: 391 <br>
33222 + 925 + 391 == 34538 TRUE -check- <br>

### Overall Summary of Alternative Taxonomic Backbone Matches
925 names are matched against the alternative catalogues/backbones and coalesced in priority of COL, ITIS, then GBIF. <br> 
COL names matched: 709 <br>
ITIS names matched: 50 <br>
GBIF names matched: 60 <br>
Remaining names unmatched: 106 names <br> 

Alignments were merged across all backbones, and assessed for duplications of acceptedNames across backbones and duplication of synoynms across backbones. These names were removed via priority of wfo > col > itis > gbif <br> 

Synonyms were then assessed for multipleMapping status in the WFO backbone, as the other backbones were judged to add little data in comparison. for those names that will map to multiple accepted names per just the genus+specificEpithet of the synoynm, I created another field called synonymWFOSpacelessAuthorship. During data prep, much like how we will require that multipleMapping AcceptedNames have an associated valid authorship, we will require this for these synonyms as well. <br> 



### Final Summary 
- # of Accepted Names: 33126
- # of Synoynms attached to these Accepted Names: 125175
- # of Synonyms with multiple mappings possible (WFO only, aka these require authorship): 8960
To do on the Horizon: 
- ~~Need to attach the synonyms to each resolution: wfo, col, itis, and gbif.~~ 
- ~~Need to make a dataframe that takes all of these and builds a relation table~~
- ~~Need to double check that no synonyms double up across the catalogues (or even within catalogs after looking at the organization...)~~
- Need to attach Native/Non-Native statuses for each species and its corresponding names. Might be a bit challenging as this requires either A) Aligning with WCVP's taxonomy then checking N/NonN status, or B) ... unsure actually
- Need to create a count loop that sums up the number of records expected for each species. This is complicated given the code currently. Try gbifdb IF possible...


### 06/19/2024 CHANGES
We are adapting wcvp instead of WFO in order to accomidate phylogenetic issues. Rewriting, but noting things I have quesstions about during the adaptation. 
- Biotypes, not exactly sure whether these should or shouldnt be applicable for harmonization. Keeping for now but check with other's about this designation status
- Seems that the cases of multiple mappings within WCVP that I've checked are due to subspecific levels of desigation, rather than fully alternative names. 
- Appears that the class illegitimate for the field taxon_status does not actually mean its unusable, it just means that name should be classified alternatively? 
Current Report for wcvp alignment with North american name list 
Accepted: 30179 Names
Illegitimate: 117 Names
Invalid: 5 Names
Misapplied: 4 Names
Orthographic: 4 Names
Synonym: 3573 Names
Unplaced: 66 Names
NA: 355 Names 
Of these NA names, coalescent resolution harmonized 252 of these names, leaving a remaining 103 unresolved. All togetherthat is 169 names without placement from our NA list. 

Finalized Accepted Name Count: 33239 Names, Finalized Synonym Count 120,744

06/21/2024
Edge case issue with alignment: 
Abies lasiocarpa var. arizonica is being mapped to both Abies grandis and Abies lasiocarpa with identical synonym authorship. Abies grandis shouldnt be mapped to at all, so something is problematic in this regard. 
 
