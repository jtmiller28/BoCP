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

Current Status: On line 249: Checking multipleAcceptableUsageIDs for multiple paths. A bit complicated...for some reason im getting that ALL synonyms are of that status (they should map to a single on with a sgrouped argument...check back in on this next. <br>

To do on the Horizon: 
- Need to attach the synonyms to each resolution: wfo, col, itis, and gbif. 
- Need to make a dataframe that takes all of these and builds a relation table
- Need to double check that no synonyms double up across the catalogues (or even within catalogs after looking at the organization...)
- Need to attach Native/Non-Native statuses for each species and its corresponding names. Might be a bit challenging as this requires either A) Aligning with WCVP's taxonomy then checking N/NonN status, or B) ... unsure actually
- Need to create a count loop that sums up the number of records expected for each species. This is complicated given the code currently. Try gbifdb IF possible...  