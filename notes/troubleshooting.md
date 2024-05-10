# A Notebook for current phases of troubleshooting

## 5/10/2024 
Current status of project: Names aligned against WFO and subsequent catalogues (COL, ITIS, GBIF in order of priority). 
Number of Original Names in Ulloa Ulloa et al. 2017: 34538
Number of Names after removal of nonparsable entities + Hybrids: 34375 (-163 Names) 
Number of Names matched to WFO: 30872 Accepted status, 4612 Synonym status 
Number of Names Unmatched to WFO: 1583 Unchecked or NA status
*Note* the sum of the above being greater than the number of original names is due some names capable of mapping to different taxonomicStatuses depending on the underlying Authorship
Number of instances of names containing multiple mappings: 3923 instances
Number of Names affected by multiple mapping in taxonomicStatus/nameAlignment: 1679 Names
Number of Names considered 'recoverable' given authorship present in the data: 1451 Names  
Number of Names 'unrecoverable' even with authorship present in the data: -228 Names 

Overall Summary for WFO names
Original Names: 34538 
Names Matched: 33222
Names Unmatched: 925 
Names lost due to unparsable/hybrid + unresolvableMultipleMappings: 391
33222 + 925 + 391 == 34538 TRUE -check- 

Overall Summary of Alternative Taxonomic Backbone Matches
925 names are matched against the alternative catalogues/backbones and coalesced in priority of COL, ITIS, then GBIF. 
COL names matched: 709
ITIS names matched: 50
GBIF names matched: 60
Remaining names unmatched: 106 names 

Current Status: On line 249: Checking multipleAcceptableUsageIDs for multiple paths. A bit complicated...for some reason im getting that ALL synonyms are of that status (they should map to a single on with a sgrouped argument...check back in on this next.

To do on the Horizon: 
- Need to attach the synonyms to each resolution: wfo, col, itis, and gbif. 
- Need to make a dataframe that takes all of these and builds a relation table
- Need to double check that no synonyms double up across the catalogues (or even within catalogs after looking at the organization...)
- Need to attach Native/Non-Native statuses for each species and its corresponding names. Might be a bit challenging as this requires either A) Aligning with WCVP's taxonomy then checking N/NonN status, or B) ... unsure actually
- Need to create a count loop that sums up the number of records expected for each species. This is complicated given the code currently. Try gbifdb IF possible...  
