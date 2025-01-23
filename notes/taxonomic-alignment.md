### Taxonomic Alignment notes 
Goal: Combine NCBI taxonomy with WCVP's alignment in order to attach DNA sequences (where avaliable) to the names. 

To reduce the complexity of the NCBI taxonomic backbone the following filters were applied:
* Removal of names containing " cf " as this designates uncertainty. 
* Removal of names containing " x " as these are hybrids which we're electing to not deal with
* Removal of names containing " sect " as these are not dealing with species level alignment 
* Removal of names containing "/" as these express uncertainty

Relational tables were built per taxonomic database following these steps: <br>
NCBI
1) Names were parsed into name and author, names that failed to parse or failed to parse past the genus level were removed.
2) Names with the taxonomic status "in-part" were removed due to their flexible mapping nature (multiple alignments by definition.)
3) Names were checked to see if they aligned with unique accepted names, those that cause did not were removed after determining they contributed to a total of 0.02% of the names available for matching. 
4) The following variations of taxonomic status were checked, and changed to "synonym": "common name", "equivalent name", "includes", "type material", "genbank common name".
5) The taxonomic status "authority" was a desigation that contained authorship within the name, these were changed to either "accepted with authorship" or "synonym with authorship" based on whether the provided name matched the ncbi aligned name, and whether the field ncbiAuthor had data from parsing.
6) Redundant names of statuses were removed from implementing step (5)
7) Names were reorganized to add name field (any possible name in the database) and its aligned name

WCVP
1) Names were reorganized to build a table with all possible names in the database and their accepted name counterparts. 
2) Unplaced names were removed entirely as they cannot arrive at a correct alignment

Names were then checked to see if all names are unique in what aligned name they arrive at, and marked by the boolean field multipleNameAlignmentsPossible TRUE where a name is NOT unique. After testing on these databases, it was determined that only ~2% of the total names contained these non-unique mapping properties in wcvp, and only 0.02% of the names in the filtered ncbi database. These names were removed, as the complexity associated with correctly mapping them requires the standardization of authorship, an often missing piece in most occurrence data. One exception was made in the removal of non-unique mapping names, in the case that one version of a name that was non-unique in its mapping was an accepted name by wcvp alignment it was retained. This decision was made based on my experience that these names are much more common than their alternative counterparts and it would have ill effects on the downstream product to remove an accepted name from the list, but should be noted that this is one source of error.

Combining WCVP and NCBI
wcvpAligned accepted name was connected to all versions of the wcvp and ncbi available names. the wcvpAligned names were also simplified into wcvpAlignedNameParent, containing only the first two word's in the character string of wcvpAlignedName in order to remove infraspecific information. These combined names were then assessed for multiple ncbi ids,by grouping the wcvpAlignedNmae and looking at whether there were unique ncbiAcceptedNameUsageIDs per wcvpAlginedName grouping. All boolean values were false, therefore I backfilled ncbiAcceptedNameUsageIDs to all versions of the wcvpAlignedName. The output was then simplified by reducing fields to the necessary ones for end interpretation and coalesced for distinct names based on priority of wcvp, then ncbi. 
 
Field descriptions:
wcvpAlignedGenus = the genus from the aligned name according to wcvp
wcvpAlignedNameParent = the first 2 words of the wcvpAlignedName
wcvpAlignedName = the aligned name according to wcvp
wcvpAlginedNameAuthors = the aligned name authors according to wcvp
wcvpAlignedNameStatus = the aligned name status (should always be "Accepted")
name = any name found in either the wcvp or ncbi backbones.
wcvpNameAuthors = name's authors (if present, wcvp only)
nameStatus = the designated status for any name: can be "Synonym"
"Accepted" etc.
source = where the name came from (either "ncbi" or "wcvp")
ncbiAcceptedNameUsageID = the NCBI ID that will coordinate with 
phylogenies. All names attached to an aligned name will have these
backfilled
multiple_ncbi_ids = a boolean check to see if there are any conflicts
in NCBI ID with the wcvp's backbone.
