# Real-RNA validation corpus

`rfam_real_24.json` contains one biological sequence and the projected consensus secondary structure from each of 24 diverse Rfam SEED families. Rfam SEED alignments are manually curated representative alignments whose structure annotations are sourced from publications, specialist databases, or expert curation. Twenty selected families report a published structure source; four report a curated prediction source.

- Source: <https://rfam.org/>
- API: `https://rfam.org/family/{accession}/alignment`
- License: CC0 1.0
- Rebuild: `python3 scripts/fetch-rfam-corpus.py`

The rebuild script chooses the 20-500 nt IUPAC sequence with the smallest gap fraction in each SEED alignment, removes alignment gaps, and projects `SS_cons` pairing onto that sequence. The accession, sequence region identifier, description, structure provenance, retrieval date, and selection rule remain in the JSON file.

# Pseudoknot corpus

`pseudoknot_real_24.json` contains the first 24 records of the published Andronescu–Pop–Condon S-Test short-pseudoknot set. The upstream text identifies RNA STRAND/PseudoBase/literature provenance per record. Sequences are normalized to upper-case RNA (`T` is accepted and normalized by Ribon); reference structures retain their published crossing bracket topology.

- Source: <https://www.cs.ubc.ca/labs/algorithms/Publications/PaperMaterials/PseudoRNA/data/TES_pkshort_max200_RenJab-S20.txt>
- Selection: first 24 non-comment records, order preserved
- Runtime use: validation only; not bundled into the Typst package
