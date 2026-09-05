# Manifest

`manifest_reference.tsv` is the source manifest for the complete six-annotator LRGASP SQANTI3 basic-statistics workflow.

It contains 90 rows representing:

- 6 annotators: LoRTIA-Plus, bambu, FLAIR, IsoQuant, NAGATA, StringTie2
- 5 chemistries: ONT-CapTrap, ONT-cDNA, ONT-dRNA, PacBio, PacBio-CapTrap
- 3 cell lines: H1, H1-endo, WTC11

Columns:

```text
Chemistry    Cell-line    GTF
```

Important legacy note: despite its name, the `GTF` column stores paths to SQANTI3 `classification.txt` files in this reused manifest.

The archived manifest preserves the original absolute paths used in the analysis. These paths are machine-specific. For a full rerun on another system, edit or regenerate a working copy of the manifest so that the `GTF` entries point to the corresponding SQANTI3 classification files on that system. Do not modify the archived manifest if exact provenance is required.
