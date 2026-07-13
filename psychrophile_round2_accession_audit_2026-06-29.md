# Psychrophile Round 2 Accession Audit

Date checked: 2026-06-29

This note records the exact evidence directly verified for the 12 psychrophile
rows added in the latest expansion round.

## Evidence standard used

- Psychrophile call: local TEMPURA snapshot at
  `thesis_data_collection/200617_TEMPURA.csv`
- Sequence mapping: live NCBI Assembly E-utilities checks
- Acceptance rule for this round: strict psychrophile `Topt_ave <= 15` in the
  local TEMPURA snapshot plus an exact-strain or exact type-material assembly
  match in NCBI Assembly metadata

## Verified rows

| Taxon | Curated strain | TEMPURA `Topt_ave` | Accepted accession | Direct NCBI Assembly match basis |
| --- | --- | ---: | --- | --- |
| Shewanella atlantica | HAW-EB5 | 10.0 | `GCF_003966265.1` | Assembly summary says `assembly from type material`; biosource lists `strain=HAW-EB5`, `culture-collection=MCCC:1A12968` |
| Shewanella canadensis | HAW-EB2 | 10.0 | `GCF_003966225.1` | Assembly summary says `assembly from type material`; biosource lists `strain=HAW-EB2`, `culture-collection=MCCC:1A12967` |
| Lacinutrix jangbogonensis | PAMC 27137 | 10.0 | `GCF_000797445.1` | Assembly summary says `assembly from type material`; biosource lists `strain=PAMC 27137` |
| Cryobacterium aureum | Hh31 | 11.0 | `GCF_002954245.1` | Assembly summary says `assembly from type material`; biosource lists `strain=Hh31`, `culture-collection=CGMCC:1.11213` |
| Cryobacterium levicorallinum | Hh34 | 11.0 | `GCF_004402405.1` | Assembly summary says `assembly from type material`; biosource lists `strain=Hh34`, `culture-collection=CGMCC:1.11211` |
| Photobacterium carnosum | TMW 2.2021 | 12.5 | `GCF_002849605.1` | Assembly summary says `assembly from type material`; biosource lists `strain=TMW 2.2021` |
| Gillisia hiemivivida | IC154 | 13.0 | `GCF_007997295.1` | Assembly summary says `assembly from type material`; biosource lists `strain=IC154`, `culture-collection=ACAM:1061` |
| Sulfurimonas gotlandica | GD1 | 15.0 | `GCF_000242915.1` | Assembly summary says `assembly from type material`; biosource lists `strain=GD1` |
| Porticoccus hydrocarbonoclasticus | MCTG13d | 15.0 | `GCF_000744735.1` | Assembly summary says `assembly from type material`; biosource lists `strain=MCTG13d` |
| Neptunomonas antarctica | S3-22 | 15.0 | `GCF_001305295.1` | Assembly summary says `assembly from type material`; biosource lists `strain=S3-22` |
| Lentisphaera profundi | SAORIC-696 | 15.0 | `GCF_028728065.1` | Assembly summary says `assembly from type material`; biosource lists `strain=SAORIC-696` |
| Arsukibacterium ikkense | GCM72 | 15.0 | `GCF_000987775.1` | Assembly summary says `assembly from type material`; biosource lists `strain=GCM72` |

## Notes

- I did not add primary-paper claims for these 12 rows in this cleanup pass.
- I only retained claims re-verified from the local TEMPURA snapshot and the
  live NCBI Assembly summaries.
- This file is meant to document the provenance of the latest psychrophile
  additions without over-claiming beyond what was checked.
