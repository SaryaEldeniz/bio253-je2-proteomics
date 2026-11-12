# Step 3: identifier map and annotation

Generated: 2026-11-12 13:24:39 UTC

Step 2 was skipped at the user's request because the data were already validated.

Master map: 4189 entries; 55076 searchable alias records; 5/5 mapping checks passed.

## Coverage of the current DEA feature universe

Phospho counts below are unique parent proteins, not phosphosite counts. Both universes use the current ancestor TSB-vs-SN contrast.

| Layer | Tested parent entries | CDS | Good JE2 maps | Unique UniProt | CDS STRING aliases | CDS with terms | CDS with physical edges |
|---|---:|---:|---:|---:|---:|---:|---:|
| total | 2324 | 2229 | 2107 | 2080 | 2229 | 2143 | 1640 |
| phospho | 1146 | 1109 | 1068 | 1063 | 1109 | 1097 |  937 |

STRING alias resolution, functional-term coverage and physical-edge coverage are separate quantities.

## Classification correction

Metadata identifies 68 total-DEA contaminant entries not flagged by its is_contaminant column.
The primary ordinary-protein analysis uses entry_class == cds. Unmatched, weak and excluded rows remain in the local audit.

## Example public identifiers

These are user-editable starting candidates, not newly selected significant hits.

| Label | JE2 locus | UniProt reference ortholog | Mapping quality |
|---|---|---|---|
| Stk1 | B7H15_06355 | Q2FZ64 | good |
| Stp1 | B7H15_06350 | Q2FZ65 | good |
| TcyA | B7H15_13375 | Q2FVL2 | good |
| TcyB | B7H15_13370 | Q2FVL3 | good |
| TcyC | B7H15_13365 | Q2FVL4 | good |
| TcyP-like | B7H15_02105 | Q2G0Z4 | good |
| AhpC | B7H15_02095 | P0A0B7 | good |
| KatA | B7H15_07015 | Q2FYU7 | good |
| NorA | B7H15_03880 | Q2G0A2 | good |

## Use and limits

UniProt assignments are reference-proteome ortholog mappings to S. aureus NCTC8325. Good means E-value <= 1e-5 and length ratio 0.9-1.1.
Keep JE2 sequence coordinates for phosphosites. Preserve many-to-one mappings; ORA uses unique target IDs and ranked analysis needs an explicit collision rule.
Local physical-network partners are associations from the supplied export, not evidence of direct kinase substrates.

See docs/step03_identifiers.md and docs/web_resources.md for interactive R use and official API references.

## Input provenance

- Mapping zip MD5: `049dcd5bb6ac0650879099cbaaf8576e`
- Local STRING aliases MD5: `3958b07986103fe4c060790b2116805d`
